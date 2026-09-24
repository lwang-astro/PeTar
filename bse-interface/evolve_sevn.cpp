//
// Created by Daniel Marín Pina
//
#include <csignal>
#include <cmath>
#include <IO.h>
#include <string>
#include <vector>

#include "sevn.h"
#include "bse_interface.h"

constexpr double VERY_LARGE_TIME{1e50};
constexpr double year_to_day(365.25);

std::map<std::string, std::string> input_params_SEVN;
std::map<std::string, std::pair<std::string, std::string>> default_params_sevn;

void construct_default_sevn_params() {
    SEVNpar sevn_par;
    for (const auto& v : sevn_par.get_num_map()) {
        default_params_sevn[v.first] = {
            std::to_string(v.second.first),
            v.second.second
        };
    }

    for (const auto& v : sevn_par.get_str_map()) {
        default_params_sevn[v.first] = v.second;
    }

    for (const auto& v : sevn_par.get_bool_map()) {
        default_params_sevn[v.first] = {
            v.second.first ? "true" : "false",
            v.second.second
        };
    }

    std::string unneeded_options[] = {
        "max_z", "max_z_he", "max_zams", "max_zams_he", "min_z", "min_z_he", "min_zams", "min_zams_he", "myself",
        "name_prefix", "o", "omode", "Z", "rseed", "list", "hard_kappa", "hard_mass_average", "hard_rhoc",
        "hard_sigma", "hard_xi", "hardmode", "ibmode", "scol", "bcol", "dtout", "check_stalling", "check_stalling_time",
        "initerror_stop", "io_literal_phases", "io_logfile", "io_logphase_bse", "nthreads", "spin", "tf", "tini",
        "eddington_factor", "rlo_max_nuclearmt", "sn_min_vkick", "star_lambda_pureHe"
    };
    for (auto& s : unneeded_options)
        default_params_sevn.erase(s);

    default_params_sevn.at("snmode").first = "rapid";
}

class PeTarSEVNIO : public IO {
public:
    PeTarSEVNIO(std::map<std::string, std::string>* params_sevn) {
        //IO class
        std::vector<std::string> args;
        std::vector<char*> argv;

        args.push_back(""); // Unused, refers to the paths
        for (auto& v : *params_sevn) {
            args.push_back("-" + v.first);
            args.push_back(v.second);
        }
        args.push_back("-io_logfile"); // Force no SEVN logging
        args.push_back("false");

        for (auto& s : args)
            argv.push_back(&s[0]);

        std::cout << "Loading SEVN tables" << std::endl;
        load(argv.size(), argv.data());
        std::cout << "SEVN tables loaded" << std::endl;
    }
};

PeTarSEVNIO& petarsevnio(std::map<std::string, std::string>* params_sevn) {
    static PeTarSEVNIO* instance = [params_sevn]() {
        return new PeTarSEVNIO(params_sevn); // Intentionally leaked
    }();
    return *instance;
}

/**
 * Transform a number to a string including the precision.
 * This is a replacement of the to_string method from the standard library
 * that has a constant precision of 6 (so for example it returns 0 for  numbers <1e-7).
 * @param T number
 * @param n precision
 */
template <typename T>
std::string number_to_string_with_precision(const T number, const int n = std::numeric_limits<T>::max_digits10) {
    std::ostringstream out;
    out.precision(n);
    out << number;
    return std::move(out).str();
}

std::string compose_mass_input(StarParameter* star) {
    assert(star->kw != 15);
    assert(star->RemnantType_SEVN != Lookup::Remnants::Empty);
    assert(star->Mzams_SEVN >= 0.0);

    std::string special_type;

    if (star->RemnantType_SEVN == Lookup::BH) special_type = "BH";
    else if (star->RemnantType_SEVN == Lookup::NS_CCSN) special_type = "NS";
    else if (star->RemnantType_SEVN == Lookup::NS_ECSN) special_type = "NSEC";
    else if (star->RemnantType_SEVN == Lookup::HeWD) special_type = "HEWD";
    else if (star->RemnantType_SEVN == Lookup::COWD) special_type = "COWD";
    else if (star->RemnantType_SEVN == Lookup::ONeWD) special_type = "ONEWD";


    std::string mass_input;
    //Mass
    if (!special_type.empty()) {
        mass_input = number_to_string_with_precision(star->mt) + special_type; // Initialise as remnants: BH, NS, WD
    } else if (star->tphys == 0.0) {
        mass_input = number_to_string_with_precision(star->Mzams_SEVN); // Initialise at ZAMS
    } else if (star->MHE_SEVN == star->mt) {
        double MCO = star->MCO_SEVN;
        if (star->MHE_SEVN == MCO) {
            MCO = std::nextafter(MCO, 0.0); // Kind of a hack, to initialise as pure CO stars
        }
        mass_input = "(" + number_to_string_with_precision(star->Mzams_SEVN) + "," +
            number_to_string_with_precision(star->mt) + "," +
            number_to_string_with_precision(MCO) + ")HE"; // Initialise as pure HE star
    } else {
        mass_input = "(" + number_to_string_with_precision(star->Mzams_SEVN) + "," +
            number_to_string_with_precision(star->mt) + ","
            + number_to_string_with_precision(star->MHE_SEVN) + "," + number_to_string_with_precision(
                star->MCO_SEVN) +
            ")"; // Default initialised
    }

    return mass_input;
}

std::string construct_plife_full(StarParameter* star) {
    if (star->Phase_SEVN == Lookup::Phases::Remnant) {
        return "%100:7";
    }
    return "%" + number_to_string_with_precision(star->Plife_SEVN * 100) + ":" + std::to_string(star->Phase_SEVN);
}

Star construct_sevnstar(StarParameter* star, const double* tphysf, const double* z,
                        std::map<std::string, std::string>* params_sevn,
                        size_t ID = 0) {

    if (!petarsevnio(params_sevn).tablesloaded) {
        throw std::runtime_error("Error in stellar evolution with SEVN: No tables loaded");
    }

    double dt = *tphysf - star->tphys;

    std::string _massinput, _z, _spin, _dt;

    // Compose the mass input (it can have suffixes) taking into account errors
    try {
        _massinput = compose_mass_input(star);
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';
        throw std::runtime_error("Error in stellar evolution with SEVN: Composing mass input");
    }

    _z = number_to_string_with_precision(*z);
    _spin = number_to_string_with_precision(star->ospin);
    _dt = number_to_string_with_precision(dt);

    std::string plife_full = construct_plife_full(star);

    // Set the input parameters
    std::vector<std::string> init_params{_massinput, _z, _spin, params_sevn->at("snmode"), plife_full, _dt, "events"};

    Star sevnstar = [&]() {
        try {
            return Star(&(petarsevnio(params_sevn)), init_params, ID, false);
        } catch (const std::exception& err) {
            //Catch all errors
            std::cerr << err.what() << '\n';
            std::cerr << "init_params = {";
            for (auto param : init_params) std::cerr << "\"" << param << "\", ";
            std::cerr << "}" << std::endl;
            throw std::runtime_error("Error in stellar evolution with SEVN: Error initialising star");
        }
    }();

    return sevnstar;
}

Binstar construct_sevnbinary(StarParameter* star1, StarParameter* star2, const double* semi_rsun, const double* ecc,
                             const double* tphysf, const double* z, std::map<std::string, std::string>* params_sevn) {

    if (!petarsevnio(params_sevn).tablesloaded) {
        throw std::runtime_error("Error in binary evolution with SEVN: No tables loaded");
    }

    double dt = *tphysf - star1->tphys;

    std::string _dt, _z, _semi_rsun, _ecc;

    _dt = number_to_string_with_precision(dt);
    _semi_rsun = number_to_string_with_precision(*semi_rsun);
    _ecc = number_to_string_with_precision(*ecc);

    std::string _massinput1 = compose_mass_input(star1);
    std::string _z1 = number_to_string_with_precision(*z);
    std::string _spin1 = number_to_string_with_precision(star1->ospin);
    std::string plife_full1 = construct_plife_full(star1);

    std::string _massinput2 = compose_mass_input(star2);
    std::string _z2 = number_to_string_with_precision(*z);
    std::string _spin2 = number_to_string_with_precision(star2->ospin);
    std::string plife_full2 = construct_plife_full(star2);

    std::vector<std::string> init_params{
        _massinput1, _z1, _spin1, params_sevn->at("snmode"), plife_full1,
        _massinput2, _z2, _spin2, params_sevn->at("snmode"), plife_full2,
        _semi_rsun, _ecc, _dt, "events"
    };

    //Set fake IDs
    size_t ID = 2;
    Binstar sevnbinary = [&]() {
        try {
            return Binstar(&(petarsevnio(params_sevn)), init_params, ID);
        } catch (const std::exception& err) {
            //Catch all errors
            std::cerr << err.what() << '\n';
            throw std::runtime_error("Error in binary evolution with SEVN: Error initialising binary");
        }
    }();

    return sevnbinary;
}

double get_timestep(StarParameter* star, double* z, std::map<std::string, std::string>* params_sevn) {
    try {
        SevnLogging svnlog = SevnLogging();
        svnlog.set_level(params_sevn->at("log_level"));
        Star sevnstar = construct_sevnstar(star, &VERY_LARGE_TIME, z, params_sevn);
        return sevnstar.getp(Timestep::ID);
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';
        throw std::runtime_error("Error in stellar evolution with SEVN: Error while computing stellar timestep");
    }
}

double get_timestep(StarParameter* star1, StarParameter* star2, const double* semi_rsun, double* ecc, double* z,
                    std::map<std::string, std::string>* params_sevn) {
    try {
        SevnLogging svnlog = SevnLogging();
        svnlog.set_level(params_sevn->at("log_level"));
        Binstar sevnbinary = construct_sevnbinary(star1, star2, semi_rsun, ecc, &VERY_LARGE_TIME, z, params_sevn);
        return sevnbinary.getp(BTimestep::ID);
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';
        throw std::runtime_error("Error in binary evolution with SEVN: Error while computing binary timestep");
    }
}

void record_star_properties(StarParameter* star, StarParameterOut* out, Star* sevnstar) {
    star->kw = sevnstar->getp(PhaseBSE::ID);
    star->mt = sevnstar->getp(Mass::ID);
    star->r = sevnstar->getp(Radius::ID);
    star->mc = sevnstar->getp(MCO::ID) > 0.0 ? sevnstar->getp(MCO::ID) : sevnstar->getp(MHE::ID);
    star->rc = sevnstar->getp(MCO::ID) > 0.0 ? sevnstar->getp(RCO::ID) : sevnstar->getp(RHE::ID);
    star->ospin = sevnstar->getp(Spin::ID);
    star->epoch = sevnstar->get_current_tphase();
    star->tphys += sevnstar->getp(Worldtime::ID);
    star->lum = sevnstar->getp(Luminosity::ID);
    star->Mzams_SEVN = sevnstar->get_zams();
    star->Plife_SEVN = (sevnstar->getp(Phase::ID) != Lookup::Phases::Remnant) ? sevnstar->plife() : 1.0;
    star->Phase_SEVN = sevnstar->getp(Phase::ID);
    star->MHE_SEVN = sevnstar->getp(MHE::ID);
    star->MCO_SEVN = sevnstar->getp(MCO::ID);
    star->RemnantType_SEVN = static_cast<int>(sevnstar->getp(RemnantType::ID));

    assert((star->RemnantType_SEVN == -1) == (star->kw == 15));

    out->menv = sevnstar->getp(Qconv::ID) * sevnstar->getp(Mass::ID);
    out->renv = sevnstar->getp(Depthconv::ID) * sevnstar->getp(Radius::ID);
    out->tm = sevnstar->get_next_tphase();

    if (sevnstar->vkick[3] > 0.0) {
        out->vkick[3] = sevnstar->vkick[3];
        out->vkick[0] = sevnstar->vkick[0];
        out->vkick[1] = sevnstar->vkick[1];
        out->vkick[2] = sevnstar->vkick[2];

        double vk2 = pow(sevnstar->vkick[0], 2) + pow(sevnstar->vkick[1], 2) + pow(sevnstar->vkick[2], 2);
        double vk2_alt = pow(sevnstar->vkick[3], 2);
        assert(abs(vk2 - vk2_alt)/vk2 < 1e-10);
    }

    if (star->kw == 15) {
        star->mt = 0.0;
        star->r = 0.0;
        star->mc = 0.0;
        star->rc = 0.0;
        star->ospin = 0.0;
        star->epoch = 0.0;
        star->lum = 0.0;
        star->MHE_SEVN = 0.0;
        star->MCO_SEVN = 0.0;
        star->RemnantType_SEVN = Lookup::Empty;
        for (double& vki : out->vkick) vki = 0.0;
    }

    assert(!std::isnan(out->vkick[0]));
    assert(!std::isnan(out->vkick[1]));
    assert(!std::isnan(out->vkick[2]));
    assert(!std::isnan(out->vkick[3]));
    assert(!std::isnan(star->mt));

    assert(!std::isinf(out->vkick[0]));
    assert(!std::isinf(out->vkick[1]));
    assert(!std::isinf(out->vkick[2]));
    assert(!std::isinf(out->vkick[3]));
    assert(!std::isinf(star->mt));
}

void evolv1_SEVN(StarParameter* star, StarParameterOut* out, double* tphysf, double* z,
                 std::map<std::string, std::string>* params_sevn) {
    SevnLogging svnlog = SevnLogging();
    svnlog.set_level(params_sevn->at("log_level"));
    if (star->Mzams_SEVN < 0.0) throw std::runtime_error("Error in stellar evolution with SEVN: Negative M_ZAMS");

    // In some rare instances, due to floating-point error, the timestep can be a (very small) negative number
    if ((*tphysf < 0.0) || (*tphysf - star->tphys < 0.0)) return;

    // If the timestep is 0, initialise and return
    if ((*tphysf == 0.0) || (*tphysf - star->tphys == 0.0)) {
        if (star->Plife_SEVN == 0.0) {
            double zero_timestep = 0.0;
            Star sevnstar = construct_sevnstar(star, &zero_timestep, z, params_sevn);
            record_star_properties(star, out, &sevnstar);
        }
        return;
    }

    if (star->kw == 15) {
        star->tphys = *tphysf;
        star->mt = 0.0;
        star->r = 0.0;
        star->mc = 0.0;
        star->rc = 0.0;
        star->ospin = 0.0;
        star->epoch = 0.0;
        star->lum = 0.0;
        star->MHE_SEVN = 0.0;
        star->MCO_SEVN = 0.0;
        star->RemnantType_SEVN = Lookup::Empty;
        for (double& vki : out->vkick) vki = 0.0;
        return;
    }

    Star sevnstar = construct_sevnstar(star, tphysf, z, params_sevn);

    try {
        for (;;) {
            sevnstar.evolve();
            if (sevnstar.breaktrigger()) {
                break;
            }
        }
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';

        double dt = *tphysf - star->tphys;

        std::string _massinput, _z, _spin, _dt;
        _massinput = compose_mass_input(star);
        _z = number_to_string_with_precision(*z);
        _spin = number_to_string_with_precision(star->ospin);
        _dt = number_to_string_with_precision(dt);
        std::string plife_full = construct_plife_full(star);

        // Set the input parameters
        std::cerr << "std::vector<std::string> init_params{\""
            << _massinput << "\", \""
            << _z << "\", \""
            << _spin << "\", \""
            << params_sevn->at("snmode") << "\", \""
            << plife_full << "\", \""
            << _dt << "\", \""
            << "events\"};" << std::endl;

        throw std::runtime_error("Error in stellar evolution with SEVN: Error while performing time evolution");
    }

    record_star_properties(star, out, &sevnstar);

}

std::vector<int> get_bse_event_id(Binstar* sevnbinary) {
    //  binary_type{
    // "Unset", //0
    // "Initial", //1
    // "Type_change", //2
    // "Start_Roche", //3
    // "End_Roche", //4
    // "Contact", //5
    // "Start_Symbiotic", //6
    // "End_Symbiotic", //7
    // "Common_envelope", //8
    // "Giant", //9
    // "Coalescence", //10
    // "Blue_straggler", //11
    // "No_remain", //12
    // "Disrupt" //13
    // }

    int eventID = round(sevnbinary->getp(BEvent::ID));
    std::vector<int> bse_event_id_list;

    if (sevnbinary->getstar(0)->getp(PhaseBSE::ID) == 15 and sevnbinary->getstar(1)->getp(PhaseBSE::ID) == 15) {
        return {12};
    }

    if (eventID == Lookup::EventsList::NoEvent) return {1};
    if (eventID == Lookup::EventsList::ChangePhase) return {2};
    if (eventID == Lookup::EventsList::ChangeRemnant) return {2};
    if (eventID == Lookup::EventsList::QHE) return {1}; // TODO: Not classified as an event in BSE, so not recorded
    if (eventID == Lookup::EventsList::GWBegin) return {1};
    if (eventID == Lookup::EventsList::RLOBegin) return {3};
    if (eventID == Lookup::EventsList::RLOEnd) return {4};
    if (eventID == Lookup::EventsList::Collision) return {5};
    if (eventID == Lookup::EventsList::CE) return {8};
    if (eventID == Lookup::EventsList::Merger) return {10};
    if (eventID == Lookup::EventsList::CE_Merger)return {8, 10};
    if (eventID == Lookup::EventsList::RLOB_Merger) return {3, 10};
    if (eventID == Lookup::EventsList::RLOB_CE) return {3, 8};
    if (eventID == Lookup::EventsList::RLOB_CE_Merger) return {3, 8, 10};
    if (eventID == Lookup::EventsList::Collision_Merger) return {5, 10};
    if (eventID == Lookup::EventsList::Collision_CE) return {5, 8};
    if (eventID == Lookup::EventsList::Collision_CE_Merger) return {5, 8, 10};
    if (eventID == Lookup::EventsList::Swallowed) return {10};
    if (eventID == Lookup::EventsList::RLOB_Swallowed) return {3, 10};
    if (eventID == Lookup::EventsList::GW_Merger) return {1};
    // GW_Merger are classified as NoEvent because they always trigger a Merger. Only the Merger is an event
    if (eventID == Lookup::EventsList::SNBroken) return {13};
    if (eventID == Lookup::EventsList::SNIa) return {13};
    if (eventID == Lookup::EventsList::RLOB_CE_Swallowed) return {3, 8, 10};
    if (eventID == Lookup::EventsList::CE_Swallowed) return {8, 10};

    throw std::runtime_error("Error in stellar evolution with SEVN: Unknown binary event");
}

void init_sevn_event(double (&bse_event)[33][9]) {
    for (int i = 0; i < BinaryEvent::getEventNMax(); i++) {
        for (int j = 0; j < 33; j++) {
            bse_event[j][i] = -1;
        }
    }
}

bool log_sevn_event(const double* tphys, Star* star0, Star* star1, const double radro0, double radro1,
                    const vector<int>& event_id_vec, int bevent,
                    double semi_rsun, double ecc, double (&bse_event)[33][9], int* current_event_idx) {
    bool retval = false;

    for (const int event_id : event_id_vec) {
        if (event_id <= 1 and *current_event_idx != BinaryEvent::getEventIndexInit()) {
            // No event and not logging initial status
            return false;
        }

        assert(BinaryEvent::getEventNMax()==BinaryEvent::getEventIndexInit());
        // The following code only works if the index of the initial event is the last one
        int nmax = BinaryEvent::getEventNMax() - 1;
        if (*current_event_idx == BinaryEvent::getEventIndexInit()) {
            nmax = BinaryEvent::getEventIndexInit();
        }
        for (int j = *current_event_idx; j <= nmax; j++) {
            bse_event[0][j] = *tphys + star0->getp(Worldtime::ID);
            bse_event[1][j] = star0->getp(Mass::ID);
            bse_event[2][j] = star1->getp(Mass::ID);
            bse_event[3][j] = star0->getp(PhaseBSE::ID);
            bse_event[4][j] = star1->getp(PhaseBSE::ID);
            bse_event[5][j] = semi_rsun;
            bse_event[6][j] = ecc;
            bse_event[7][j] = radro0;
            bse_event[8][j] = radro1;
            bse_event[9][j] = (j == *current_event_idx) ? event_id : -1;
            bse_event[10][j] = star0->getp(Luminosity::ID);
            bse_event[11][j] = star1->getp(Luminosity::ID);
            bse_event[12][j] = star0->getp(Radius::ID);
            bse_event[13][j] = star1->getp(Radius::ID);
            bse_event[14][j] = star0->getp(MCO::ID) > 0.0 ? star0->getp(MCO::ID) : star0->getp(MHE::ID);
            bse_event[15][j] = star1->getp(MCO::ID) > 0.0 ? star1->getp(MCO::ID) : star1->getp(MHE::ID);
            bse_event[16][j] = star0->getp(MCO::ID) > 0.0 ? star0->getp(RCO::ID) : star0->getp(RHE::ID);
            bse_event[17][j] = star1->getp(MCO::ID) > 0.0 ? star1->getp(RCO::ID) : star1->getp(RHE::ID);
            bse_event[18][j] = star0->getp(Spin::ID);
            bse_event[19][j] = star1->getp(Spin::ID);
            bse_event[20][j] = star0->get_zams();
            bse_event[21][j] = star1->get_zams();
            bse_event[22][j] = star0->getp(MHE::ID);
            bse_event[23][j] = star1->getp(MHE::ID);
            bse_event[24][j] = star0->getp(MCO::ID);
            bse_event[25][j] = star1->getp(MCO::ID);
            bse_event[26][j] = (star0->getp(Phase::ID) != Lookup::Phases::Remnant) ? star0->plife() : 1.0;
            bse_event[27][j] = (star1->getp(Phase::ID) != Lookup::Phases::Remnant) ? star0->plife() : 1.0;
            bse_event[28][j] = round(star0->getp(Phase::ID));
            bse_event[29][j] = round(star1->getp(Phase::ID));
            bse_event[30][j] = round(star0->getp(RemnantType::ID));
            bse_event[31][j] = round(star1->getp(RemnantType::ID));
            bse_event[32][j] = bevent;
        }

        if (*current_event_idx == BinaryEvent::getEventIndexInit()) {
            // Update the event index if it's not the initial one
            *current_event_idx = std::min(BinaryEvent::getEventNMax() - 1, *current_event_idx + 1);
        }

        retval |= BSEManager::isDisrupt(event_id);
        retval |= BSEManager::isMerger(event_id);
        retval |= BSEManager::isNoRemnant(event_id);
    }
    return retval;
}

void log_sevn_event(StarParameter* star1, StarParameter* star2, const double* tphys, const double* z, double* semi_rsun,
                    const double* ecc, double (&bse_event)[33][9], std::map<std::string, std::string>* params_sevn,
                    int* current_event_idx, int event_type) {
    vector<int> event_id_vec{event_type};

    try {
        if (*semi_rsun > 0.0 and *ecc >= 0.0 and *ecc < 1.0) {
            Binstar sevnbinary = construct_sevnbinary(star1, star2, semi_rsun, ecc, &VERY_LARGE_TIME, z, params_sevn);
            Star* sevnstar1 = sevnbinary.getstar(0);
            Star* sevnstar2 = sevnbinary.getstar(1);
            log_sevn_event(tphys, sevnstar1, sevnstar2,
                           sevnstar1->getp(Radius::ID) / sevnstar1->getp(RL0::ID),
                           sevnstar2->getp(Radius::ID) / sevnstar2->getp(RL0::ID),
                           event_id_vec, Lookup::EventsList::NoEvent, *semi_rsun,
                           *ecc, bse_event, current_event_idx);
            // We can not know what the previous bevent is. The temporary fix is to pass NoEvent, because
            // otherwise we would need to keep track of the bevent and significantly change the structure of the code
        } else {
            Star sevnstar1 = construct_sevnstar(star1, &VERY_LARGE_TIME, z, params_sevn, 0);
            Star sevnstar2 = construct_sevnstar(star2, &VERY_LARGE_TIME, z, params_sevn, 1);


            int logg_idx = BinaryEvent::getEventIndexInit();
            log_sevn_event(tphys, &sevnstar1, &sevnstar2, std::nan(""), std::nan(""),
                           event_id_vec, Lookup::EventsList::NoEvent, *semi_rsun, *ecc, bse_event, &logg_idx);
            // We can not know what the previous bevent is. The temporary fix is to pass NoEvent, because
            // otherwise we would need to keep track of the bevent and significantly change the structure of the code
        }
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';
        throw std::runtime_error("Error in binary evolution with SEVN: Error while logging SEVN event");
    }
}

void evolv2_SEVN(StarParameter* star1, StarParameter* star2, StarParameterOut* out1, StarParameterOut* out2,
                 double* tphysf, double* z, double* period_days, double* semi_rsun, double* ecc,
                 double (&bse_event)[33][9], std::map<std::string, std::string>* params_sevn) {
    SevnLogging svnlog = SevnLogging();
    svnlog.set_level(params_sevn->at("log_level"));

    assert(star1->Mzams_SEVN > 0.0);
    assert(star2->Mzams_SEVN > 0.0);

    Binstar sevnbinary = construct_sevnbinary(star1, star2, semi_rsun, ecc, tphysf, z, params_sevn);

    int current_event_idx = 0;

    // In some rare instances, due to floating-point error, the timestep can be a (very small) negative number
    if (*tphysf < 0.0) return;
    if (*tphysf - star1->tphys < 0.0) return;
    if (*tphysf - star2->tphys < 0.0) return;

    // If the timestep is 0, initialise and return
    if ((*tphysf == 0.0) || (*tphysf - star1->tphys == 0.0) || (*tphysf - star2->tphys == 0.0)) {
        if ((star1->Plife_SEVN == 0.0) || (star2->Plife_SEVN == 0.0)) {
            Star* sevnstar1 = sevnbinary.getstar(0);
            Star* sevnstar2 = sevnbinary.getstar(1);
            log_sevn_event(&star1->tphys,
                           sevnstar1, sevnstar2,
                           sevnstar1->getp(Radius::ID) / sevnstar1->getp(RL0::ID),
                           sevnstar2->getp(Radius::ID) / sevnstar2->getp(RL0::ID),
                           get_bse_event_id(&sevnbinary), round(sevnbinary.getp(BEvent::ID)),
                           sevnbinary.getp(Semimajor::ID), sevnbinary.getp(Eccentricity::ID),
                           bse_event, &current_event_idx);
            record_star_properties(star1, out1, sevnbinary.getstar(0));
            record_star_properties(star2, out2, sevnbinary.getstar(1));
        }
        return;
    }


    if (star1->RemnantType_SEVN == Lookup::Empty) {
        star1->tphys = *tphysf;
        throw std::runtime_error("Error in stellar evolution with SEVN: Massless remnant");
        return;
    }
    if (star2->RemnantType_SEVN == Lookup::Empty) {
        star2->tphys = *tphysf;
        throw std::runtime_error("Error in stellar evolution with SEVN: Massless remnant");
        return;
    }

    assert(abs(star1->tphys - star2->tphys) < 1e-10);

    try {
        for (;;) {
            sevnbinary.evolve();
            Star* sevnstar1 = sevnbinary.getstar(0);
            Star* sevnstar2 = sevnbinary.getstar(1);
            bool is_breaking_event = log_sevn_event(&star1->tphys,
                                                    sevnstar1, sevnstar2,
                                                    sevnstar1->getp(Radius::ID) / sevnstar1->getp(RL0::ID),
                                                    sevnstar2->getp(Radius::ID) / sevnstar2->getp(RL0::ID),
                                                    get_bse_event_id(&sevnbinary), round(sevnbinary.getp(BEvent::ID)),
                                                    sevnbinary.getp(Semimajor::ID), sevnbinary.getp(Eccentricity::ID),
                                                    bse_event, &current_event_idx);

            if (is_breaking_event or sevnbinary.breaktrigger()) {
                break;
            }
        }
    } catch (const std::exception& err) {
        // Catch all errors
        std::cerr << err.what() << '\n';

        double dt = *tphysf - star1->tphys;
        std::string _dt, _z, _semi_rsun, _ecc;

        _dt = number_to_string_with_precision(dt);
        _semi_rsun = number_to_string_with_precision(*semi_rsun);
        _ecc = number_to_string_with_precision(*ecc);

        std::string _massinput1 = compose_mass_input(star1);
        std::string _z1 = number_to_string_with_precision(*z);
        std::string _spin1 = number_to_string_with_precision(star1->ospin);
        std::string plife_full1 = construct_plife_full(star1);

        std::string _massinput2 = compose_mass_input(star2);
        std::string _z2 = number_to_string_with_precision(*z);
        std::string _spin2 = number_to_string_with_precision(star2->ospin);
        std::string plife_full2 = construct_plife_full(star2);

        std::vector<std::string> init_params{
            _massinput1, _z1, _spin1, params_sevn->at("snmode"), plife_full1,
            _massinput2, _z2, _spin2, params_sevn->at("snmode"), plife_full2,
            _semi_rsun, _ecc, _dt, "events"
        };

        std::cout << "std::vector<std::string> init_params{\""
            << _massinput1 << "\", \"" << _z1 << "\", \"" << _spin1 << "\", \"" << params_sevn->at("snmode")
            << "\", \"" << plife_full1 << "\", \"" << _massinput2 << "\", \"" << _z2 << "\", \"" << _spin2 << "\", \""
            << params_sevn->at("snmode") << "\", \"" << plife_full2 << "\", \"" << _semi_rsun << "\", \"" <<
            _ecc << "\", \"" << _dt << "\", \"" << "events\"};" << std::endl;

        throw std::runtime_error("Error in binary evolution with SEVN: Error while performing time evolution");
    }

    record_star_properties(star1, out1, sevnbinary.getstar(0));
    record_star_properties(star2, out2, sevnbinary.getstar(1));

    *period_days = sevnbinary.getp(Period::ID) * year_to_day;
    *ecc = sevnbinary.getp(Eccentricity::ID);
    *semi_rsun = sevnbinary.getp(Semimajor::ID);

    if (star1->kw == 15 || star2->kw == 15 || sevnbinary.broken) {
        *period_days = 0.0;
        *ecc = -1.0;
        *semi_rsun = 0.0;
    }
}

void merge_SEVN(StarParameter* star1, StarParameter* star2, StarParameterOut* out1, StarParameterOut* out2,
                double* z, double* semi_rsun, double* ecc, bool log_bse_event, double (&bse_event)[33][9],
                std::map<std::string, std::string>* params_sevn) {
    SevnLogging svnlog = SevnLogging();
    svnlog.set_level(params_sevn->at("log_level"));

    constexpr double SMALL_VALUE = 1e-10;

    double input_sma = *semi_rsun;
    double input_ecc = *ecc;

    if ((*semi_rsun <= 0) && (*ecc >= 1)) {
        // This approach just converts every merge into an equivalent parabolic merge (same periapsis), perhaps to
        // be improved in a future version of PeTar-SEVN
        input_ecc = 1 - SMALL_VALUE;
        input_sma = std::max(SMALL_VALUE, *semi_rsun * (1.0 - *ecc) / (1.0 - input_ecc));
    } else if (!((*semi_rsun > 0) && (0 <= *ecc) && (*ecc < 1))) {
        // Parabolic merge (with minor numerical issues)
        input_ecc = 1 - SMALL_VALUE;
        input_sma = SMALL_VALUE;
    }

    Binstar sevnbinary = construct_sevnbinary(star1, star2, &input_sma, &input_ecc, &VERY_LARGE_TIME, z,
                                              params_sevn);
    KollisionHurley coll;

    if (coll.check_collision(&sevnbinary)) {
        coll.evolve(&sevnbinary);
        sevnbinary.evolve();

        bool is_merger_successful = false;
        for (int i = 0; i < 3; i++) {
            is_merger_successful |= (sevnbinary.getstar(0)->getp(PhaseBSE::ID) == 15);
            is_merger_successful |= (sevnbinary.getstar(1)->getp(PhaseBSE::ID) == 15);

            if (is_merger_successful) break;
            sevnbinary.evolve();
        }

        if (!is_merger_successful) {
            throw std::runtime_error("Error in collision with SEVN: No merger after collision");
        }

        record_star_properties(star1, out1, sevnbinary.getstar(0));
        record_star_properties(star2, out2, sevnbinary.getstar(1));

        if (std::isnan(*ecc))
            assert(sevnbinary.broken);
        if (std::isnan(*semi_rsun))
            assert(sevnbinary.broken);

        *ecc = -1.0;
        *semi_rsun = 0.0;

        Star* sevnstar1 = sevnbinary.getstar(0);
        Star* sevnstar2 = sevnbinary.getstar(1);

        if (log_bse_event) {
            int current_event_idx = 0;
            log_sevn_event(&star1->tphys,
                           sevnstar1, sevnstar2,
                           sevnstar1->getp(Radius::ID) / sevnstar1->getp(RL0::ID),
                           sevnstar2->getp(Radius::ID) / sevnstar2->getp(RL0::ID),
                           get_bse_event_id(&sevnbinary), round(sevnbinary.getp(BEvent::ID)),
                           sevnbinary.getp(Semimajor::ID), sevnbinary.getp(Eccentricity::ID),
                           bse_event, &current_event_idx);
        }
    }
}


void printconst(std::map<std::string, std::string>* params_sevn) {
    std::cout << "----- SEVN parameter list: -----" << std::endl;
    std::cout << "Some of these options (used for SEVN I/O) might be unused in PeTar" << std::endl;
    for (auto& v : *params_sevn) {
        std::cout << v.first << ": " << v.second << std::endl;
    }
}