#pragma once

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include "Common/Float.h"
#include "static_variables.hpp"
#ifdef GALPY
#include "galpy_interface.h"
#include "status.hpp"
#include "soft_ptcl.hpp"
#endif

//! IO parameters manager for gas drag in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsGasDrag{
public:
    IOParamsContainer input_par_store;
    IOParams<long long int> mode; // option to switch perturbation
    IOParams<double> sound_speed;
    IOParams<double> coulomb_log;
    IOParams<double> polytropic_constant; // K
    IOParams<double> polytropic_exponent; // gamma
    IOParams<double> ifunc_mach_lower;
    IOParams<double> ifunc_mach_upper;
    IOParams<long long int> ifunc_smooth_order;
#ifdef GALPY
    IOParams<long long int> galpy_gaspot_index;
    IOParams<double> scale_density;
#else
    IOParams<double> gas_density;
    IOParams<double> decay_time;
#endif
    IOParams<double> gravitational_constant;
    IOParams<std::string> fname_par;

    bool print_flag;
    IOParamsGasDrag(): input_par_store(),
                       mode  (input_par_store, 0,          "gdf-hard-mode", "GDF external force for hard integration; 0, not used; 1, gas dynamical friction (Ostriker 1999, Rozner et al. 2022); 2, gas dynamical friction only in radial direction"),
                       sound_speed  (input_par_store, 0.0, "gdf-sound-speed",  "sound speed in units of PeTar input, if given 0, calculate by sqrt(P/rho), based on hydrostatic equilibrium"),
                       coulomb_log  (input_par_store, 3.1, "gdf-coulomb-log",  "coulomb logarithm"),
                       polytropic_constant (input_par_store, 1.0, "gdf-K", "Polytropic constant in units of PeTar input, used to evaluate Pressure P = K rho^gamma"),
                       polytropic_exponent (input_par_store, 4.0/3.0, "gdf-gamma", "Polytropic exponent, used to evaluate Pressure P = K rho^gamma"),
                       ifunc_mach_lower (input_par_store, 0.9, "gdf-ifunc-mach-lower", "lower Mach boundary for Ifunc Hermite interpolation; recommended range: about 0.85 (or smaller) for derivative order = 3, typical 0.9 for order = 2"),
                       ifunc_mach_upper (input_par_store, 1.1, "gdf-ifunc-mach-upper", "upper Mach boundary for Ifunc Hermite interpolation; recommended range: about 1.15 (or larger) for derivative order = 3, typical 1.1 for order = 2"),
                       ifunc_smooth_order (input_par_store, 2, "gdf-ifunc-smooth-order", "Ifunc Hermite smooth derivative order at boundaries (2 or 3); recommended: 2 for robustness, 3 with a wider Mach range (e.g. 0.85-1.15 or wider)"),
#ifdef GALPY
                       galpy_gaspot_index(input_par_store, -1, "gdf-gaspot-index",  "galpy potential set index for gas component, used for obtaining gas density", "None"),
                       scale_density(input_par_store, 1/G_ASTRO, "gdf-scale-density", "scale factor for galpy potential density","1/G"),
#else
                       gas_density  (input_par_store, 1.0, "gdf-gas-density",  "gas density in units of PeTar input"),
                       decay_time   (input_par_store, 0.0, "gdf-decay-time",  "gas density decay time scale in units of PeTar input, if 0, no decay"),
#endif
                       gravitational_constant (input_par_store, 1.0, "G", "Gravitational constant", NULL, false),
                       fname_par    (input_par_store, "input.par", "p", "Input parameter file for external force (this option should be used first before any other options)",NULL,false),
                       print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] print_format_info: if true, print the format information
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int ext_flag=-1;
        const struct option long_options[] = {
            {mode.key,    required_argument, &ext_flag, 0},
#ifdef GALPY
            {galpy_gaspot_index.key, required_argument, &ext_flag, 1},
            {scale_density.key, required_argument, &ext_flag, 2},
#else
            {gas_density.key, required_argument, &ext_flag, 1},
            {decay_time.key,  required_argument, &ext_flag, 2},
#endif
            {sound_speed.key, required_argument, &ext_flag, 3},
            {coulomb_log.key, required_argument, &ext_flag, 4},
            {polytropic_constant.key, required_argument, &ext_flag, 5},
            {polytropic_exponent.key, required_argument, &ext_flag, 6},
            {ifunc_mach_lower.key, required_argument, &ext_flag, 7},
            {ifunc_mach_upper.key, required_argument, &ext_flag, 8},
            {ifunc_smooth_order.key, required_argument, &ext_flag, 9},
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-G:p:h", long_options, &option_index)) != -1)
            switch (copt) {
            case 0:
                switch (ext_flag) {
                case 0:
                    mode.value = atoi(optarg);
                    if(print_flag) mode.print(std::cout);
                    opt_used+=2;
                    break;
#ifdef GALPY
                case 1:
                    galpy_gaspot_index.value = atoi(optarg);
                    if(print_flag) galpy_gaspot_index.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    scale_density.value = atof(optarg);
                    if(print_flag) scale_density.print(std::cout);
                    opt_used+=2;
                    break;
#else
                case 1:
                    gas_density.value = atof(optarg);
                    if(print_flag) gas_density.print(std::cout);
                    opt_used+=2;
                    break;
                case 2:
                    decay_time.value = atof(optarg);
                    if(print_flag) decay_time.print(std::cout);
                    opt_used+=2;
                    break;
#endif
                case 3:
                    sound_speed.value = atof(optarg);
                    if(print_flag) sound_speed.print(std::cout);
                    opt_used+=2;
                    break;
                case 4:
                    coulomb_log.value = atof(optarg);
                    if(print_flag) coulomb_log.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    polytropic_constant.value = atof(optarg);
                    if(print_flag) polytropic_constant.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    polytropic_exponent.value = atof(optarg);
                    if(print_flag) polytropic_exponent.print(std::cout);
                    opt_used+=2;
                    break;
                case 7:
                    ifunc_mach_lower.value = atof(optarg);
                    if(print_flag) ifunc_mach_lower.print(std::cout);
                    opt_used+=2;
                    break;
                case 8:
                    ifunc_mach_upper.value = atof(optarg);
                    if(print_flag) ifunc_mach_upper.print(std::cout);
                    opt_used+=2;
                    break;
                case 9:
                    ifunc_smooth_order.value = atoi(optarg);
                    if(print_flag) ifunc_smooth_order.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'G':
                gravitational_constant.value = atof(optarg);
                if(print_flag) gravitational_constant.print(std::cout);
                opt_used += 2;
                assert(gravitational_constant.value>0.0);
                break;
            case 'p':
                fname_par.value = optarg;
                if(print_flag) {
                    std::string fgalpy_par = fname_par.value+".exthard";
                    FILE* fpar_in;
                    if( (fpar_in = fopen(fgalpy_par.c_str(),"r")) == NULL) {
                        fprintf(stderr,"Error: Cannot open file %s.\n", fgalpy_par.c_str());
                        abort();
                    }
                    input_par_store.readAscii(fpar_in);
                    fclose(fpar_in);
                }
                opt_used+=2;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                input_par_store.mpi_broadcast();
                PS::Comm::barrier();
#endif
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"----- GDF external perturbation for hard integration options: -----"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }

        if(print_flag) std::cout<<"----- Finish reading input options of GDF external perturbation for hard integration -----\n";

        return opt_used;
    }
};

class GasDragForce{
public:
    int mode;
#ifdef GALPY
    int galpy_gaspot_index;
    GalpyManager* galpy_manager;
    Status* status;
    Float scale_density;
#else
    Float gas_density;
    Float gas_density_init;
    Float decay_time;
    Float time;
#endif
    Float sound_speed;
    Float coulomb_log;
    Float polytropic_constant;
    Float polytropic_exponent;
    Float gravitational_constant;
    bool calc_sound_speed;
    int ifunc_smooth_order;
    int ifunc_poly_ncoef;
    Float mach_inter_lower;
    Float mach_inter_upper;
    Float ifunc_coef[8];
    Float difunc_coef[7];

    GasDragForce(): mode(0),
#ifdef GALPY
                    galpy_gaspot_index(-1), galpy_manager(NULL), status(NULL), scale_density(1.0),
#else
                    gas_density(1.0), gas_density_init(1.0), decay_time(0.0), time(0.0),
#endif
                    sound_speed(0.0), coulomb_log(3.1), polytropic_constant(1.0), polytropic_exponent(4.0/3.0), gravitational_constant(1.0), calc_sound_speed(true),
                    ifunc_smooth_order(2), ifunc_poly_ncoef(6), mach_inter_lower(0.9), mach_inter_upper(1.1), ifunc_coef{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, difunc_coef{0.0,0.0,0.0,0.0,0.0,0.0,0.0} {}

    Float calcIfuncSubsonic(const Float _mach) const {
        return 0.5*std::log((1.0+_mach)/(1.0-_mach)) - _mach;
    }

    Float calcDIfuncSubsonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        return mach2/(1.0-mach2);
    }

    Float calcDDIfuncSubsonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        const Float den = 1.0 - mach2;
        return 2.0*_mach/(den*den);
    }

    Float calcDDDIfuncSubsonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        const Float den = 1.0 - mach2;
        return 2.0*(1.0+3.0*mach2)/(den*den*den);
    }

    Float calcIfuncSupersonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        return 0.5*std::log(1.0-1.0/mach2) + coulomb_log;
    }

    Float calcDIfuncSupersonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        return 1.0/(mach2*_mach - _mach);
    }

    Float calcDDIfuncSupersonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        const Float den = mach2*_mach - _mach;
        return -(3.0*mach2 - 1.0)/(den*den);
    }

    Float calcDDDIfuncSupersonic(const Float _mach) const {
        const Float mach2 = _mach*_mach;
        const Float den = mach2*_mach - _mach;
        const Float den2 = den*den;
        const Float den3 = den2*den;
        const Float dden = 3.0*mach2 - 1.0;
        return -6.0*_mach/den2 + 2.0*dden*dden/den3;
    }

    Float polynomialDerivativeBasis(const int _power, const int _derivative_order, const Float _x) const {
        if (_power<_derivative_order) return 0.0;
        Float coeff = 1.0;
        for (int i=0; i<_derivative_order; i++) coeff *= (_power-i);
        Float xpow = 1.0;
        for (int i=0; i<_power-_derivative_order; i++) xpow *= _x;
        return coeff*xpow;
    }

    Float evaluatePolynomial(const Float _x, const Float* _coef, const int _n_coef) const {
        ASSERT(_n_coef>0);
        Float value = _coef[_n_coef-1];
        for (int i=_n_coef-2; i>=0; i--) value = value*_x + _coef[i];
        return value;
    }

    void hermite_interpolation(const Float _x0, const Float _x1, const Float* _f0, const Float* _f1, const int _smooth_order, Float* _coef, int& _n_coef) {
        ASSERT(_smooth_order>=2 && _smooth_order<=3);
        const int n_eq = 2*(_smooth_order+1);
        const int n_col = n_eq + 1;
        Float mat[8][9] = {{0.0}};

        for (int d=0; d<=_smooth_order; d++) {
            const int row_l = d;
            const int row_r = d + _smooth_order + 1;
            for (int p=0; p<n_eq; p++) {
                mat[row_l][p] = polynomialDerivativeBasis(p, d, _x0);
                mat[row_r][p] = polynomialDerivativeBasis(p, d, _x1);
            }
            mat[row_l][n_eq] = _f0[d];
            mat[row_r][n_eq] = _f1[d];
        }

        for (int i=0; i<n_eq; i++) {
            int pivot = i;
            Float pivot_abs = std::abs(mat[i][i]);
            for (int j=i+1; j<n_eq; j++) {
                const Float cand_abs = std::abs(mat[j][i]);
                if (cand_abs > pivot_abs) {
                    pivot = j;
                    pivot_abs = cand_abs;
                }
            }
            ASSERT(pivot_abs>0.0);
            if (pivot!=i) {
                for (int k=i; k<n_col; k++) {
                    const Float tmp = mat[i][k];
                    mat[i][k] = mat[pivot][k];
                    mat[pivot][k] = tmp;
                }
            }

            const Float diag = mat[i][i];
            for (int k=i; k<n_col; k++) mat[i][k] /= diag;

            for (int j=0; j<n_eq; j++) {
                if (j==i) continue;
                const Float fac = mat[j][i];
                if (fac==0.0) continue;
                for (int k=i; k<n_col; k++) mat[j][k] -= fac*mat[i][k];
            }
        }

        _n_coef = n_eq;
        for (int i=0; i<n_eq; i++) _coef[i] = mat[i][n_eq];
    }

    void updateIfuncCoefficients(const bool _print_flag=false, const int _smooth_order=2) {
        const int smooth_order = (_smooth_order==3) ? 3 : 2;
        ifunc_smooth_order = smooth_order;

        const Float ifunc_l = calcIfuncSubsonic(mach_inter_lower);
        const Float difunc_l = calcDIfuncSubsonic(mach_inter_lower);
        const Float ddifunc_l = calcDDIfuncSubsonic(mach_inter_lower);
        const Float dddifunc_l = calcDDDIfuncSubsonic(mach_inter_lower);
        const Float ifunc_r = calcIfuncSupersonic(mach_inter_upper);
        const Float difunc_r = calcDIfuncSupersonic(mach_inter_upper);
        const Float ddifunc_r = calcDDIfuncSupersonic(mach_inter_upper);
        const Float dddifunc_r = calcDDDIfuncSupersonic(mach_inter_upper);

        const Float left[4] = {ifunc_l, difunc_l, ddifunc_l, dddifunc_l};
        const Float right[4] = {ifunc_r, difunc_r, ddifunc_r, dddifunc_r};
        for (int i=0; i<8; i++) ifunc_coef[i] = 0.0;
        for (int i=0; i<7; i++) difunc_coef[i] = 0.0;
        hermite_interpolation(mach_inter_lower, mach_inter_upper, left, right, smooth_order, ifunc_coef, ifunc_poly_ncoef);

        for (int i=0; i<ifunc_poly_ncoef-1; i++) {
            difunc_coef[i] = (i+1)*ifunc_coef[i+1];
        }

        if (_print_flag) {
            std::cout<<std::setprecision(18);
            std::cout<<"GasDrag Ifunc Hermite coefficients (smooth order="<<ifunc_smooth_order<<", coulomb_log="<<coulomb_log<<")"<<std::endl;
            std::cout<<"  Ifunc:";
            for (int i=0; i<ifunc_poly_ncoef; i++) std::cout<<" c"<<i<<"="<<ifunc_coef[i];
            std::cout<<std::endl;
            std::cout<<"  dIfunc:";
            for (int i=0; i<ifunc_poly_ncoef-1; i++) std::cout<<" c"<<i<<"="<<difunc_coef[i];
            std::cout<<std::endl;
        }
    }

#ifdef GALPY
    //! initial parameters for perturbation
    /*!
      @param[in] _input: input parameter
      @param[in] _galpy_manager: galpy manager pointer
      @param[in] _status: system status for information of time and pcm position and velocity offsets used for converting to galactic frame
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsGasDrag& _input, GalpyManager& _galpy_manager, Status& _status, bool _print_flag=false) {
        mode = _input.mode.value;
        galpy_gaspot_index = _input.galpy_gaspot_index.value;
        galpy_manager = &_galpy_manager;
        status = &_status;
        scale_density = _input.scale_density.value;
        sound_speed = _input.sound_speed.value;
        coulomb_log = _input.coulomb_log.value;
        polytropic_constant = _input.polytropic_constant.value;
        polytropic_exponent = _input.polytropic_exponent.value;
        mach_inter_lower = _input.ifunc_mach_lower.value;
        mach_inter_upper = _input.ifunc_mach_upper.value;
        ifunc_smooth_order = _input.ifunc_smooth_order.value;
        gravitational_constant = _input.gravitational_constant.value;
        if (sound_speed>0.0) calc_sound_speed = false;
        else calc_sound_speed = true;
        updateIfuncCoefficients(_print_flag, ifunc_smooth_order);
    }

#else
    //! initial parameters for perturbation
    void initial(const IOParamsGasDrag& _input, const Float _time, const bool _print_flag=false) {
        mode = _input.mode.value;
        gas_density_init = _input.gas_density.value;
        decay_time  = _input.decay_time.value;
        sound_speed = _input.sound_speed.value;
        coulomb_log = _input.coulomb_log.value;
        polytropic_constant = _input.polytropic_constant.value;
        polytropic_exponent = _input.polytropic_exponent.value;
        mach_inter_lower = _input.ifunc_mach_lower.value;
        mach_inter_upper = _input.ifunc_mach_upper.value;
        ifunc_smooth_order = _input.ifunc_smooth_order.value;
        gravitational_constant = _input.gravitational_constant.value;
        if (sound_speed>0.0) calc_sound_speed = false;
        else calc_sound_speed = true;
        updateIfuncCoefficients(_print_flag, ifunc_smooth_order);
        updateTime(_time);
    }

    //! update time and gas density
    void updateTime(const Float _time) {
        time = _time;
        if (decay_time>0.0) gas_density = gas_density_init * exp(-time/decay_time);
        else gas_density = gas_density_init;
    }
#endif

    //! External force for one particle in hard part
    /*!
      Gas dynamical friction
      Due to the acceleration dependence, this function must be used at the end of acceleration calculation
      (Ostriker 1999, https://ui.adsabs.harvard.edu/abs/1999ApJ...513..252O,
      Rozner 2022, https://arxiv.org/abs/2212.00807)

      @param[out] _acc0: acceleration
      @param[out] _acc1: jerk
      @param[in] _particle: particle data
      @param[in] _calc_acc1: if true, calculate jerk

      Return: the next integration time step (default: maximum floating point number)
    */
    template<class Tp>
    Float calcAccJerkExternal(Float* _acc0, Float* _acc1, const Tp& _particle, const bool _calc_acc1, const FPSoft& _center, const PS::S64 _center_id) {
        if (mode==0)
            return NUMERIC_FLOAT_MAX;

        // ignore center particle
        if (_particle.id == _center_id)
            return NUMERIC_FLOAT_MAX;

        auto& mass = _particle.mass;
        auto& pos = _particle.pos;
        auto& vel = _particle.vel;

        const Float PI = 4.0*atan(1.0);
        Float G = gravitational_constant;
        Float G2 = G*G;

        Float pos_rel[3] = {pos[0], pos[1], pos[2]};
        Float vel_rel[3] = {vel[0], vel[1], vel[2]};
#ifdef GALPY
        // in galactic frame, required by galpy and used to calculate radial direction
        if (_center_id>0) {
            // refer to center position and velocity
            pos_rel[0] -= _center.pos[0];
            pos_rel[1] -= _center.pos[1];
            pos_rel[2] -= _center.pos[2];

            vel_rel[0] -= _center.vel[0];
            vel_rel[1] -= _center.vel[1];
            vel_rel[2] -= _center.vel[2];
        }
        else {
            // refer to gas potential center position and velocity
            Float pot_pos[3];
            galpy_manager->getSetPos(galpy_gaspot_index, pot_pos);
            pos_rel[0] += status->pcm.pos[0] - pot_pos[0];
            pos_rel[1] += status->pcm.pos[1] - pot_pos[1];
            pos_rel[2] += status->pcm.pos[2] - pot_pos[2];

            // in gas potential center reference
            Float pot_vel[3];
            galpy_manager->getSetVel(galpy_gaspot_index, pot_vel);
            vel_rel[0] += status->pcm.vel[0] - pot_vel[0];
            vel_rel[2] += status->pcm.vel[1] - pot_vel[1];
            vel_rel[3] += status->pcm.vel[2] - pot_vel[2];
        }

        // in galactic frame, required by galpy and used to calculate radial direction
        Float pos_g[3] = {pos[0] + status->pcm.pos[0],
                          pos[1] + status->pcm.pos[1],
                          pos[2] + status->pcm.pos[2]};
        Float gas_density = scale_density*galpy_manager->calcSetDensity(galpy_gaspot_index, status->time, pos_g, &pos[0]);
#else
        if (_center_id>0) {
            // refer to center position and velocity
            pos_rel[0] -= _center.pos[0];
            pos_rel[1] -= _center.pos[1];
            pos_rel[2] -= _center.pos[2];

            vel_rel[0] -= _center.vel[0];
            vel_rel[1] -= _center.vel[1];
            vel_rel[2] -= _center.vel[2];
        }
#endif

        Float r_rel = std::sqrt(pos_rel[0]*pos_rel[0] + pos_rel[1]*pos_rel[1] + pos_rel[2]*pos_rel[2]);

        // subtract keplerian velocity orbiting around the center
        if (_center_id>0) {
            // circular velocity
            Float v_circle = std::sqrt(G*_center.mass/r_rel);
            // get angular momentum direction
            Float r_cross_v[3] = {pos_rel[1]*vel_rel[2] - pos_rel[2]*vel_rel[1],
                                  pos_rel[2]*vel_rel[0] - pos_rel[0]*vel_rel[2],
                                  pos_rel[0]*vel_rel[1] - pos_rel[1]*vel_rel[0]};
            // get tangent velocity direction
            Float r_cross_v_cross_r[3] = {r_cross_v[1]*pos_rel[2] - r_cross_v[2]*pos_rel[1],
                                          r_cross_v[2]*pos_rel[0] - r_cross_v[0]*pos_rel[2],
                                          r_cross_v[0]*pos_rel[1] - r_cross_v[1]*pos_rel[0]};
            // normalization factor
            Float r2vsq = std::sqrt(r_cross_v_cross_r[0]*r_cross_v_cross_r[0] + r_cross_v_cross_r[1]*r_cross_v_cross_r[1] + r_cross_v_cross_r[2]*r_cross_v_cross_r[2]);

            // remove circular velocity
            vel_rel[0] -= v_circle*r_cross_v_cross_r[0]/r2vsq;
            vel_rel[1] -= v_circle*r_cross_v_cross_r[1]/r2vsq;
            vel_rel[2] -= v_circle*r_cross_v_cross_r[2]/r2vsq;
        }

        if (calc_sound_speed)
            sound_speed = std::sqrt(polytropic_constant*std::pow(gas_density, polytropic_exponent-1));

        Float v2 = vel_rel[0]*vel_rel[0] + vel_rel[1]*vel_rel[1] + vel_rel[2]*vel_rel[2];
        Float v = std::sqrt(v2);
        Float cs2 = sound_speed*sound_speed;

        Float mach = v/sound_speed;
        Float Ifunc, dIfunc;
        if (mach<mach_inter_lower) {
            Ifunc = calcIfuncSubsonic(mach);
            dIfunc = calcDIfuncSubsonic(mach);
        }
        else if (mach>=mach_inter_lower && mach<mach_inter_upper) {
            // Hermite interpolation with configurable smooth order (C2 or C3).
            // Recommended setup: order=2 with 0.9-1.1 for robustness; order=3 with 0.85-1.15 for smoother high-order derivatives.
            Ifunc = evaluatePolynomial(mach, ifunc_coef, ifunc_poly_ncoef);
            dIfunc = evaluatePolynomial(mach, difunc_coef, ifunc_poly_ncoef-1);
        }
        else{
            Ifunc = calcIfuncSupersonic(mach);
            dIfunc = calcDIfuncSupersonic(mach);
        }


        // More reasonable formulation with v^2+cs^2 in the denominator,
        // which behaves better in the subsonic regime.
        // The original Ostriker 1999 formulation with v^3 in the denominator is recovered in the supersonic limit.
        Float v2_cs2 = v2 + cs2;
        Float c1 = -4*PI*G2*mass*gas_density*v/(v2_cs2*v2_cs2)*Ifunc;
//        Float v3 = v2*v;
//        Float c1 = -4*PI*G2*mass*gas_density/v3*Ifunc;

        if (mode==1) {
            // GDF force
            _acc0[0] += c1*vel_rel[0];
            _acc0[1] += c1*vel_rel[1];
            _acc0[2] += c1*vel_rel[2];
        }
        else{
            // GDF force only in radial direction
            if (r_rel>0) {
                Float vel_r = (vel_rel[0]*pos_rel[0] + vel_rel[1]*pos_rel[1] + vel_rel[2]*pos_rel[2])/r_rel;
                _acc0[0] += c1*vel_r*pos_rel[0]/r_rel;
                _acc0[1] += c1*vel_r*pos_rel[1]/r_rel;
                _acc0[2] += c1*vel_r*pos_rel[2]/r_rel;
            }
        }

        if (_calc_acc1) {
            ASSERT(_acc1!=NULL);
            Float vdota = vel_rel[0]*_acc0[0] + vel_rel[1]*_acc0[1] + vel_rel[2]*_acc0[2];

            // d(v/(v^2+cs^2)^2)/dt = (cs^2 - 3v^2)/(v^2+cs^2)^5(v^2) v dot a
            Float c2 = c1*(cs2 - 3*v2)/(v2_cs2*v2)*vdota;

            // d(1/v^3)/dt = -3/v^5 v dot a
            //Float c2 = -3*c1/v2*vdota;

            // d(v/ds)/dt  = v dot a / (v*ds)
            Float c3 = c1/Ifunc*dIfunc*vdota/(v*sound_speed);

            if (mode==1) {
                // GDF force derivative
                _acc1[0] += (c2+c3)*vel_rel[0] + c1*_acc0[0];
                _acc1[1] += (c2+c3)*vel_rel[1] + c1*_acc0[1];
                _acc1[2] += (c2+c3)*vel_rel[2] + c1*_acc0[2];
            }
            else if (mode==2) {
                if (r_rel>0) {
                    // GDF force derivative only in radial direction, assuming pos_rel is constant. For time-dependent pos_rel, additional term of time derivative of pos_rel should be added.
                    Float vel_r = (vel_rel[0]*pos_rel[0] + vel_rel[1]*pos_rel[1] + vel_rel[2]*pos_rel[2])/r_rel;
                    Float a_r = (_acc1[0]*pos_rel[0] + _acc1[1]*pos_rel[1] + _acc1[2]*pos_rel[2])/r_rel;
                    Float pos_fac = (a_r + (v2 - 2*vel_r*vel_r)/r_rel);
                    _acc1[0] += ((c2+c3)*vel_r*pos_rel[0] + c1*(pos_fac*pos_rel[0] + vel_r*vel_rel[0]))/r_rel;
                    _acc1[1] += ((c2+c3)*vel_r*pos_rel[1] + c1*(pos_fac*pos_rel[1] + vel_r*vel_rel[1]))/r_rel;
                    _acc1[2] += ((c2+c3)*vel_r*pos_rel[2] + c1*(pos_fac*pos_rel[2] + vel_r*vel_rel[2]))/r_rel;
                }
            }
        }

        return NUMERIC_FLOAT_MAX;
    }

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(mode>=0 && mode<=2);
        if (mode>0) {
            ASSERT(sound_speed>=0.0);
            ASSERT(mach_inter_lower>0.0);
            ASSERT(mach_inter_upper>mach_inter_lower);
            ASSERT(ifunc_smooth_order==2 || ifunc_smooth_order==3);
#ifdef GALPY
            ASSERT(galpy_gaspot_index>=0);
            ASSERT(galpy_manager!=NULL);
            ASSERT(status!=NULL);
#endif
        }
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"----- External perturbation for hard integration parameters -----\n"
             <<"external mode: "<<mode<<std::endl
             <<"sound speed: "<<sound_speed<<std::endl
             <<"coulomb log: "<<coulomb_log<<std::endl
             <<"polytropic constant: "<<polytropic_constant<<std::endl
             <<"polytropic exponent: "<<polytropic_exponent<<std::endl
             <<"ifunc mach lower: "<<mach_inter_lower<<std::endl
             <<"ifunc mach upper: "<<mach_inter_upper<<std::endl
             <<"ifunc smooth order: "<<ifunc_smooth_order<<std::endl;
#ifdef GALPY
        _fout<<"galpy gaspot index: "<<galpy_gaspot_index<<std::endl
             <<"scale density: "<<scale_density<<std::endl;
#else
        _fout<<"gas density init: "<<gas_density_init<<std::endl
             <<"decay time: "<<decay_time<<std::endl;
#endif
        _fout<<"gravitational constant: "<<gravitational_constant<<std::endl
             <<"----- Finish reading input options of external perturbation for hard integration -----\n";
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
    }

    void printColumnBinary(std::ostream& _fout) const {
        _fout.write(reinterpret_cast<const char*>(this), sizeof(*this));
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = fread(this, sizeof(*this), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void readBinary(std::istream& _fin) {
        _fin.read(reinterpret_cast<char*>(this), sizeof(*this));
        if (!_fin) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1.\n";
            abort();
        }
    }
};
