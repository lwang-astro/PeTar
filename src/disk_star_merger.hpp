#pragma once
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <getopt.h>
#include "io.hpp"
#include <cassert>
#include "Common/Float.h"

//! Class for managing mergers
enum class StarType:int {none = -1, smbh = 0, bh = 1, star = 2, seed = 3, star_remnant = 4, bh_remnant = 5};

//! IO parameters manager for external perturbation in hard integration
/*! For initializing the COMMON block variables from the commander option.
  The description of each parameter is also provided.
 */
class IOParamsDiskStarMerger{
public:
    IOParamsContainer input_par_store;
    IOParams<double> merger_mass_loss_rate; //!< mass loss rate after merger
    IOParams<double> stellar_radius_power_index; //!< stellar radius power index
    IOParams<double> stellar_radius_scale; //!< stellar radius scale
    IOParams<double> merger_time_delay; //!< time delay for merger to increase mass and change radius
    IOParams<double> initial_mass; //!< initial mass of star seed
    IOParams<double> target_mass; //!< equilibrium mass
    IOParams<double> mass_growth_factor; //!< mass growth factor (c) for increasing mass to equlibrium (dM/dt = c M^2)
    IOParams<double> mass_loss_rate; //!< mass decrease rate dM/dt for decreasing mass to equlibrium 
    IOParams<long long int> redistribute_star_mode; //!< redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity to opposite side of the center

    bool print_flag; //!< print flag
    //! Constructor
    IOParamsDiskStarMerger(): input_par_store(),
                              merger_mass_loss_rate(input_par_store, 0.0, "merger-mass-loss-rate", "mass loss rate for merger"),
                              stellar_radius_power_index(input_par_store, 0.6, "stellar-radius-power", "stellar radius power index 'n', rs = s M^n"),
                              stellar_radius_scale(input_par_store, 0.0046, "stellar-radius-scale", "stellar radius scale 's', rs = s M^n"),
                              merger_time_delay(input_par_store, 0.0, "merger-time-delay", "time delay for merger to increase mass"),
                              initial_mass(input_par_store, 10.0, "initial-mass", "initial mass for star approaching equilibrium"),
                              target_mass(input_par_store, 300.0, "target-mass", "mass for star approaching equilibrium"),
                              mass_growth_factor(input_par_store, 0.0, "mass-growth-factor", "mass growth factor (c) for increasing mass to equlibrium (dM/dt = c M^2)"),
                              mass_loss_rate(input_par_store, 0.0, "mass-loss-rate", "mass loss rate dM/dt for decreasing mass to equlibrium"), 
                              redistribute_star_mode(input_par_store, 1, "redistribute-star-mode", "redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity to opposite side of the center; 2: redistribute star by choosing next type 3 star"),
                              print_flag(false) {}

    //! reading parameters from GNU option API
    /*!
      @param[in] argc: number of options
      @param[in] argv: string of options
      @param[in] opt_used_pre: already used option number from previous reading, use to correctly count the remaining argument number
      \return -1 if help is used; else the used number of argv
     */
    int read(int argc, char *argv[], const bool print_format_info=true, const int opt_used_pre=0) {
        static int merger_flag=-1;
        const struct option long_options[] = {
            {merger_mass_loss_rate.key, required_argument, &merger_flag, 0},  
            {stellar_radius_power_index.key, required_argument, &merger_flag, 1},  
            {stellar_radius_scale.key, required_argument, &merger_flag, 2},
            {merger_time_delay.key, required_argument, &merger_flag, 3},  
            {initial_mass.key, required_argument, &merger_flag, 4},
            {target_mass.key, required_argument, &merger_flag, 5},
            {mass_growth_factor.key, required_argument, &merger_flag, 6},
            {mass_loss_rate.key, required_argument, &merger_flag, 7},
            {redistribute_star_mode.key, required_argument, &merger_flag, 8},
            {"help",      no_argument,       0, 'h'},
            {0,0,0,0}
        };

        int opt_used=opt_used_pre;
        int copt;
        int option_index;
        optind = 0;
        while ((copt = getopt_long(argc, argv, "-z:p:h", long_options, &option_index)) != -1) 
            switch (copt) {
            case 0:
                switch (merger_flag) {
                case 0:
                    merger_mass_loss_rate.value = atof(optarg);
                    if(print_flag) merger_mass_loss_rate.print(std::cout);
                    opt_used+=2;
                    break;            
                case 1:
                    stellar_radius_power_index.value = atof(optarg);
                    if(print_flag) stellar_radius_power_index.print(std::cout);
                    opt_used+=2;
                    break;            
                case 2:
                    stellar_radius_scale.value = atof(optarg);
                    if(print_flag) stellar_radius_scale.print(std::cout);
                    opt_used+=2;
                    break;
                case 3:
                    merger_time_delay.value = atof(optarg);
                    if(print_flag) merger_time_delay.print(std::cout);
                    opt_used+=2;
                    break;
                case 4:
                    initial_mass.value = atof(optarg);
                    if(print_flag) initial_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 5:
                    target_mass.value = atof(optarg);
                    if(print_flag) target_mass.print(std::cout);
                    opt_used+=2;
                    break;
                case 6:
                    mass_growth_factor.value = atof(optarg);
                    if(print_flag) mass_growth_factor.print(std::cout);
                    opt_used+=2;
                    break;
                case 7:
                    mass_loss_rate.value = atof(optarg);
                    if(print_flag) mass_loss_rate.print(std::cout);
                    opt_used+=2;
                    break;
                case 8:
                    redistribute_star_mode.value = atoi(optarg);
                    if(print_flag) redistribute_star_mode.print(std::cout);
                    opt_used+=2;
                    break;
                default:
                    break;
                }
                break;
            case 'h':
                if(print_flag){
                    std::cout<<"Disk star merger parameters, options:"<<std::endl;
                    input_par_store.printHelp(std::cout, print_format_info);
                }
                return -1;
            case '?':
                opt_used +=2;
                break;
            default:
                break;
            }
        
        if(print_flag) std::cout<<"----- Finish reading input options of disk star merger parameters -----\n";

        return opt_used;
    }
};

class DiskStarMergerManager {
public:
    Float merger_mass_loss_rate; //!< mass loss rate
    Float stellar_radius_power_index; //!< radius amplifier
    Float stellar_radius_scale; //!< radius scale
    Float merger_time_delay; //!< time delay for merger to increase mass
    Float initial_mass; //!< initial mass
    Float target_mass; //!< equilibrium mass
    Float mass_growth_factor; //!< mass growth factor (c) for increasing mass to equlibrium (dM/dt = c M^2)
    Float mass_loss_rate; //!< mass decrease rate dM/dt for decreasing mass to equlibrium
    int redistribute_star_mode; //!< redistribute star mode, 0: no redistribute; 1: redistribute star position and velocity to opposite side of the center
    
    //! Constructor
    DiskStarMergerManager(): merger_mass_loss_rate(0.0), stellar_radius_power_index(0.6), stellar_radius_scale(0.0046), merger_time_delay(0.0), target_mass(0.0), mass_growth_factor(0.0), mass_loss_rate(0.0), redistribute_star_mode(1) {}

    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        assert(merger_mass_loss_rate>=0.0);
        assert(stellar_radius_power_index>=0.0);
        assert(stellar_radius_scale>=0.0);
        assert(merger_time_delay>=0.0);
        assert(initial_mass>=0.0);
        assert(target_mass>=0.0);
        assert(mass_growth_factor>=0.0);
        assert(mass_loss_rate>=0.0);
        assert(redistribute_star_mode>=0 && redistribute_star_mode<=2);
        return true;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"merger_mass_loss_rate : "<<merger_mass_loss_rate<<std::endl
             <<"stellar_radius_power_index : "<<stellar_radius_power_index<<std::endl
             <<"stellar_radius_scale : "<<stellar_radius_scale<<std::endl
             <<"merger_time_delay : "<<merger_time_delay<<std::endl
             <<"initial_mass : "<<initial_mass<<std::endl
             <<"target_mass : "<<target_mass<<std::endl
             <<"mass_growth_factor : "<<mass_growth_factor<<std::endl
             <<"mass_loss_rate : "<<mass_loss_rate<<std::endl
             <<"redistribute_star_mode : "<<redistribute_star_mode<<std::endl;
    }

    //! initial parameters for disk star mergers
    /*!
      @param[in] _input: input parameter
      @param[in] _print_flag: printing flag
     */
    void initial(const IOParamsDiskStarMerger& _input, const bool _print_flag=false) {
        merger_mass_loss_rate = _input.merger_mass_loss_rate.value;
        stellar_radius_power_index = _input.stellar_radius_power_index.value;
        stellar_radius_scale = _input.stellar_radius_scale.value;
        merger_time_delay = _input.merger_time_delay.value;
        initial_mass = _input.initial_mass.value;
        target_mass = _input.target_mass.value;
        mass_growth_factor = _input.mass_growth_factor.value;
        mass_loss_rate = _input.mass_loss_rate.value;
        redistribute_star_mode = _input.redistribute_star_mode.value;
    }

    //! calcMergerProperties
    /*! Calculate the new properties of the merger
        @param[in,out] p1: particle 1, will be merger 
        @param[in] p2: particle 2, will be zero mass particle
        @param[in] time: current time
        */
    template <class TParticle>
    void calcMergerProperties(TParticle* p1, TParticle* p2, const Float& time) {

        TParticle *pm, *p0; // set final merger to star or first BH
        if (p1->star.getType()==StarType::bh && p2->star.getType()==StarType::star) {
            pm = p2;
            p0 = p1;
        }
        else if(p1->star.getType()==StarType::seed) {
            pm = p2;
            p0 = p1;
        }
        else if(p1->star.getType()!=StarType::smbh && p2->star.getType()==StarType::smbh) {
            pm = p2;
            p0 = p1;
        }
        else {
            pm = p1;
            p0 = p2;
        }

        Float mcm = p1->mass + p2->mass;
        for (int k=0; k<3; k++) {
            pm->pos[k] = (p1->mass*p1->pos[k] + p2->mass*p2->pos[k])/mcm;
            pm->vel[k] = (p1->mass*p1->vel[k] + p2->mass*p2->vel[k])/mcm;
            p0->pos[k] = pm->pos[k]*(1+1e-8)+1e-12;
            p0->vel[k] = pm->vel[k]*(1+1e-8)+1e-12;
        }
        
        // only increase mass and change radius after time delay
        if (time > merger_time_delay) {
            Float new_mass = mcm * (1 - merger_mass_loss_rate);    
            pm->dm += new_mass - pm->mass;
            pm->mass = new_mass;
            pm->radius = std::pow(new_mass, stellar_radius_power_index);
        }

        p0->dm -= p0->mass;
        p0->mass = 0.0;
        p0->radius = 0.0;

        if (p0->star.getType() == StarType::bh) {
            pm->star.n_merger_bh++;
            p0->star.setType(StarType::bh_remnant);
        }
        else if (p0->star.getType() == StarType::star) {
            pm->star.n_merger_star++;
            pm->star.n_merger_bh += p0->star.n_merger_bh;
            p0->star.setType(StarType::star_remnant);
        }

        pm->star.last_merger_time = time;
        pm->star.last_mass_change_time = time;
        p0->star.last_mass_change_time = time;
        p0->star.n_merger_star = 0;
        p0->star.n_merger_bh = 0;
    }

    //! calculate mass change
    /*! Calculate mass change, for mass < target mass, increase mass; for mass > equilibrium mass, decrease mass
        Increase mass formula:  
        dM/dt = c M^2
        c: mass growth factor
        @param[in] p: particle
        @param[in] time: current time

        \return 0: no change; 1: modified mass
        */
    template <class TParticle>
    int calcMassChange(TParticle* p, const Float& time) {
        // if type is star, evolve mass to equilibrium mass    
        if (p->star.getType()==StarType::star) {
            Float dt = time - p->star.last_mass_change_time;
            if (dt>0) {
                // increase mass
                if (p->mass < target_mass && mass_growth_factor > 0) {
                    Float new_mass = p->mass + mass_growth_factor*p->mass*p->mass*dt;
                    if (new_mass > target_mass) {
                        new_mass = target_mass;
                    }
                    p->dm += new_mass - p->mass;
                    p->mass = new_mass;
                    p->radius = stellar_radius_scale * std::pow(new_mass, stellar_radius_power_index);
                    p->star.last_mass_change_time = time;
                    return 1;
                }
                // decrease mass
                else if (p->mass > target_mass && mass_loss_rate > 0) {
                    Float new_mass = p->mass - mass_loss_rate * dt;
                    if (new_mass < target_mass) {
                        new_mass = target_mass;
                    }
                    p->dm += new_mass - p->mass;
                    p->mass = new_mass;
                    p->radius = stellar_radius_scale * std::pow(new_mass, stellar_radius_power_index);
                    p->star.last_mass_change_time = time;
                    return 1;
                }
            }         
        }

        return 0;
    }
    

    //! redistribute star position and velocity
    /*!
        redistribute star position and velocity to opposite side of the center     
        @param[in,out] p: particle to redistribute
        @param[in] center: center particle
        @param[in] system: particle system for search star seed
        @param[in] n_system: number of particles in system
        @return 0: no redistribute or need to delete star; 1: redistribute star
    */
    template <class TParticle>
    int redistributeStar(TParticle* p, TParticle* center, TParticle* system, const int n_system) {
        // if type is star, redistribute 
        if (p->star.getType() == StarType::star_remnant) {
            if (redistribute_star_mode == 1) {
                Float pos[3];
                Float vel[3];
                for (int k=0; k<3; k++) {
                    pos[k] = p->pos[k] - center->pos[k];
                    vel[k] = p->vel[k] - center->vel[k];
                }
                Float r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
                if (r>0) {
                    for (int k=0; k<3; k++) {
                        p->pos[k] = center->pos[k] - pos[k];
                        p->vel[k] = center->vel[k] - vel[k];
                    }
                    p->star.setType(StarType::star);
                    p->mass = initial_mass;
                    p->dm = 0.0;
                    p->radius = stellar_radius_scale * std::pow(initial_mass, stellar_radius_power_index);
                    return 1;
                }
            }
            else if (redistribute_star_mode == 2) {
                // redistribute star by choosing a star seed with a similar distance to center
                Float pos[3] = {p->pos[0] - center->pos[0], 
                                p->pos[1] - center->pos[1], 
                                p->pos[2] - center->pos[2]};
                Float r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
                // Search next 10 star seeds and pickup one with the closest r
                int index_seed = -1;
                Float dr = NUMERIC_FLOAT_MAX;
                int n_found = 0;
                for (int i=0; i<n_system; i++) {
                    auto pi = &system[i];
                    if (pi->star.getType() == StarType::seed) {
                        Float pos_seed[3] = {pi->pos[0] - center->pos[0], 
                                             pi->pos[1] - center->pos[1], 
                                             pi->pos[2] - center->pos[2]};
                        Float ri = std::sqrt(pos_seed[0]*pos_seed[0] + pos_seed[1]*pos_seed[1] + pos_seed[2]*pos_seed[2]);
                        Float dri = std::abs(r - ri);
                        if (dri<dr) {
                            dr = dri;
                            index_seed = i;
                        }
                        n_found++;
                        if (n_found>=10) break;
                    }
                }
                if (index_seed >= 0) {
                    auto pi = &system[index_seed];
                    pi->star.setType(StarType::star);
                    pi->mass = initial_mass;
                    pi->dm = 0.0;
                    pi->radius = stellar_radius_scale * std::pow(initial_mass, stellar_radius_power_index);
                    pi->star.last_mass_change_time = p->star.last_mass_change_time;
                        
                    pi->star.n_merger_star = 0;
                    pi->star.n_merger_bh = 0;
                    return 0;
                }
            }
        }
        return 0;
    } 

};


//! class for disk star merger parameters of individual stars
class StarParameter{
public:
    long long int type; //!< type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
    long long int n_merger_star; //!< times of merger with star
    long long int n_merger_bh; //!< times of merger with black hole
    Float last_mass_change_time; //!< time delay for mass approach target
    Float last_merger_time; //!< last merger time

    //! Constructor
    StarParameter():  type(-1),
                      n_merger_star(0),
                      n_merger_bh(0),
                      last_mass_change_time(0.0),
                      last_merger_time(0.0)                    
                      {}

    //! set type
    /*!
      @param[in] _type: type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
     */
    void setType(const StarType _type) {
        type = static_cast<long long int>(_type);
    }

    //! get type
    /*!
      @return type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
     */
    StarType getType() const {
        return static_cast<StarType>(type);
    }

    //! initial parameters for disk star merger
    /*!
      @param[in] _type: type of object; 0: supermassive black hole; 1: black hole; 2: star; 3: star seed; 4: star zero mass remnant; 5: black hole zero mass remnant
      @param[in] _last_mass_change_time: time delay for mass approach target
     */
    void initial(const StarType _type, const Float& _last_mass_change_time=0.0) {
        type = static_cast<long long int>(_type);
        n_merger_star = 0;
        n_merger_bh = 0;
        last_mass_change_time = _last_mass_change_time;
        last_merger_time = 0;
    }

    //! write class data with ASCII format
    void writeAscii(FILE* fp) const{
        fprintf(fp, "%lld %lld %lld %26.17e %26.17e\n", type, n_merger_star, n_merger_bh, last_mass_change_time, last_merger_time);
    }

    //! read class data with ASCII format
    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lld %lld %lld %lf %lf ", 
                          &type, &n_merger_star, &n_merger_bh, &last_mass_change_time, &last_merger_time);
        if(rcount<5) {
            std::cerr<<"Error: Data reading fails! requiring data number is 5, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! for print in one line
    void print(std::ostream & fout) const{
        fout<<" type= "<<type
            <<" n_merger_star= "<<n_merger_star
            <<" n_merger_bh= "<<n_merger_bh
            <<" last_mass_change_time= "<<last_mass_change_time
            <<" last_merger_time= "<<last_merger_time;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"type"
             <<std::setw(_width)<<"n_merger_star"
             <<std::setw(_width)<<"n_merger_bh"
             <<std::setw(_width)<<"last_mass_change_time"
             <<std::setw(_width)<<"last_merger_time";
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<type
             <<std::setw(_width)<<n_merger_star
             <<std::setw(_width)<<n_merger_bh
             <<std::setw(_width)<<last_mass_change_time
             <<std::setw(_width)<<last_merger_time;
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        _fout<<std::setw(_offset)<<" "<<++counter<<". type: type of object; 0: SMBH; 1: stellar-mass BH; 2: star; 3: star seed; 4: star remnant; 5: BH remnant\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". n_merger_star: times of merger with star\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". n_merger_bh: times of merger with black hole\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". last_mass_change_time: last mass change time\n";
        _fout<<std::setw(_offset)<<" "<<++counter<<". last_merger_time: last merger time\n";
        return counter;
    }
    
};
