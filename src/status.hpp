#pragma once

#include "particle_base.hpp"
#include "energy.hpp"

//! class for measure the status of the system
class Status {
public:
    PS::F64 time;  // time of system
    PS::S64 n_real_loc;  // number of real particles in the local MPI process
    PS::S64 n_real_glb;  // number of real particles in all MPI processes
    PS::S64 n_all_loc;   // number of all (real+aritficial) particles in the local MPI process
    PS::S64 n_all_glb;   // number of all (real+aritficial) particles in all MPI process
    PS::S32 n_remove_glb; // number of removed particles in all MPI processes
    PS::S32 n_escape_glb; // number of escaped particles in all MPI processes
    PS::F64 half_mass_radius; // half mass radius of the system (not calculated)
    EnergyAndMomemtum energy; // energy of the system
    struct ParticleCM{ // local structure for system center
        PS::F64 mass;
        PS::F64vec pos;
        PS::F64vec vel;
#ifdef SMOOTH_CM_USING_VEL_PRED
        PS::F64vec acc;
        PS::F64 time_record;
#endif
#ifdef SMOOTH_CM_USING_RECORES
        std::vector<PS::F64> pos_bk[3], vel_bk[3], time_bk;
#endif
        bool is_center_shift_flag;

        ParticleCM(): mass(0.0), pos(PS::F64vec(0.0)), vel(PS::F64vec(0.0)), 
#ifdef SMOOTH_CM_USING_VEL_PRED
                      acc(PS::F64vec(0.0)), time_record(0.0), 
#endif
#ifdef SMOOTH_CM_USING_RECORES
                      pos_bk{}, vel_bk{}, time_bk(), 
#endif
                      is_center_shift_flag(false) {}

        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumnTitle(std::ofstream & _fout, const PS::S32 _width=20) {
            _fout<<std::setw(_width)<<"CM.mass"
                 <<std::setw(_width)<<"CM.pos.x"
                 <<std::setw(_width)<<"CM.pos.y"
                 <<std::setw(_width)<<"CM.pos.z"
                 <<std::setw(_width)<<"CM.vel.x"
                 <<std::setw(_width)<<"CM.vel.y"
                 <<std::setw(_width)<<"CM.vel.z";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ofstream & _fout, const PS::S32 _width=20) {
            _fout<<std::setw(_width)<<mass
                 <<std::setw(_width)<<pos.x
                 <<std::setw(_width)<<pos.y
                 <<std::setw(_width)<<pos.z
                 <<std::setw(_width)<<vel.x
                 <<std::setw(_width)<<vel.y
                 <<std::setw(_width)<<vel.z;
        }

        //! print title and values in one lines
        /*! print titles and values in one lines
          @param[out] _fout: std::ostream output object
        */
        void print(std::ostream & _fout) {
            _fout<<"C.M.: mass: "<<mass
                 <<" pos: "<<pos
                 <<" vel: "<<vel;
        }

        void clear() {
            is_center_shift_flag = false;
            mass = 0.0;
            pos = PS::F64vec(0.0);
            vel = PS::F64vec(0.0);
#ifdef SMOOTH_CM_USING_VEL_PRED
            acc = PS::F64vec(0.0);
            time_record = 0.0;
#endif
#ifdef SMOOTH_CM_USING_RECORES
            for (int k=0; k<3; k++) {
                pos_bk[k].clear();
                vel_bk[k].clear();
            }
            time_bk.clear();
#endif
        }

    } pcm;

    Status(): time(0.0), n_real_loc(0), n_real_glb(0), n_all_loc(0), n_all_glb(0), n_remove_glb(0), n_escape_glb(0), half_mass_radius(0), energy(), pcm() {}

    //! shift particle system center to c.m. frame 
    /*! set is_center_shift_flag to true
      @param[in,out] _tsys: particle system
      @param[in] _n: number of particle
    */
    template <class Tsoft>
    void shiftToCenterOfMassFrame(Tsoft* _tsys, const PS::S64 _n) {
        if (!pcm.is_center_shift_flag) {
            for (int i=0; i<_n; i++) {
                _tsys[i].pos -= pcm.pos;
                _tsys[i].vel -= pcm.vel;
            }
        }
        pcm.is_center_shift_flag = true;
    }

    //! shift particle system center back to original frame
    /*! set is_center_shift_flag to false
      @param[in,out] _tsys: particle system
      @param[in] _n: number of particle
    */
    template <class Tsoft>
    void shiftToOriginFrame(Tsoft* _tsys, const PS::S64 _n) {
        if (pcm.is_center_shift_flag) {
            for (int i=0; i<_n; i++) {
                _tsys[i].pos += pcm.pos;
                _tsys[i].vel += pcm.vel;
            }
        }
        pcm.is_center_shift_flag = false;
    }

    //! calculate the center of system 
    /*!
      @param[in] _tsys: particle system
      @param[in] _n: number of particle
      @param[in] _mode: calculation mode, 1: center-of-the-mass; 2: number (no mass) weighted center; 3: soft potential weighted center
     */
    template <class Tsoft>
    void calcCenterOfMass(Tsoft* _tsys, const PS::S64 _n, int _mode=3) {
        PS::F64 mass = 0.0;
        PS::F64vec pos_cm = PS::F64vec(0.0);
        PS::F64vec vel_cm = PS::F64vec(0.0);
#ifdef SMOOTH_CM_USING_VEL_PRED
        PS::F64vec acc_cm = PS::F64vec(0.0);
#endif

        if (_mode==1) { // center of the mass
//#pragma omp declare reduction(+:PS::F64vec:omp_out += omp_in) initializer (omp_priv=PS::F64vec(0.0))
//#pragma omp parallel for reduction(+:mass,pos_cm,vel_cm)
            for (int i=0; i<_n; i++) {
                auto& pi = _tsys[i];
                PS::F64 mi = pi.mass;
                mass  += mi;
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(pi.vel.x));
                assert(!std::isnan(pi.vel.y));
                assert(!std::isnan(pi.vel.z));
#endif
                pos_cm += mi*pi.pos;
                vel_cm += mi*pi.vel;
#ifdef SMOOTH_CM_USING_VEL_PRED
                acc_cm += mi*pi.acc;
#endif
            }        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pcm.mass = PS::Comm::getSum(mass);
            pcm.pos  = PS::Comm::getSum(pos_cm);
            pcm.vel  = PS::Comm::getSum(vel_cm);
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = PS::Comm::getSum(acc_cm);
#endif
#else
            pcm.mass = mass;
            pcm.pos  = pos_cm;
            pcm.vel  = vel_cm;
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = acc_cm;
#endif
#endif
            if (pcm.mass>0) {
                pcm.pos /= pcm.mass;
                pcm.vel /= pcm.mass;
#ifdef SMOOTH_CM_USING_VEL_PRED
                pcm.acc /= pcm.mass;
#endif
            }
        }
        else if (_mode==2) { // no mass weighted center
//#pragma omp declare reduction(+:PS::F64vec:omp_out += omp_in) initializer (omp_priv=PS::F64vec(0.0))
//#pragma omp parallel for reduction(+:mass,pos_cm,vel_cm)
            for (int i=0; i<_n; i++) {
                auto& pi = _tsys[i];
                PS::F64 mi = pi.mass;
                mass  += mi;
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(pi.vel.x));
                assert(!std::isnan(pi.vel.y));
                assert(!std::isnan(pi.vel.z));
#endif
                pos_cm += pi.pos;
                vel_cm += pi.vel;
#ifdef SMOOTH_CM_USING_VEL_PRED
                acc_cm += pi.acc;
#endif
            }        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pcm.mass = PS::Comm::getSum(mass);
            pcm.pos  = PS::Comm::getSum(pos_cm);
            pcm.vel  = PS::Comm::getSum(vel_cm);
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = PS::Comm::getSum(acc_cm);
#endif
            PS::S64 n_glb = PS::Comm::getSum(_n);
#else
            pcm.mass = mass;
            pcm.pos  = pos_cm;
            pcm.vel  = vel_cm;
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = acc_cm;
#endif
            PS::S64 n_glb = _n;
#endif
            if (n_glb>0) {
                pcm.pos /= PS::F64(n_glb);
                pcm.vel /= PS::F64(n_glb);
#ifdef SMOOTH_CM_USING_VEL_PRED
                pcm.acc /= PS::F64(n_glb);
#endif
            }
        }
        else if (_mode==3) { // soft potential
            PS::F64 pot_tot = 0.0;
            for (int i=0; i<_n; i++) {
                auto& pi = _tsys[i];
                PS::F64 mi = pi.mass;
                PS::F64 poti = pi.pot_soft;
#ifdef EXTERNAL_POT_IN_PTCL
                poti -= pi.pot_ext; // remove external potential
#endif
                mass  += mi;
#ifdef NAN_CHECK_DEBUG
                assert(!std::isnan(pi.vel.x));
                assert(!std::isnan(pi.vel.y));
                assert(!std::isnan(pi.vel.z));
#endif
                pos_cm += poti*pi.pos;
                vel_cm += poti*pi.vel;
#ifdef SMOOTH_CM_USING_VEL_PRED
                acc_cm += poti*pi.acc;
#endif
                pot_tot += poti;
            }        

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pcm.mass = PS::Comm::getSum(mass);
            pcm.pos  = PS::Comm::getSum(pos_cm);
            pcm.vel  = PS::Comm::getSum(vel_cm);
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = PS::Comm::getSum(acc_cm);
#endif
            PS::F64 pot_tot_glb = PS::Comm::getSum(pot_tot);
#else
            pcm.mass = mass;
            pcm.pos  = pos_cm;
            pcm.vel  = vel_cm;
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.acc  = acc_cm;
#endif
            PS::F64 pot_tot_glb = pot_tot;
#endif
            if (pot_tot_glb!=0) {
                pcm.pos /= pot_tot_glb;
                pcm.vel /= pot_tot_glb;
#ifdef SMOOTH_CM_USING_VEL_PRED
                pcm.acc /= pot_tot_glb;
#endif
            }
        }
    }

    //! calculate the center of system and shift particle systems to center frame
    /*!
      @param[in] _tsys: particle system
      @param[in] _n: number of particle
      @param[in] _mode: calculation mode, 1: center-of-the-mass; 2: number (no mass) weighted center; 3: soft potential weighted center (default: 3).
      @param[in] _initial_flag: if true, it is the initial calculation of c.m. without shift (default: false).
     */
    template <class Tsoft>
    void calcAndShiftCenterOfMass(Tsoft* _tsys, const PS::S64 _n, const int _mode=3, const bool _initial_flag=false) {
        if (_initial_flag) {
            if(pcm.is_center_shift_flag) {
                std::cerr<<"Error: particle system is in c.m. frame, cannot initial c.m."<<std::endl;
                abort();
            }
            pcm.clear();
            pcm.is_center_shift_flag=true;
#ifdef SMOOTH_CM_USING_VEL_PRED
            pcm.time_record = time;
#endif
        }
        assert(pcm.is_center_shift_flag);

        ParticleCM pcm_bk = pcm;
        calcCenterOfMass(_tsys, _n, _mode);
        
#ifdef SMOOTH_CM_USING_VEL_PRED
        // position jump check if current time > recored time
        PS::F64 dt = time - pcm.time_record;
        if (dt>0.0) {
            // new position prediction from backup pcm
            PS::F64vec vel_pred = pcm_bk.vel + 0.5*pcm_bk.acc*dt; // half step to get median velocity
            PS::F64vec dpos_pred = vel_pred*dt;
            PS::F64 dpos_rel_pred = std::sqrt(dpos_pred[0]*dpos_pred[0] + dpos_pred[1]*dpos_pred[2] + dpos_pred[1]*dpos_pred[2]);

#ifdef RECORD_CM_IN_HEADER
            // new calculated pcm is in old c.m. frame
            PS::F64vec dpos = pcm.pos;
#else
            PS::F64vec dpos = pcm.pos - pcm_bk.pos;
#endif
            PS::F64 dpos_rel = std::sqrt(dpos[0]*dpos[0] + dpos[1]*dpos[1] + dpos[2]*dpos[2]);
        
            if (dpos_rel > 10 * dpos_rel_pred) {
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
                PS::S32 my_rank = PS::Comm::getRank();
                if (my_rank == 0) {
                    std::cerr<<"Large jump of c.m. position detected, used predicted position instead, new pos shift: "<<pcm.pos<<"; predicted pos shift: "<<dpos_pred<<"; old pos: "<<pcm_bk.pos<<"; old vel: "<<pcm_bk.vel<<std::endl;
                }
#else
                std::cerr<<"Large jump of c.m. position detected, used predicted position instead, new pos shift: "<<pcm.pos<<"; predicted pos shift: "<<dpos_pred<<"; old pos: "<<pcm_bk.pos<<"; old vel: "<<pcm_bk.vel<<std::endl;
#endif
#ifdef RECORD_CM_IN_HEADER
                // only include differences of predicted and old positions and velocities
                pcm.pos = 10*dpos_pred;
                pcm.vel = 10*pcm_bk.acc*dt;
#else
                // record predicted position and velocity
                pcm.pos = pcm_bk.pos + 10*dpos_pred;
                pcm.vel = pcm_bk.vel + 10*pcm_bk.acc*dt;
#endif
            }
            pcm.time_record = time;
        }
#endif // SMOOTH_CM_USING_VEL_PRED

#ifdef SMOOTH_CM_USING_RECORES
        bool record_flag = true;

        // do not record if time_bk == time
        if (pcm.time_bk.size()>0) 
            record_flag = (pcm.time_bk.back() < time);
            
        if (record_flag) {
            for (int k=0; k<3; k++) {
                pcm.pos_bk[k].push_back(pcm.pos[k]);
                pcm.vel_bk[k].push_back(pcm.vel[k]);
            }
            pcm.time_bk.push_back(time);
        }

        std::cerr<<"Record size: "<<pcm.time_bk.size()<<"; time: "<<pcm.time_bk.back()
        <<"; pos: "<<pcm.pos_bk[0].back()<<" "<<pcm.pos_bk[1].back()<<" "<<pcm.pos_bk[2].back()
        <<std::endl;
        
        // when recorded 16 times, calculate the median value and update pcm
        if (pcm.time_bk.size() == 16) {
            
            for (int k=0; k<3; k++) {
                std::sort(pcm.pos_bk[k].begin(), pcm.pos_bk[k].end());
                std::sort(pcm.vel_bk[k].begin(), pcm.vel_bk[k].end());
            
                // Use the median value of sorted x, y, z
                pcm.pos[k] = 0.5*(pcm.pos_bk[k][8] + pcm.pos_bk[k][9]);
                pcm.vel[k] = 0.5*(pcm.vel_bk[k][8] + pcm.vel_bk[k][9]);

                // reset
                pcm.pos_bk[k].clear();
                pcm.vel_bk[k].clear();
            }
            pcm.time_bk.clear();
        }

        // for initialization or at 16th records, update pcm
        if (_initial_flag || pcm.time_bk.size() ==0) {
#endif
            // correct particle
            for (int i=0; i<_n; i++) {
                _tsys[i].pos -= pcm.pos;
                _tsys[i].vel -= pcm.vel;
            }
            pcm.pos += pcm_bk.pos;
            pcm.vel += pcm_bk.vel;
#ifdef SMOOTH_CM_USING_RECORES
            std::cerr<<"Update: time: "<<time<<"; pcm: "<<pcm.pos<<std::endl;
        }
        else {
            pcm.pos = pcm_bk.pos;
            pcm.vel = pcm_bk.vel;
        }
#endif
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ofstream & _fout, const PS::S32 _width=20) {
        _fout<<std::setw(_width)<<"Time"
             <<std::setw(_width)<<"N_real_loc"
             <<std::setw(_width)<<"N_real_glb"
             <<std::setw(_width)<<"N_all_loc"
             <<std::setw(_width)<<"N_all_glb"
             <<std::setw(_width)<<"N_rm_glb"
             <<std::setw(_width)<<"N_esc_glb";
        energy.printColumnTitle(_fout, _width);
        pcm.printColumnTitle(_fout, _width);
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ofstream & _fout, const PS::S32 _width=20) {
        _fout<<std::setw(_width)<<time
             <<std::setw(_width)<<n_real_loc
             <<std::setw(_width)<<n_real_glb
             <<std::setw(_width)<<n_all_loc
             <<std::setw(_width)<<n_all_glb
             <<std::setw(_width)<<n_remove_glb
             <<std::setw(_width)<<n_escape_glb;
        energy.printColumn(_fout, _width);
        pcm.printColumn(_fout, _width);
    }

    //! print title and values in one lines
    /*! print titles and values in one lines
      @param[out] _fout: std::ostream output object
      @param[in] _precision: floating data precision except time (15 digital)
    */
    void print(std::ostream & _fout, const PS::S32 _precision=7) {
        _fout<<"Time: "<<std::setprecision(15)<<time
             <<std::setprecision(_precision);
        _fout<<"  N_real(loc): "<<n_real_loc
             <<"  N_real(glb): "<<n_real_glb
             <<"  N_all(loc): "<<n_all_loc
             <<"  N_all(glb): "<<n_all_glb
             <<"  N_remove(glb): "<<n_remove_glb
             <<"  N_escape(glb): "<<n_escape_glb
             <<std::endl;
        energy.print(_fout);
        pcm.print(_fout);
        _fout<<std::endl;
    }
};
