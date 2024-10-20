#pragma once
#include "Common/binary_tree.h"
#include <iostream>

//! orbital sampling method for counter force from binary
class OrbitalSamplingManager{
    PS::S32 n_split_;    // oribital particle splitting number
    PS::F64* decca_list_;   // d ecca per orbital particle 
    PS::F64* dsin_ecca_list_;  // d sin(ecca) per orbital particle
    PS::F64  decca_; // ecca interval

    inline PS::F64 calcMeanAnomaly(const PS::F64 _ecca, const PS::F64 _sin_ecca, const PS::F64 _ecc) {
        return _ecca - _ecc*_sin_ecca;
    }
    
public:
    PS::F64 gravitational_constant; ///> gravitational constant
    
    //! initializer
    OrbitalSamplingManager(): decca_list_(NULL), dsin_ecca_list_(NULL), decca_(0.0), gravitational_constant(-1.0) {}

    //! check paramters
    bool checkParams() {
        ASSERT(n_split_>=0);
        ASSERT(gravitational_constant>0);
        if (n_split_>0) {
            ASSERT(decca_list_!=NULL);
            ASSERT(dsin_ecca_list_!=NULL);
            ASSERT(decca_>0.0);
        }
        return true;
    }

    //! create orbit sampling particles
    /*! each component have n_split_ samples with equal interval of eccentric anomaly.
        Mass is weighted by mean anomaly and summation is the same as binary mass.
        @param[in] _ptcl_artificial: particle array to store the sample particles, 2*n_split_ will be used
        @param[in] _bin: binary orbit 
     */
    template <class Tptcl>
    void createSampleParticles(Tptcl* _ptcl_artificial,
                                    COMM::BinaryTree<Tptcl> &_bin) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        PS::F64 m_check[2]={0.0,0.0};
#endif
        PS::F64 inverse_twopi = 1.0/(decca_*n_split_);
        for (int i=0; i<n_split_; i++) {
            // center_of_mass_shift(*(Tptcl*)&_bin,p,2);
            // generate particles at different orbitial phase
            COMM::Binary::orbitToParticle(_ptcl_artificial[2*i], _ptcl_artificial[2*i+1], _bin, decca_*i, gravitational_constant);

//#ifdef ARTIFICIAL_PARTICLE_DEBUG
//            PS::F64 semi_i,ecc_i,r_i,rv_i;
//            COMM::Binary::particleToSemiEcc(semi_i, ecc_i, r_i, rv_i, _ptcl_artificial[2*i], _ptcl_artificial[2*i+1], gravitational_constant);
//            if (abs(semi_i-_bin.semi)>1e-10) {
//                std::cerr<<"semi_i "<<semi_i<<" bin.semi "<<_bin.semi<<std::endl;
//                abort();
//            }
//            assert(abs(ecc_i-_bin.ecc)<1e-10);
//#endif

            //// use velocity to weight mass (not accurate)
            //PS::F64vec dvvec= _ptcl_artificial[2*i].vel - _ptcl_artificial[2*i+1].vel;
            //PS::F64 odv = 1.0/std::sqrt(dvvec*dvvec);

            PS::F64 mass_member[2]= {_bin.m1, _bin.m2};
            for (int j=0; j<2; j++) {
                Tptcl* pj = &_ptcl_artificial[2*i+j];

                // set mass
                PS::F64 dmean_anomaly = calcMeanAnomaly(decca_list_[i], dsin_ecca_list_[i], _bin.ecc);
                pj->mass = mass_member[j] * dmean_anomaly * inverse_twopi;

                // center_of_mass_correction 
                pj->pos += _bin.pos;
                pj->vel += _bin.vel;
#ifdef PETAR_USE_MPFRC
                pj->pos_mp = pj->pos;
#endif
#ifdef ARTIFICIAL_PARTICLE_DEBUG
                assert(mass_member[j]>0);
                assert(pj->mass >0);
                m_check[j] += pj->mass;
#endif
            }
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            PS::F64vec pos_i_cm = (_ptcl_artificial[2*i].mass*_ptcl_artificial[2*i].pos+_ptcl_artificial[2*i+1].mass*_ptcl_artificial[2*i+1].pos)/(_ptcl_artificial[2*i].mass + _ptcl_artificial[2*i+1].mass);
            assert(abs((pos_i_cm - _bin.pos)*(pos_i_cm - _bin.pos))<1e-10);
#endif
            //mnormal += odv;
        }

        // normalized the mass of each particles to keep the total mass the same as c.m. mass
        //PS::F64 mfactor = 1.0/mnormal;
        //for (int i=i_start_orb; i<n_pair; i++) 
        //    for (int j=0; j<2; j++) 
        //        _ptcl_artificial[2*i+j].mass *= mfactor;


#ifdef ARTIFICIAL_PARTICLE_DEBUG
        if(n_split_>0) {
            assert(abs(m_check[0]-_bin.m1)<1e-10);
            assert(abs(m_check[1]-_bin.m2)<1e-10);
        }
#endif
    }

    //! set particle split number
    /*! @param[in] _n_split: particle split number of orbital samples, total artificial particle number is 2*_n_split+TidalTensor::n_point+1, first TidalTensor::n_point are tidal tensor particles.
     */
    void setParticleSplitN(const PS::S32& _n_split) {
        ASSERT(_n_split>=0);
        n_split_ =  _n_split;

        // get interval of ecc anomaly
        if (n_split_>0) {
            decca_ = 8.0*atan(1.0)/n_split_;

            // initial array
            if (decca_list_    !=NULL) delete [] decca_list_;
            if (dsin_ecca_list_!=NULL) delete [] dsin_ecca_list_;

            decca_list_     = new PS::F64[n_split_];
            dsin_ecca_list_ = new PS::F64[n_split_];

            // calculate ecca boundary 
            PS::F64 ecca_list    [n_split_+1];
            PS::F64 sin_ecca_list[n_split_+1];
            for (PS::S32 i=0; i<=n_split_; i++) {
                ecca_list[i]     = decca_*i-0.5*decca_;
                sin_ecca_list[i] = std::sin(ecca_list[i]);
            }
            // get ecca interval per orbital particle
            for (PS::S32 i=0; i<n_split_; i++) {
                decca_list_[i]     = ecca_list    [i+1] - ecca_list    [i];
                dsin_ecca_list_[i] = sin_ecca_list[i+1] - sin_ecca_list[i];
            }
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            PS::F64 decca_sum=0.0;
            //std::cerr<<"dEcca:";
            for (PS::S32 i=0; i<n_split_; i++) {
                //std::cerr<<" "<<decca_list_[i]-dsin_ecca_list_[i];
                decca_sum += decca_list_[i] - 0.5*dsin_ecca_list_[i];
                assert(decca_list_[i]-dsin_ecca_list_[i]>=0.0);
            }
            ASSERT(abs(decca_sum-decca_*n_split_)<1e-10);
#endif            
        } 
    }

    //! get particle split number
    PS::S32 getParticleSplitN() const{
        return n_split_;
    }

    //! get particle number 
    PS::S32 getParticleN() const {
        return n_split_*2;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        fwrite(&n_split_,       sizeof(PS::S32), 1, _fp);
        fwrite(&gravitational_constant, sizeof(PS::F64), 1, _fp);
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        PS::S32 n_split;
        size_t rcount = fread(&n_split,  sizeof(PS::S32), 1, _fin);
        rcount += fread(&gravitational_constant, sizeof(PS::F64), 1, _fin);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        setParticleSplitN(n_split);
    }    

    //! destructor 
    ~OrbitalSamplingManager() {
        if (decca_list_!=NULL) {
            delete [] decca_list_;
            decca_list_=NULL;
        }
        if (dsin_ecca_list_!=NULL) {
            delete [] dsin_ecca_list_;
            dsin_ecca_list_=NULL;
        }
    }

    //! operator = 
    /*! Copy function will remove the local data and also copy the particle data or the link
     */
    OrbitalSamplingManager& operator = (const OrbitalSamplingManager& _manager) {
        if (decca_list_!=NULL) {
            delete [] decca_list_;
            decca_list_=NULL;
        }
        if (dsin_ecca_list_!=NULL) {
            delete [] dsin_ecca_list_;
            dsin_ecca_list_=NULL;
        }
        setParticleSplitN(_manager.n_split_);
        gravitational_constant = _manager.gravitational_constant;
        return *this;
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"n_split   : "<<n_split_<<std::endl
             <<"G:        : "<<gravitational_constant<<std::endl;
    }
};
