#pragma once

#include "Common/Float.h"
#include "changeover.hpp"

//! hermite interaction class 
class HermiteInteraction{
public:
    Float eps_sq; // softening parameter
    Float G;      // gravitational constant
    ChangeOver* changeover; ///> changover control

    // constructor
    HermiteInteraction(): eps_sq(Float(-1.0)), G(Float(-1.0)) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(G>0.0);
        ASSERT(changeover!=NULL);
        return true;
    }        
    
    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq: "<<eps_sq<<std::endl
             <<"G     : "<<G<<std::endl;
    }    

    //! calculate acceleration and jerk of one pair
    /*! \return the distance square of the pair
     */
    template<class Tpi, class Tpj>
    inline Float calcAccJerkPair(H4::ForceH4& _fi,
                                 const Tpi& _pi,
                                 const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float dv[3] = {_pj.vel[0] - _pi.vel[0],
                             _pj.vel[1] - _pi.vel[1],
                             _pj.vel[2] - _pi.vel[2]};
        const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
        const Float r = sqrt(dr2_eps);
        ASSERT(r>0.0);
        const Float rinv = 1.0/r;
        const PS::F64 k = changeover->calcAcc0W(r);
        const PS::F64 kdot = changeover->calcAcc1W(r);
          
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float mor3 = G*_pj.mass*rinv3; 
        const Float mor3k = mor3*k;
        const Float mor3kd = mor3*drdv*kdot;
        const Float acc0[3] = {mor3k*dr[0], mor3k*dr[1], mor3k*dr[2]};
        const Float acc1[3] = {mor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + mor3kd*dr[0],
                               mor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + mor3kd*dr[1],
                               mor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + mor3kd*dr[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        return dr2;
    }

    //! calculate potential from j to i
    /*! \return the potential
     */
    template<class Tpi, class Tpj>
    inline Float calcPotPair(const Tpi& _pi,
                             const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float r = sqrt(dr2_eps);
        ASSERT(r>0.0);
        const Float rinv = 1.0/r;
        const PS::F64 k = changeover->calcPotW(r);
        
        return -G*_pj.mass*rinv*k;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
        fwrite(this, sizeof(*this),1,_fp);
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
};
