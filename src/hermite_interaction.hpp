#pragma once

#include "Common/Float.h"
#include "changeover.hpp"

//! hermite interaction class 
class HermiteInteraction{
public:
    Float eps_sq; // softening parameter
    Float gravitational_constant;      // gravitational constant

    // constructor
    HermiteInteraction(): eps_sq(Float(-1.0)), gravitational_constant(Float(-1.0)) {}

    //! check whether parameters values are correct
    /*! \return true: all correct
     */
    bool checkParams() {
        ASSERT(eps_sq>=0.0);
        ASSERT(gravitational_constant>0.0);
        return true;
    }        
    
    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"eps_sq: "<<eps_sq<<std::endl
             <<"G     : "<<gravitational_constant<<std::endl;
    }    

    //! calculate separation square between i and j particles
    /*!
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tpi, class Tpj>
    inline Float calcR2Pair(const Tpi& _pi,
                            const Tpj& _pj) {
        const Float dr[3] = {_pj.pos[0]-_pi.pos[0], 
                             _pj.pos[1]-_pi.pos[1],
                             _pj.pos[2]-_pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        return dr2;
    }

    //! calculate acceleration and jerk of one pair
    /*!
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tpi, class Tpj>
    inline Float calcAccJerkPairSingleSingle(H4::ForceH4& _fi,
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
        const Float drdot = drdv*rinv;
        const Float kp = ChangeOver::calcPotWTwo(_pi.changeover,_pj.changeover, r);
        const Float k = ChangeOver::calcAcc0WTwo(_pi.changeover, _pj.changeover, r);
        const Float kdot = ChangeOver::calcAcc1WTwo(_pi.changeover, _pj.changeover, r, drdot);
          
        const Float rinv2 = rinv*rinv;

        const Float gmor = gravitational_constant*_pj.mass*rinv;
        const Float gmor3 = gmor*rinv2; 
        const Float gmor3k = gmor3*k;
        const Float gmor3kd = gmor3*kdot;
        const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
        const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                               gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                               gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        _fi.pot += - gmor*kp;

        return dr2;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      \return the minimum distance square of i and member j
     */
    template<class Tpi, class Tgroup>
    inline Float calcAccJerkPairSingleGroupMember(H4::ForceH4& _fi,
                                                  const Tpi& _pi,
                                                  const Tgroup& _gj) {
        const int n_member = _gj.particles.getSize();
        auto* member_adr = _gj.particles.getOriginAddressArray();
        Float r2_min = NUMERIC_FLOAT_MAX;
        for (int i=0; i<n_member; i++) {
            const auto& pj = *member_adr[i];
            ASSERT(_pi.id!=pj.id);
            const Float dr[3] = {pj.pos[0]-_pi.pos[0], 
                                 pj.pos[1]-_pi.pos[1],
                                 pj.pos[2]-_pi.pos[2]};
            Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float dr2_eps = dr2 + eps_sq;
            const Float dv[3] = {pj.vel[0] - _pi.vel[0],
                                 pj.vel[1] - _pi.vel[1],
                                 pj.vel[2] - _pi.vel[2]};
            const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
            const Float r = sqrt(dr2_eps);
            ASSERT(r>0.0);
            const Float rinv = 1.0/r;
            const Float drdot = drdv*rinv;
            const Float kp = ChangeOver::calcPotWTwo(_pi.changeover, pj.changeover, r);
            const Float k = ChangeOver::calcAcc0WTwo(_pi.changeover, pj.changeover, r);
            const Float kdot = ChangeOver::calcAcc1WTwo(_pi.changeover, pj.changeover, r, drdot);
          
            const Float rinv2 = rinv*rinv;

            const Float gmor = gravitational_constant*pj.mass*rinv;
            const Float gmor3 = gmor*rinv2; 
            const Float gmor3k = gmor3*k;
            const Float gmor3kd = gmor3*kdot;
            const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
            const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                                   gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                                   gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
            _fi.acc0[0] += acc0[0];
            _fi.acc0[1] += acc0[1];
            _fi.acc0[2] += acc0[2];

            _fi.acc1[0] += acc1[0];
            _fi.acc1[1] += acc1[1];
            _fi.acc1[2] += acc1[2];

            _fi.pot += - gmor*kp;

            r2_min = std::min(r2_min, dr2);
        }
        return r2_min;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      @param[in]: _pj: predicted group j cm
      \return the distance square of the pair
     */
    template<class Tpi, class Tgroup, class Tpcmj>
    inline Float calcAccJerkPairSingleGroupCM(H4::ForceH4& _fi,
                                              const Tpi& _pi,
                                              const Tgroup& _gj,
                                              const Tpcmj& _pj) {
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
        const Float drdot = drdv*rinv;

        Float mkp = 0.0;
        Float mk = 0.0;
        Float mkdot = 0.0;
        auto* ptcl_mem = _gj.particles.getDataAddress();
        ASSERT(_gj.particles.cm.mass==_pj.mass);
        for (int i=0; i<_gj.particles.getSize(); i++) {
            mkp   += ptcl_mem[i].mass * ChangeOver::calcPotWTwo(_pi.changeover, ptcl_mem[i].changeover, r);
            mk    += ptcl_mem[i].mass * ChangeOver::calcAcc0WTwo(_pi.changeover, ptcl_mem[i].changeover, r);
            mkdot += ptcl_mem[i].mass * ChangeOver::calcAcc1WTwo(_pi.changeover, ptcl_mem[i].changeover, r, drdot);
        }
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float gmor3k = gravitational_constant*mk*rinv3;
        const Float gmor3kd = gravitational_constant*mkdot*rinv3;
        const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
        const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                               gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                               gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        _fi.pot += - gravitational_constant*mkp*rinv;

        return dr2;
    }

    //! calculate acceleration and jerk of one pair group cm and single
    /*!
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _pj: particle j
      \return the distance square of the pair
     */
    template<class Tgroup, class Tpcmi, class Tpj>
    inline Float calcAccJerkPairGroupCMSingle(H4::ForceH4& _fi,
                                              const Tgroup& _gi,
                                              const Tpcmi& _pi,
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
        const Float drdot = drdv*rinv;

        Float kp = 0.0;
        Float k = 0.0;
        Float kdot = 0.0;
        auto* ptcl_mem = _gi.particles.getDataAddress();
        ASSERT(_gi.particles.cm.mass==_pi.mass);
        for (int i=0; i<_gi.particles.getSize(); i++) {
            kp   += ptcl_mem[i].mass * ChangeOver::calcPotWTwo(_pi.changeover, ptcl_mem[i].changeover, r);
            k    += ptcl_mem[i].mass * ChangeOver::calcAcc0WTwo(_pj.changeover, ptcl_mem[i].changeover, r);
            kdot += ptcl_mem[i].mass * ChangeOver::calcAcc1WTwo(_pj.changeover, ptcl_mem[i].changeover, r, drdot);
        }
        kp   /= _pi.mass;
        k    /= _pi.mass;
        kdot /= _pi.mass;

        const Float rinv2 = rinv*rinv;

        const Float gmor = gravitational_constant*_pj.mass*rinv;
        const Float gmor3 = gmor*rinv2; 
        const Float gmor3k = gmor3*k;
        const Float gmor3kd = gmor3*kdot;
        const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
        const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                               gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                               gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        _fi.pot += - gmor*kp;

        return dr2;
    }


    //! calculate acceleration and jerk of one pair group cm and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      \return the minimum distance square of i and member j
     */
    template<class Tpi, class Tgroup>
    inline Float calcAccJerkPairGroupCMGroupMember(H4::ForceH4& _fi,
                                                   const Tgroup& _gi,
                                                   const Tpi& _pi,
                                                   const Tgroup& _gj) {
        const int n_member = _gj.particles.getSize();
        auto* member_adr = _gj.particles.getOriginAddressArray();
        Float r2_min = NUMERIC_FLOAT_MAX;
        ASSERT(_gi.particles.cm.mass==_pi.mass);
        for (int i=0; i<n_member; i++) {
            const auto& pj = *member_adr[i];
            ASSERT(_pi.id!=pj.id);
            const Float dr[3] = {pj.pos[0]-_pi.pos[0], 
                                 pj.pos[1]-_pi.pos[1],
                                 pj.pos[2]-_pi.pos[2]};
            Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
            Float dr2_eps = dr2 + eps_sq;
            const Float dv[3] = {pj.vel[0] - _pi.vel[0],
                                 pj.vel[1] - _pi.vel[1],
                                 pj.vel[2] - _pi.vel[2]};
            const Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];
            const Float r = sqrt(dr2_eps);
            ASSERT(r>0.0);
            const Float rinv = 1.0/r;
            const Float drdot = drdv*rinv;

            Float k = 0.0;
            Float kdot = 0.0;
            Float kp = 0.0;
            auto* ptcl_mem = _gi.particles.getDataAddress();
            for (int i=0; i<_gi.particles.getSize(); i++) {
                kp   += ptcl_mem[i].mass * ChangeOver::calcPotWTwo(_pi.changeover, ptcl_mem[i].changeover, r);
                k    += ptcl_mem[i].mass * ChangeOver::calcAcc0WTwo(pj.changeover, ptcl_mem[i].changeover, r);
                kdot += ptcl_mem[i].mass * ChangeOver::calcAcc1WTwo(pj.changeover, ptcl_mem[i].changeover, r, drdot);
            }
            kp   /= _pi.mass;
            k    /= _pi.mass;
            kdot /= _pi.mass;

            const Float rinv2 = rinv*rinv;

            const Float gmor = gravitational_constant*pj.mass*rinv;
            const Float gmor3 = gmor*rinv2; 
            const Float gmor3k = gmor3*k;
            const Float gmor3kd = gmor3*kdot;
            const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
            const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                                   gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                                   gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
            _fi.acc0[0] += acc0[0];
            _fi.acc0[1] += acc0[1];
            _fi.acc0[2] += acc0[2];

            _fi.acc1[0] += acc1[0];
            _fi.acc1[1] += acc1[1];
            _fi.acc1[2] += acc1[2];

            _fi.pot += - gmor*kp;

            r2_min = std::min(r2_min, dr2);
        }
        return r2_min;
    }

    //! calculate acceleration and jerk of one pair single and resolved group
    /*! 
      @param[out]: _fi: acceleration for i particle
      @param[in]: _gi: particle group i
      @param[in]: _pi: particle i
      @param[in]: _gj: particle group j
      @param[in]: _pj: predicted group j cm
      \return the distance square of the pair
     */
    template<class Tpcmi, class Tgroup, class Tpcmj>
    inline Float calcAccJerkPairGroupCMGroupCM(H4::ForceH4& _fi,
                                               const Tgroup& _gi,
                                               const Tpcmi& _pi,
                                               const Tgroup& _gj,
                                               const Tpcmj& _pj) {
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
        const Float drdot = drdv*rinv;

        // pj.mass * k
        Float mk = 0.0;
        Float mkdot = 0.0;
        Float mkp = 0.0;
        auto* ptcl_mem_i = _gi.particles.getDataAddress();
        auto* ptcl_mem_j = _gj.particles.getDataAddress();
        for (int i=0; i<_gi.particles.getSize(); i++) {
            Float mkj = 0.0;
            Float mkpj = 0.0;
            Float mkdotj = 0.0;
            for (int j=0; j<_gj.particles.getSize(); j++) {
                mkpj   += ptcl_mem_j[j].mass * ChangeOver::calcPotWTwo(ptcl_mem_i[i].changeover, ptcl_mem_j[j].changeover, r);
                mkj    += ptcl_mem_j[j].mass * ChangeOver::calcAcc0WTwo(ptcl_mem_i[i].changeover, ptcl_mem_j[j].changeover, r);
                mkdotj += ptcl_mem_j[j].mass * ChangeOver::calcAcc1WTwo(ptcl_mem_i[i].changeover, ptcl_mem_j[j].changeover, r, drdot);
            }
            mkp   += ptcl_mem_i[i].mass * mkpj;
            mk    += ptcl_mem_i[i].mass * mkj;
            mkdot += ptcl_mem_i[i].mass * mkdotj;
        }
        mkp   /= _pi.mass;
        mk    /= _pi.mass;
        mkdot /= _pi.mass;
          
        const Float rinv2 = rinv*rinv;
        const Float rinv3 = rinv2*rinv;

        const Float gmor3k = gravitational_constant*mk*rinv3;
        const Float gmor3kd = gravitational_constant*mkdot*rinv3;
        const Float acc0[3] = {gmor3k*dr[0], gmor3k*dr[1], gmor3k*dr[2]};
        const Float acc1[3] = {gmor3k*dv[0] - 3.0*drdv*rinv2*acc0[0] + gmor3kd*dr[0],
                               gmor3k*dv[1] - 3.0*drdv*rinv2*acc0[1] + gmor3kd*dr[1],
                               gmor3k*dv[2] - 3.0*drdv*rinv2*acc0[2] + gmor3kd*dr[2]};
        _fi.acc0[0] += acc0[0];
        _fi.acc0[1] += acc0[1];
        _fi.acc0[2] += acc0[2];

        _fi.acc1[0] += acc1[0];
        _fi.acc1[1] += acc1[1];
        _fi.acc1[2] += acc1[2];

        _fi.pot += - gravitational_constant*mkp*rinv;

        return dr2;
    }

    //! calculate pair potential energy
    template<class Tpi, class Tpj>
    Float calcEnergyPotSingleSingle(const Tpi& pi, const Tpj& pj) {
        const Float dr[3] = {pj.pos[0] - pi.pos[0], 
                             pj.pos[1] - pi.pos[1],
                             pj.pos[2] - pi.pos[2]};
        Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float dr2_eps = dr2 + eps_sq;
        const Float r = sqrt(dr2_eps);
        ASSERT(r>0.0);
        const Float rinv = 1.0/r;
        const Float k = ChangeOver::calcPotWTwo(pi.changeover, pj.changeover, r);
        
        return -gravitational_constant*pi.mass*pj.mass*rinv*k;
    }

    //! calculate kinetic and potential energy of the system
    /*!
      @param[out] _energy: hermite energy
      @param[in] _particles: (all) particle list
      @param[in] _n_particle: number of particles
      @param[in] _groups: group list
      @param[in] _group_index: active groups index array for energy calculation
      @param[in] _n_group: number of active groups
      @param[in] _perturber: perturber 
     */
    template<class Tp, class Tgroup, class Tpert>
    inline void calcEnergy(H4::HermiteEnergy& _energy, const Tp* _particles, const int _n_particle, const Tgroup* _groups, const int* _group_index, const int _n_group, const Tpert& _perturber) {
        _energy.ekin = _energy.epot = _energy.epert = 0.0;
        for (int i=0; i<_n_particle; i++) {
            auto& pi = _particles[i];
            if (pi.mass==0.0) continue;
            _energy.ekin += pi.mass* (pi.vel[0]*pi.vel[0] + pi.vel[1]*pi.vel[1] + pi.vel[2]*pi.vel[2]);
            Float poti = 0.0;
            for (int j=0; j<i; j++) {
                if (i==j) continue;
                auto& pj = _particles[j];
                if (pj.mass==0.0) continue;
                const Float dr[3] = {pj.pos[0] - pi.pos[0], 
                                     pj.pos[1] - pi.pos[1],
                                     pj.pos[2] - pi.pos[2]};
                Float dr2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float dr2_eps = dr2 + eps_sq;
                const Float r = sqrt(dr2_eps);
                ASSERT(r>0.0);
                const Float rinv = 1.0/r;
                const Float k = ChangeOver::calcPotWTwo(pi.changeover, pj.changeover, r);
        
                poti += -pj.mass*rinv*k;
            }
            _energy.epot += gravitational_constant*poti*pi.mass;
        }

#ifdef SOFT_PERT
        // tidal tensor energy
        for (int k=0; k<_n_group; k++) {
            const int i = _group_index[k];
            const int n_member = _groups[i].particles.getSize();
            if (_groups[i].perturber.soft_pert!=NULL) {
                auto* pert = _groups[i].perturber.soft_pert;
                for (int j=0; j<n_member; j++) {
                    _energy.epert += _groups[i].particles[j].mass*pert->evalPot(_groups[i].particles[j].pos);
                }
            }
        }
#endif
        _energy.ekin *= 0.5;
        //_energy.epert *= 0.5;
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
