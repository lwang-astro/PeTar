#pragma once
#include<particle_simulator.hpp>
#include"kepler.hpp"

#define ID_PHASE_SHIFT 4
#define ID_PHASE_MASKER 0xF

#ifndef ARRAY_ALLOW_LIMIT
#define ARRAY_ALLOW_LIMIT 1000000000
#endif

template<class Tptcl>
class SearchGroup{
private:
    typedef std::pair<PS::S32, PS::S32> PLinker;

    // group member list
    PS::ReallocatableArray<PS::S32> group_list_;
    PS::ReallocatableArray<PS::S32> group_list_disp_;
    PS::ReallocatableArray<PS::S32> group_list_n_;

    //! Search partner to create group
    /* If r.v >0, Use bin factor sqrt(r^2 - (r*v/|v|)^2) to check, otherwise use distance to check. The mass ratio is also considered.
       The perturber acceleration from non-partner member is stored in mass_bk for the later stablility check
       @param[out] _part_list: partner index in _ptcl for each particle
       @param[out] _part_list_disp: offset to separate partner index for different particles
       @param[out] _part_list_n: number of partners for each particle
       @param[in,out] _ptcl: particle data, mass_bk is updated
       @param[in] _n: number of particles
    */
    void searchPartner(PS::ReallocatableArray<PS::S32> & _part_list,
                       PS::ReallocatableArray<PS::S32> & _part_list_disp,
                       PS::ReallocatableArray<PS::S32> & _part_list_n,
                       Tptcl *_ptcl,
                       const PS::S32 _n) {

        _part_list.clearSize();
        _part_list_disp.reserve(_n);
#ifdef HARD_DEBUG
        assert(_n<ARRAY_ALLOW_LIMIT);
#endif        
        _part_list_disp.resizeNoInitialize(_n);
        _part_list_n.reserve(_n);
        _part_list_n.resizeNoInitialize(_n);
        
        // find partner
        PS::S32 offset = 0;
        for(int i=0; i<_n; i++) {
            //PS::S32 ip = p_list[i];
            _ptcl[i].mass_bk.d = 0.0;
            _part_list_n[i] = 0;
            _part_list_disp[i] = offset;
            for(int j=0; j<_n; j++) {
                //PS::S32 jp = p_list[j];
                if(i==j) continue;
                PS::F64vec dr = _ptcl[i].pos-_ptcl[j].pos;
                PS::F64 r2 = dr*dr;
                // use simple criterion
                PS::F64 rin_min = std::min(_ptcl[i].changeover.getRin(), _ptcl[j].changeover.getRin());
                if (r2<rin_min*rin_min) {
                    _part_list.push_back(j);
                    _part_list_n[i]++;
                    offset++;
                }
                else {
                    //store perturbation force
                    _ptcl[i].mass_bk.d += _ptcl[j].mass/r2; 
                }

                //PS::F64vec dv = _ptcl[i].vel-_ptcl[j].vel;
                //PS::F64 rv = dr*dv;
                //PS::F64 r2 = dr*dr;
                //PS::F64 mi = _ptcl[i].mass;
                //PS::F64 mj = _ptcl[j].mass;
                //PS::F64 mass_factor= mi>mj ? mi/mj : mj/mi; ///> get mass ratio
                //PS::F64 mass_factor=std::max(ptcl[ip].mass,ptcl[jp].mass)*Ptcl::mean_mass_inv;
#ifdef HARD_DEBUG
                if(r2==0) {
                    std::cerr<<"Error: zero distance! i="<<i<<" j="<<j<<std::endl;
                    abort();
                }
#endif
                // incoming case
                //if (rv>0) {
                //    PS::F64 v2 = dv*dv;
                //    //  use bin factor to check
                //    if (r2*v2-rv*rv<_r_crit2*v2*mass_factor) {
                //        _part_list.push_back(j);
                //        _part_list_n[i]++;
                //        offset++;
                //    }
                //    else {
                //        //store perturbation force
                //        _ptcl[i].mass_bk += _ptcl[j].mass/r2; 
                //    }
                //}
                //else { // outgoing case
                //    if (r2<_r_crit2*mass_factor) {
                //        _part_list.push_back(j);
                //        _part_list_n[i]++;
                //        offset++;
                //    }
                //    else {
                //        //store perturbation force
                //        _ptcl[i].mass_bk += _ptcl[j].mass/r2; 
                //    }
                //}
            }
        }
    }

    void mergeCluster(PS::ReallocatableArray<PS::S32> & group_list,
                      PS::ReallocatableArray<PS::S32> & group_list_disp,
                      PS::ReallocatableArray<PS::S32> & group_list_n,
                      //PS::ReallocatableArray<PS::S32> & p_list,
                      const PS::S32 _n_ptcl,
                      PS::S32 part_list[],
                      PS::S32 part_list_disp[],
                      PS::S32 part_list_n[]) {

        //const PS::S32 _n_ptcl = p_list.size();
        // partner index with marker
        PS::ReallocatableArray<PLinker> partner_index; 
        partner_index.reserve(_n_ptcl);

        // map index from ptcl_org to partner_index
        PS::ReallocatableArray<PS::S32> reverse_list; 
#ifdef HARD_DEBUG
        assert(_n_ptcl<ARRAY_ALLOW_LIMIT);
#endif        
        reverse_list.reserve(_n_ptcl);
        reverse_list.resizeNoInitialize(_n_ptcl);

#ifdef HARD_DEBUG
        assert(group_list.size()==0);
        assert(group_list_disp.size()==0);
        assert(group_list_n.size()==0);
#endif

        for(int i=0; i<_n_ptcl; i++) {
            if(part_list_n[i]>0) {
                reverse_list[i] = partner_index.size();
                partner_index.push_back(PLinker(i,-1));
            }
#ifdef HARD_DEBUG
            else {
                reverse_list[i] = -1;
            }
#endif
        }
        PS::S32 n_tot = partner_index.size();
        //PS::S32 n_groups = 0;

        PS::S32 n_mem = 0;
        for(int i=0; i<n_tot; i++) {
            // PS::S32 k = partner_index[i].first;
            if(partner_index[i].second>=0) continue;

            PS::S32 npart = connectGroups(i,i,part_list, part_list_disp, part_list_n,partner_index,reverse_list);
            group_list_n.push_back(npart);
#ifdef HARD_DEBUG
            assert(npart>0);
#endif
            group_list_disp_.push_back(n_mem);
            n_mem += npart;
            group_list.push_back(partner_index[i].first);
            PS::S32 inext=partner_index[i].second;
            PS::S32 k=1;
            while (inext!=i) {
                group_list.push_back(partner_index[inext].first);
                inext=partner_index[inext].second;
                k++;
#ifdef HARD_DEBUG
                assert(k<=npart);
#endif
            }
#ifdef HARD_DEBUG
            assert(k==npart);
            if(npart>_n_ptcl) {
                std::cerr<<"Error: connect group particle number mismatch: npart ="<<npart<<" ; _n_ptcl = "<<_n_ptcl<<std::endl;
                abort();
            }
#endif

            //n_group.push_back(npart);
            //n_groups++;
        }
    }

    PS::S32 connectGroups(const PS::S32 ip,
                          const PS::S32 iend,
                          PS::S32 part_list[],
                          PS::S32 part_list_disp[],
                          PS::S32 part_list_n[],
                          PS::ReallocatableArray<PLinker> & partner_index,
                          PS::ReallocatableArray<PS::S32> & reverse_list) {
        PS::S32 n_connected = 0;
        PS::S32 n_reduce = 0;
        PS::U32 inow=ip;
        PS::S32 kp = partner_index[ip].first;
        std::vector<PS::U32> rlist;
        for(PS::S32 j=0; j<part_list_n[kp]; j++) {
            PS::S32 inext = reverse_list[part_list[part_list_disp[kp]+j]];
            if(partner_index[inext].second<0&&inext!=iend) {
                if(partner_index[inow].second>=0) n_reduce++;
                partner_index[inow].second = inext;
                rlist.push_back(inext);
                n_connected++;
                inow = inext;
            }
        }
        if(n_connected>0) {
            partner_index[inow].second = iend;
            n_connected++;
        }

        for(PS::U32 j=0; j<rlist.size(); j++) {
            inow = rlist[j];
            PS::U32 inext = partner_index[inow].second;
            n_connected += connectGroups(inow,inext,part_list,part_list_disp,part_list_n,partner_index,reverse_list);
        }
        return n_connected - n_reduce;
    }


public:

    void searchAndMerge(Tptcl *_ptcl_in_cluster, const PS::S32 _n_ptcl) {
        PS::ReallocatableArray<PS::S32> part_list;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_disp;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_n;      ///partner list
        
        searchPartner(part_list, part_list_disp, part_list_n, _ptcl_in_cluster, _n_ptcl);
        mergeCluster(group_list_, group_list_disp_, group_list_n_, _n_ptcl, part_list.getPointer(), part_list_disp.getPointer(), part_list_n.getPointer());
    }

    PS::S32 getNumberOfGroups() const {
        return group_list_n_.size();
    }

    PS::S32* getMemberList(const std::size_t igroup) {
        return &group_list_[group_list_disp_[igroup]];
    }

    PS::S32 getGroupListSize() const {
        return group_list_.size();
    }

    PS::S32 getNumberOfGroupMembers(const std::size_t igroup) const {
        return group_list_n_[igroup];
    }

};

