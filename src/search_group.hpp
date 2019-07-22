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


    //! Set member particle status=1, return orderd particle member index list
    /* @param[in]  _bin: binary tree root
       @param[in]  _adr_ref: ptcl_org first particle address as reference to calculate the particle index.
       @param[out] _ptcl_adr_sys: particle index list in global particle system (not _ptcl_in_cluster)
     */
    template<class Tptree>
    PS::S32 setGroupMemberPars(Tptree &_bin, 
                               const Tptcl* _adr_ref, 
                               PS::S32* _ptcl_adr_sys) {
        PS::S32 nloc = 0;
        for(int i=0; i<2; i++) {
            if(_bin.member[i]->status.d!=0) {
                _bin.member[i]->changeover = _bin.changeover;
#ifdef HARD_DEBUG
                assert(_bin.member[i]->mass>0.0);
#endif                
                //_bin.member[i]->mass_bk.d = _bin.member[i]->mass;
                nloc += setGroupMemberPars(*(Tptree*)_bin.member[i], _adr_ref, &_ptcl_adr_sys[nloc]);
            }
            else {
                //if(is_top) bin.member[i]->status = -std::abs(id_offset+bin.member[i]->id*n_split);
                //else bin.member[i]->status = -std::abs(id_offset+bid*n_split);
                _bin.member[i]->status.d = 1; // used for number count later
                // To be consistent, all these paramters are set later in findGroupsAndCreateArtificalParticles
                //_bin.member[i]->changeover = _bin.changeover;
//#ifdef HARD_DEBUG
                //assert(_bin.member[i]->mass>0.0);
//#endif                
                //_bin.member[i]->mass_bk.d = _bin.member[i]->mass;
                //_bin.member[i]->r_search = _bin.r_search;
//#ifdef SPLIT_MASS
                //_bin.member[i]->mass    = 0.0;
//#endif
                _ptcl_adr_sys[nloc] = _bin.member[i]-_adr_ref;
                nloc += 1;
            }
        }
        return nloc;
    }

    //! generate kepler sampling artifical particles
    /*  @param[in]     _i_cluster: cluster index
        @param[in]     _i_group: group index in the cluster
        @param[in,out] _ptcl_in_cluster: particle data in local cluster
        @param[out]    _ptcl_new: artifical particles that will be added
        @param[in,out] _empty_list: the list of _ptcl_in_cluster that can be used to store new artifical particles, reduced when used
        @param[out]    _group_ptcl_adr_list: group member particle index list in _ptcl_in_cluster 
        @param[in]     _bin: binary tree root
        @param[in]     _id_offset: for artifical particles, the offset of starting id.
        @param[in]     _n_split: split number for artifical particles

        also store the binary parameters in mass_bk of first 4 pairs of artifical particles
        acc, ecc
        peri, tstep (integrator step estimation ),
        inc, OMG,
        omg, ecca
        tperi, stable_factor 
        mass1, mass2

        first two artifical particles status is n_members of two components.
     */
    template<class Tptree>
    void keplerOrbitGenerator(const PS::S32 _i_cluster,
                              const PS::S32 _i_group,
                              Tptcl* _ptcl_in_cluster,
                              PS::ReallocatableArray<Tptcl> & _ptcl_new,
                              PS::ReallocatableArray<PS::S32> & _empty_list,
                              PS::S32 *_group_ptcl_adr_list,
                              Tptree &_bin,
                              const PS::F64 _r_bin,
                              const PS::S64 _id_offset,
                              const PS::S32 _n_split) {
        const PS::F64 dE = 8.0*atan(1.0)/(_n_split-4);
//      const PS::F64 dE = 8.0*atan(1.0)/_n_split;
        //const PS::F64 dt = _bin.peri/_n_split;
        if (_n_split<8) {
            std::cerr<<"N_SPLIT "<<_n_split<<" to small to save binary parameters, should be >= 8!";
            abort();
        }
        PS::S32 i_cg[2]={_i_cluster, _i_group};
        PS::F64 bindata[5][2];
        //Tptcl* plist[_n_split][2];
        /*
          acc, ecc
          peri, tstep,
          inc, OMG,
          omg, ecca
          tperi, stable_factor
         */
        bindata[0][0] = _bin.semi;
        bindata[0][1] = _bin.ecc;
        bindata[1][0] = _bin.peri;
        bindata[1][1] = _bin.tstep;
        bindata[2][0] = _bin.inc; 
        bindata[2][1] = _bin.OMG; 
        bindata[3][0] = _bin.omg; 
        bindata[3][1] = _bin.ecca;
        bindata[4][0] = _bin.tperi;
        bindata[4][1] = _bin.stable_factor;

        PS::F64* pm[_n_split][2];
        PS::F64 mfactor;
        PS::F64 mnormal=0.0;

        const PS::S32 n_members = _bin.status.d;
        Tptcl* adr_ref= _ptcl_in_cluster;
        PS::S32 nbin = setGroupMemberPars(_bin, adr_ref, _group_ptcl_adr_list);
        assert(nbin==n_members);

        // Make sure the _ptcl_new will not make new array due to none enough capacity during the following loop, otherwise the p[j] pointer will point to wrong position
        _ptcl_new.reserveEmptyAreaAtLeast(2*_n_split+1-_empty_list.size());
        // First 4 is used for tidal tensor points
        // remaining is used for sample points
        for (int i=0; i<_n_split; i++) {
            Tptcl* p[2];
            for(int j=0; j<2; j++) {
                if(_empty_list.size()>0) {
                    PS::S32 k = _empty_list.back();
                    _empty_list.decreaseSize(1);
                    p[j] = &_ptcl_in_cluster[k];
                }
                else {
                    _ptcl_new.push_back(Tptcl());
                    p[j] = &_ptcl_new.back();
                }
                p[j]->mass = _bin.member[j]->mass;
#ifdef HARD_DEBUG
                assert(_bin.member[j]->mass>0);
#endif
                p[j]->id = _id_offset + (_bin.member[j]->id)*_n_split + i;
                if(i==0) p[j]->status.d = _bin.member[j]->status.d; // store the component member number 
                else if(i==1) p[j]->status.d = i_cg[j]+1; // store the i_cluster and i_group for identify artifical particles, +1 to avoid 0 value (status.d>0)
                else p[j]->status.d = (_bin.id<<ID_PHASE_SHIFT)|i; // not used, but make status.d>0
                
                if(i>=4) 
                    pm[i][j] = &(p[j]->mass);
            }
            if (i>=4) {
                PS::S32 iph = i-4;

                for (int j=0; j<2; j++) {
                    // use member changeover, if new changeover is different, record the scale ratio 
                    p[j]->changeover = _bin.member[j]->changeover;
                    if (abs(p[j]->changeover.getRin()-_bin.changeover.getRin())>1e-10) {
                        p[j]->changeover.r_scale_next = _bin.changeover.getRin()/p[j]->changeover.getRin();
                        p[j]->r_search = std::max(p[j]->r_search, _bin.r_search);
#ifdef HARD_DEBUG
                        assert(p[j]->r_search > p[j]->changeover.getRout());
#endif 
                    }
                    else p[j]->r_search = _bin.r_search;
                }

                // center_of_mass_shift(*(Tptcl*)&_bin,p,2);
                // generate particles at different orbitial phase
                OrbParam2PosVel(p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass,
                                _bin.semi, _bin.ecc, _bin.inc, _bin.OMG, _bin.omg, dE*iph);
                //DriveKeplerOrbParam(p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass, (i+1)*dt, _bin.semi, _bin.ecc, _bin.inc, _bin.OMG, _bin.omg, _bin.peri, _bin.ecca);
            }
            else {
                // use c.m. r_search 
                for (int j=0; j<2; j++) {
                    p[j]->r_search = _bin.member[j]->r_search;
                    p[j]->changeover = _bin.changeover;
                }

                ///* Assume apo-center distance is the maximum length inside box
                //   Then the lscale=apo/(2*sqrt(2))
                // */
                // PS::F64 lscale = _bin.semi*(1+_bin.ecc)*0.35;

                // Use fixed 0.5*r_bin to determine lscale
                PS::F64 lscale = 0.16*_r_bin;

                // 8 points box 
                switch(i) {
                case 0:
                    p[0]->pos = PS::F64vec(lscale, 0,      -lscale) + _bin.pos;
                    p[1]->pos = PS::F64vec(0,      lscale, -lscale) + _bin.pos;
                    break;
                case 1:
                    p[0]->pos = PS::F64vec(-lscale, 0,      -lscale) + _bin.pos;
                    p[1]->pos = PS::F64vec(0,      -lscale, -lscale) + _bin.pos;
                    break;
                case 2:
                    p[0]->pos = PS::F64vec(lscale, 0,      lscale) + _bin.pos;
                    p[1]->pos = PS::F64vec(0,      lscale, lscale) + _bin.pos;
                    break;
                case 3:
                    p[0]->pos = PS::F64vec(-lscale, 0,      lscale) + _bin.pos;
                    p[1]->pos = PS::F64vec(0,      -lscale, lscale) + _bin.pos;
                    break;
                default:
                    std::cerr<<"Error: index >= 4!\n";
                    abort();
                }
                p[0]->vel = _bin.vel;
                p[1]->vel = _bin.vel;
                p[0]->mass = p[1]->mass = 0.0;
            }
            // binary parameters
            if(i<5) {
                p[0]->mass_bk.d = bindata[i][0];
                p[1]->mass_bk.d = bindata[i][1];
            }
            else if(i==5) {
                p[0]->mass_bk.d = p[0]->mass;
                p[1]->mass_bk.d = p[1]->mass;
            }
            else if(i==6) {
                p[0]->mass_bk.d = 0.0; // indicate the order
                p[1]->mass_bk.d = 1.0;
            }
            if(i>=4) {
                PS::F64vec dvvec= p[0]->vel - p[1]->vel;
                PS::F64 odv = 1.0/std::sqrt(dvvec*dvvec);
                for(int j=0; j<2; j++) p[j]->mass *= odv;
                //if(i==0) p[j]->mass /= 2.0*_n_split;
                //else p[j]->mass /= _n_split;
                mnormal += odv;
                center_of_mass_correction(*(Tptcl*)&_bin,p,2);
#ifdef HARD_DEBUG
                //check rsearch consistence:
                PS::F64 rsearch_bin = _bin.r_search+_bin.semi*(1+_bin.ecc);
                for(int j=0; j<2; j++) {
                    PS::F64vec dp = p[j]->pos-_bin.pos;
                    PS::F64 dr = dp*dp;
                    assert(dr<=rsearch_bin*rsearch_bin);
//                    if(p[j]->id==10477) {
//                        std::cerr<<"i="<<i<<" dr="<<sqrt(dr)<<std::endl;
//                        p[j]->print(std::cerr);
//                        std::cerr<<std::endl;
//                    }
                }
#endif
            }
        }

        mfactor = 1.0/mnormal;

        for (int i=4; i<_n_split; i++) 
            for (int j=0; j<2; j++) 
                *pm[i][j] *= mfactor;
        
        //mfactor *= 0.5;
        //*pm[0][0] *= mfactor;
        //*pm[0][1] *= mfactor;
        //mfactor /= _n_split;
        // collect the member address .

        Tptcl* pcm;
        //PS::S64 pcm_adr;
        if(_empty_list.size()>0) {
            pcm = &_ptcl_in_cluster[_empty_list.back()];
            //pcm_adr = - _empty_list.back(); // if is in empty list, put negative address
            _empty_list.decreaseSize(1);
        }
        else {
            _ptcl_new.push_back(Tptcl());
            pcm = &_ptcl_new.back();
            //pcm_adr = _ptcl_new.size()-1; // if is in new, put ptcl_new address
        }
        pcm->mass_bk.d = _bin.mass;
        pcm->mass = 0.0;
        pcm->pos = _bin.pos;
        pcm->vel = _bin.vel;
        pcm->id  = - std::abs(_bin.id);
        pcm->r_search = _bin.r_search;
        pcm->changeover = _bin.changeover;

        pcm->r_search += _bin.semi*(1+_bin.ecc);  // depend on the mass ratio, the upper limit distance to c.m. from all members and artifical particles is apo-center distance

        pcm->status.d = nbin;
    }

    //! generate artifical particles,
    /*  @param[in]     _i_cluster: cluster index
        @param[in,out] _ptcl_in_cluster: particle data
        @param[in]     _n_ptcl: total number of particle in _ptcl_in_cluster.
        @param[out]    _ptcl_artifical: artifical particles that will be added
        @param[out]    _n_groups: number of groups in current cluster
        @param[in,out] _group_list: 1-D group member index array, will be reordered by the minimum distance chain for each group
        @param[in]     _group_list_disp: offset of group boundary index in _group_list
        @param[in]     _group_list_n: number of members in each group
        @param[in,out] _empty_list: the list of _ptcl_in_cluster that can be used to store new artifical particles, reduced when used
        @param[in]     _rbin: binary detection criterion radius
        @param[in]     _rin: inner radius of soft-hard changeover function
        @param[in]     _rout: outer radius of soft-hard changeover function
        @param[in]     _dt_tree: tree time step for calculating r_search
        @param[in]     _id_offset: for artifical particles, the offset of starting id.
        @param[in]     _n_split: split number for artifical particles
     */
    template<class Tptree>
    void generateNewPtcl(const PS::S32 _i_cluster,
                         Tptcl* _ptcl_in_cluster,
                         const PS::S32 _n_ptcl,
                         //PS::ReallocatableArray<PS::S32> & p_list,
                         PS::ReallocatableArray<Tptcl> & _ptcl_artifical,
                         PS::S32 &_n_groups,
                         PS::ReallocatableArray<PS::S32> & _group_list,
                         PS::ReallocatableArray<PS::S32> & _group_list_disp,
                         PS::ReallocatableArray<PS::S32> & _group_list_n,
                         PS::ReallocatableArray<PS::S32> & _empty_list,
                         const PS::F64 _rbin,
                         const PS::F64 _rin,
                         const PS::F64 _rout,
                         const PS::F64 _dt_tree,
                         const PS::S64 _id_offset,
                         const PS::S32 _n_split){
#ifdef HARD_DEBUG
//        assert(ptcl.size()==0);
//        assert(group_ptcl.size()==0);
//        assert(_n_ptcl<=_n_ptcl_in_cluster);
//        assert(ptcl_map.size()==0);
#endif
        // reset all single status to 0
        //for (int i=0; i<p_list.size(); i++) _ptcl_in_cluster[p_list[i]].status.d = 0;

        //_n_groups = _group_list_n.size();
        PS::S32 group_ptcl_adr_list[_n_ptcl];
        PS::S32 group_ptcl_adr_offset=0;
        _n_groups = 0;
        for (int i=0; i<_group_list_n.size(); i++) {
            PS::ReallocatableArray<Tptree> bins;   // hierarch binary tree
            PS::ReallocatableArray<Tptree*> stab_bins; // stable checked binary tree
            bins.reserve(_n_groups);
            stab_bins.reserve(_n_groups);
#ifdef HARD_DEBUG
            assert(_group_list_n[i]<ARRAY_ALLOW_LIMIT);
#endif        
            bins.resizeNoInitialize(_group_list_n[i]-1);
            // build hierarch binary tree from the minimum distant neighbors

            const PS::S32 group_start = _group_list_disp[i];
            const PS::S32 group_n = _group_list_n[i];

            keplerTreeGenerator(bins.getPointer(), &_group_list[group_start], group_n, _ptcl_in_cluster, _rin, _rout, _dt_tree);
         
            // reset status to 0
            for (int j=0; j<group_n; j++) _ptcl_in_cluster[_group_list[group_start+j]].status.d=0;

            // stability check and break groups
            bool stab_flag=stabilityCheck<Tptcl>(stab_bins, bins.back(), _rbin, _rin, _rout, _dt_tree);
            if(stab_flag) stab_bins.push_back(&bins.back());
            
            for (int i=0; i<stab_bins.size(); i++) {
                keplerOrbitGenerator(_i_cluster, _n_groups, _ptcl_in_cluster, _ptcl_artifical, _empty_list, &group_ptcl_adr_list[group_ptcl_adr_offset], *stab_bins[i], _rbin, _id_offset, _n_split);
                group_ptcl_adr_offset += stab_bins[i]->status.d;
                _n_groups++;
            }
        }

        assert(group_ptcl_adr_offset<=_n_ptcl);

        // Reorder the ptcl that group member come first
        PS::S32 ptcl_list_reorder[_n_ptcl];
        for (int i=0; i<_n_ptcl; i++) ptcl_list_reorder[i] = i;
 
        // shift single after group members
        PS::S32 i_single_front=group_ptcl_adr_offset;
        PS::S32 i_group = 0;
        while (i_group<group_ptcl_adr_offset) {
            // if single find inside group_ptcl_adr_offset, exchange single with group member out of the offset
            if(_ptcl_in_cluster[i_group].status.d==0) {
                while(_ptcl_in_cluster[i_single_front].status.d==0) {
                    i_single_front++;
                    assert(i_single_front<_n_ptcl);
                }
                // Swap index
                PS::S32 plist_tmp = ptcl_list_reorder[i_group];
                ptcl_list_reorder[i_group] = ptcl_list_reorder[i_single_front];
                ptcl_list_reorder[i_single_front] = plist_tmp;
                i_single_front++; // avoild same particle be replaced
            }
            i_group++;
        }

#ifdef HARD_DEBUG
        // check whether the list is correct
        PS::S32 plist_new[group_ptcl_adr_offset];
        for (int i=0; i<group_ptcl_adr_offset; i++) plist_new[i] = group_ptcl_adr_list[i];
        std::sort(plist_new, plist_new+group_ptcl_adr_offset, [](const PS::S32 &a, const PS::S32 &b) {return a < b;});
        std::sort(ptcl_list_reorder, ptcl_list_reorder+group_ptcl_adr_offset, [](const PS::S32 &a, const PS::S32 &b) {return a < b;});
        for (int i=0; i<group_ptcl_adr_offset; i++) assert(ptcl_list_reorder[i]==plist_new[i]);
#endif        

        // overwrite the new ptcl list for group members by reorderd list
        for (int i=0; i<group_ptcl_adr_offset; i++) ptcl_list_reorder[i] = group_ptcl_adr_list[i];

        // templately copy ptcl data
        Tptcl ptcl_tmp[_n_ptcl];
        for (int i=0; i<_n_ptcl; i++) ptcl_tmp[i]=_ptcl_in_cluster[i];

        // reorder ptcl
        for (int i=0; i<_n_ptcl; i++) _ptcl_in_cluster[i]=ptcl_tmp[ptcl_list_reorder[i]];

        for (int i=0; i<_empty_list.size(); i++) {
            PS::S32 ik = _empty_list[i];
            _ptcl_in_cluster[ik].mass = 0.0;
            _ptcl_in_cluster[ik].id = -1;
            _ptcl_in_cluster[ik].status.d = -1;
        }

    }

//    void checkRoutChange(PS::ReallocatableArray<RCList> & r_out_change_list,
//                         PS::ReallocatableArray<PS::ReallocatableArray<PS::S32>> & group_list,
//                         Tptcl* ptcl){
//        for (int i=0; i<group_list.size(); i++) {
//            PS::S64 r_out_max = 0.0;
//            for (int j=0; j<group_list[i].size(); i++) {
//                PS::S32 k = group_list[i][j];
//                r_out_max = std::max(r_out_max,ptcl[k].r_out);
//            }
//            for (int j=0; j<group_list[i].size(); i++) {
//                PS::S32 k = group_list[i][j];
//                if(ptcl[k].r_out != r_out_max) {
//                    r_out_change_list.push_back(RCList(k,ptcl[k].r_out));
//#ifdef HARD_DEBUG
//                    std::cerr<<"Rout change detected, p["<<k<<"].r_out: "<<ptcl[k].r_out<<" -> "<<r_out_max<<std::endl;
//#endif
//                    ptcl[k].r_out = r_out_max;
//                }
//            }
//        }
//    }


public:

    void searchAndMerge(Tptcl *_ptcl_in_cluster, const PS::S32 _n_ptcl) {
        PS::ReallocatableArray<PS::S32> part_list;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_disp;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_n;      ///partner list
        
        searchPartner(part_list, part_list_disp, part_list_n, _ptcl_in_cluster, _n_ptcl);
        mergeCluster(group_list_, group_list_disp_, group_list_n_, _n_ptcl, part_list.getPointer(), part_list_disp.getPointer(), part_list_n.getPointer());
    }

    //! generate artifical particles,
    /*  @param[in]     _i_cluster: cluster index
        @param[in,out] _ptcl_in_cluster: particle data
        @param[in]     _n_ptcl: total number of particle in _ptcl_in_cluster.
        @param[out]    _ptcl_artifical: artifical particles that will be added
        @param[out]    _n_groups: number of groups in current cluster
        @param[in]     _rbin: binary detection criterion radius
        @param[in]     _rin: inner radius of soft-hard changeover function
        @param[in]     _rout: outer radius of soft-hard changeover function
        @param[in]     _dt_tree: tree time step for calculating r_search
        @param[in]     _id_offset: for artifical particles, the offset of starting id.
        @param[in]     _n_split: split number for artifical particles
    */
    void generateList(const PS::S32 _i_cluster,
                      Tptcl *_ptcl_in_cluster, 
                      const PS::S32 _n_ptcl,
                      PS::ReallocatableArray<Tptcl> & _ptcl_artifical,
                      PS::S32 &_n_groups,
                      const PS::F64 _rbin,
                      const PS::F64 _rin,
                      const PS::F64 _rout,
                      const PS::F64 _dt_tree,
                      const PS::S64 _id_offset,
                      const PS::S32 _n_split) {
        if (_n_split>(1<<ID_PHASE_SHIFT)) {
            std::cerr<<"Error! ID_PHASE_SHIFT is too small for phase split! shift bit: "<<ID_PHASE_SHIFT<<" n_split: "<<_n_split<<std::endl;
            abort();
        }
        PS::ReallocatableArray<PS::S32> emtpy_list;
        generateNewPtcl<PtclTree<Tptcl>>(_i_cluster, _ptcl_in_cluster, _n_ptcl, _ptcl_artifical, _n_groups, group_list_, group_list_disp_, group_list_n_, emtpy_list, _rbin, _rin, _rout, _dt_tree, _id_offset, _n_split);
        //searchPerturber(pert_list_, _ptcl_in_cluster, _n_ptcl);
    }

    PS::S32 getNumOfGroups() const {
        return group_list_n_.size();
    }

    const PS::S32* getGroup(const std::size_t igroup) const {
        return &group_list_[group_list_disp_[igroup]];
    }

    PS::S32 getGroupListSize() const {
        return group_list_.size();
    }

    PS::S32 getGroupN(const std::size_t igroup) const {
        return group_list_n_[igroup];
    }

};


//! Group paramter class
/* get artifical particle group parameters
 */
class GroupPars{
public:
    PS::S32 id;         ///> group id corresponding to the first member
    PS::S32 i_cluster;  ///> cluster index
    PS::S32 i_group;    ///> group index in cluster
    PS::S32 n_members;  ///> number of members
    PS::S32 n_members_1st; ///> number of members in first component
    PS::S32 n_members_2nd; ///> number of members in second component
    PS::S32 offset_cm;   ///> c.m. index offset in group
    PS::S32 offset_orb;  ///> orbital partical offset in group
    PS::S32 offset_tt;   ///> tital tensor partical offset in group
    PS::S32 n_ptcl_artifical; ///> artifical particle number

    GroupPars(const PS::S32 _n_split): id(-10), i_cluster(-1), i_group(-1), n_members(0), n_members_1st(0), n_members_2nd(0), offset_cm(2*_n_split), 
                                       offset_orb(8), 
                                       offset_tt(0), n_ptcl_artifical(2*_n_split+1) {}
    GroupPars() {init(0);}

    void init(const PS::S32 _n_split) {
        id=-10;
        i_cluster=-1;
        i_group=-1;
        n_members=0;
        n_members_1st=0;
        n_members_2nd=0;
        offset_cm=2*_n_split;
        offset_orb=8;
        offset_tt=0;
        n_ptcl_artifical=2*_n_split+1;
    }

    // assume the binary information stored in artifical star mass_bk
    template <class Tptcl>
    void getBinPars(Binary &bin, const Tptcl* _ptcl_artifical) {
        bin.semi   = _ptcl_artifical[0].mass_bk.d;
        bin.ecc  = _ptcl_artifical[1].mass_bk.d;
        bin.peri = _ptcl_artifical[2].mass_bk.d;
        bin.tstep= _ptcl_artifical[3].mass_bk.d;
        bin.inc  = _ptcl_artifical[4].mass_bk.d;
        bin.OMG  = _ptcl_artifical[5].mass_bk.d;
        bin.omg  = _ptcl_artifical[6].mass_bk.d;
        bin.ecca = _ptcl_artifical[7].mass_bk.d;
        bin.tperi= _ptcl_artifical[8].mass_bk.d;
        bin.stable_factor= _ptcl_artifical[9].mass_bk.d;
        bin.m1   = _ptcl_artifical[10].mass_bk.d;
        bin.m2   = _ptcl_artifical[11].mass_bk.d;
    }

    //! return group parameters
    /* @param[in] _ptcl_artifical: ptcl artifical particle group
    */
    template <class Tptcl>
    void getGroupIndex(Tptcl* _ptcl_artifical) {
        n_members_1st = _ptcl_artifical[0].status.d;
        n_members_2nd = _ptcl_artifical[1].status.d;
        i_cluster = _ptcl_artifical[2].status.d-1;
        i_group   = _ptcl_artifical[3].status.d-1;
        n_members = _ptcl_artifical[offset_cm].status.d;
        id        =-_ptcl_artifical[offset_cm].id;
    }
    
};
