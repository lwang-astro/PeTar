#pragma once
#include<particle_simulator.hpp>
#include<unordered_map>
#include<map>
#include"ptcl.hpp"
#include"Common/binary_tree.h"

//extern const PS::F64 SAFTY_OFFSET_FOR_SEARCH;

#ifndef ARRAY_ALLOW_LIMIT
#define ARRAY_ALLOW_LIMIT 1000000000
#endif

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
template <class Ttree>
void SetRankComm(const PS::F64ort pos_domain[], Ttree &tree,
                 PS::ReallocatableArray<PS::S32> & rank_neighbor,
                 const PS::F64 eps_sq = 0.0){
    rank_neighbor.clearSize();
//    static const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = 1.01;
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
//    const PS::F64ort pos_my_domain = pos_domain[my_rank];
    PS::F64ort pos_my_domain = tree.getOuterBoundaryOfLocalTree();
    PS::F64ort pos_my_domain_in = tree.getInnerBoundaryOfLocalTree();
    PS::F64ort * pos_domain_out = new PS::F64ort[n_proc];
    MPI_Allgather(&pos_my_domain, 1, PS::GetDataType<PS::F64ort>(),
                  pos_domain_out, 1, PS::GetDataType<PS::F64ort>(), MPI_COMM_WORLD);
    PS::F64ort * pos_domain_in = new PS::F64ort[n_proc];
    MPI_Allgather(&pos_my_domain_in, 1, PS::GetDataType<PS::F64ort>(),
                  pos_domain_in,  1, PS::GetDataType<PS::F64ort>(), MPI_COMM_WORLD);
//    const PS::F64 r_crit_sq = (EPJSoft::r_search * EPJSoft::r_search + eps_sq)*SAFTY_FACTOR_FOR_SEARCH_SQ*1.01;
//    pos_domain[my_rank].getOuterBoundaryOfLocalTree
    rank_neighbor.clearSize();
    for(int i=0; i<n_proc; i++){
        if(i==my_rank) continue;
        //else if( r_crit_sq >= pos_my_domain.getDistanceMinSQ(pos_domain[i]) ){
        else if( pos_my_domain.contained(pos_domain_in[i]) || pos_my_domain_in.contained(pos_domain_out[i])){
            rank_neighbor.push_back(i);
        }
//        rank_neighbor.push_back(i);
    }
    delete [] pos_domain_in;
    delete [] pos_domain_out;
}
#endif

struct Cluster{
    PS::S32 id_;
    PS::S32 n_ptcl_; // include outer particle
    PS::S32 n_ptcl_stored_; // exclude outer particle
    PS::S32 adr_head_; // for cluster_isolated, adr_sys_multi_cluster_isolated_
    PS::S32 rank_;
    Cluster(): id_(-1), n_ptcl_(0), n_ptcl_stored_(0), adr_head_(-1), rank_(-1){}
    Cluster(const PS::S32 _id, const PS::S32 _n_ptcl, const PS::S32 _n_ptcl_stored, 
            const PS::S32 _adr_head, const PS::S32 _rank):
        id_(_id), n_ptcl_(_n_ptcl), n_ptcl_stored_(_n_ptcl_stored), 
        adr_head_(_adr_head), rank_(_rank){}
    void clear(){
        id_ = -1; n_ptcl_ = 0; n_ptcl_stored_ = 0; adr_head_ = -1; rank_ = -1;
    }
    void dump(){
        std::cout<<" id="<<id_<<" n_ptcl="<<n_ptcl_
                 <<" n_ptcl_stored="<<n_ptcl_stored_
                 <<" adr_head="<<adr_head_
                 <<" rank="<<rank_<<std::endl;
    }
};

struct PtclCluster{
    PS::S32 id_;
    PS::S32 adr_sys_;
    PS::S32 adr_ngb_head_; // adr of id_ngb_multi_cluster
    PS::S32 n_ngb_;
    bool flag_searched_;
    PtclCluster * next_;
    PS::S32 rank_org_;
    PS::S32 id_cluster_;
    PtclCluster(): id_(-10), adr_sys_(-1), adr_ngb_head_(-1), n_ngb_(0), 
                   flag_searched_(false), next_(NULL), rank_org_(-1), id_cluster_(id_){}
    PtclCluster(const PS::S32 _id,    const PS::S32 _adr_sys,    const PS::S32 _adr_ngb_head, 
                const PS::S32 _n_ngb, const bool _flag_searched, PtclCluster * const _next, 
                PS::S32 _rank_org):
        id_(_id), adr_sys_(_adr_sys), adr_ngb_head_(_adr_ngb_head),  
        n_ngb_(_n_ngb), flag_searched_(_flag_searched), next_(_next), rank_org_(_rank_org), 
        id_cluster_(_id){}
    void clear(){
        id_ = -10; adr_sys_ = -1; adr_ngb_head_ = -1; n_ngb_ = 0; 
        flag_searched_ = false; next_ = NULL; rank_org_ = -1; id_cluster_ = id_;
    }
    void dump(){
        std::cout<<" id="<<id_<<" adr_sys="<<adr_sys_<<" adr_ngb_head_="<<adr_ngb_head_
                 <<" n_ngb="<<n_ngb_<<" flag_searched_="<<flag_searched_
                 <<" rank_org="<<rank_org_<<" id_cluster_="<<id_cluster_<<std::endl;
    }
    void setIdClusterImpl(PS::ReallocatableArray<PtclCluster> & ptcl,
                          PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & adr_ngb){
        for(PS::S32 i=0; i<n_ngb_; i++){
            id_cluster_ = (ptcl[adr_ngb[adr_ngb_head_+i].second].id_cluster_ < id_cluster_) 
                ? ptcl[adr_ngb[adr_ngb_head_+i].second].id_cluster_ : id_cluster_;
        }
    }
};

class PtclOuter{
public:
    PS::S32 id_;
    PS::S32 id_ngb_;
    PS::S32 rank_org_;
    PtclOuter(): id_(-10), id_ngb_(-1), rank_org_(-1){}
    PtclOuter(const PS::S32 _id, const PS::S32 _id_ngb, const PS::S32 _rank_org): 
        id_(_id), id_ngb_(_id_ngb), rank_org_(_rank_org){}
    void clear(){
        id_ = -10; id_ngb_ = -1; rank_org_ = -1;
    }
    void dump(){
        std::cout<<" id="<<id_<<" id_ngb_="<<id_ngb_<<" rank_org_="<<rank_org_<<std::endl;
    }
};

class Mediator{
public:
    PS::S32 id_;
    PS::S32 adr_sys_;
    PS::S32 adr_pcluster_;
    PS::S32 id_cluster_;
    PS::S32 rank_send_;
    Mediator():id_(-10), adr_sys_(-1), adr_pcluster_(-1), id_cluster_(-1), rank_send_(-1){} 
    Mediator(const PS::S32 _id, const PS::S32 _adr_sys, 
             const PS::S32 _adr_pcluster, const PS::S32 _id_cluster,
             const PS::S32 _rank_send):
        id_(_id), adr_sys_(_adr_sys), adr_pcluster_(_adr_pcluster), id_cluster_(_id_cluster),
        rank_send_(_rank_send){}
    void clear(){
        id_ = -10; adr_sys_ = -1; adr_pcluster_ = -1; id_cluster_ = -1; rank_send_ = -1;
    }
    void dump(){
        std::cout<<" id="<<id_<<" adr_sys="<<adr_sys_
                 <<" adr_pcluster_="<<adr_pcluster_<<" id_cluster="<<id_cluster_
                 <<" rank_send_="<<rank_send_<<std::endl;
    }
};

class PtclComm: public Ptcl{
public:
    PS::S32 id_cluster;

    template <class Tp>
    PtclComm(const Tp &p): Ptcl(p), id_cluster(-1) {}
    PtclComm() {}

    void print(std::ostream & fout){
        Ptcl::print(fout);
        fout<<" id_cluster="<<id_cluster<<std::endl;
    }
};

/*
  class exchangeInfo{
  PS::ReallocatableArray<PS::S32> n_send_;
  PS::ReallocatableArray<PS::S32> n_recv_;
  PS::ReallocatableArray<PS::S32> n_disp_send_;
  PS::ReallocatableArray<PS::S32> n_disp_recv_;
  void set(const PS::S32 n_proc_send, const PS::S32 n_proc_recv){
  n_send_.resizeNoInitialize(n_proc_send);
  n_disp_send_.resizeNoInitialize(n_proc_send+1);
  n_recv_.resizeNoInitialize(n_proc_recv);
  n_disp_recv_.resizeNoInitialize(n_proc_recv+1);
  }
  };
*/


class SearchCluster{
public:
    PS::ReallocatableArray<PS::S32> n_ptcl_in_multi_cluster_isolated_;
    PS::ReallocatableArray<PS::S32> n_ptcl_in_multi_cluster_isolated_offset_;
    PS::ReallocatableArray<Mediator> mediator_sorted_id_cluster_;
    PS::ReallocatableArray<PtclComm> ptcl_recv_;
    PS::ReallocatableArray<PS::S32> adr_sys_multi_cluster_isolated_;
private:
    PS::ReallocatableArray<Cluster> cluster_comm_;
    std::unordered_map<PS::S32, PS::S32> id_to_adr_pcluster_;
    // 1st: adr of self 2nd: adr of ngb
    PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > adr_ngb_multi_cluster_;
    PS::ReallocatableArray<PS::S32> * adr_sys_one_cluster_;
    PS::ReallocatableArray<PtclCluster> * ptcl_cluster_;
    PS::ReallocatableArray<PS::S32> id_cluster_send_;
    PS::ReallocatableArray<PS::S32> id_cluster_recv_;
    PS::ReallocatableArray<PS::S32> adr_pcluster_send_;
    PS::ReallocatableArray<PS::S32> adr_pcluster_recv_;
    PS::ReallocatableArray<PS::S32> n_cluster_send_;
    PS::ReallocatableArray<PS::S32> n_cluster_recv_;
    PS::ReallocatableArray<PS::S32> n_cluster_disp_send_;
    PS::ReallocatableArray<PS::S32> n_cluster_disp_recv_;
    PS::ReallocatableArray<PS::S32> rank_send_cluster_;
    PS::ReallocatableArray<PS::S32> rank_recv_cluster_;
    PS::S32 n_pcluster_self_node_;
    PS::ReallocatableArray<PS::S32> adr_sys_ptcl_send_;
    PS::ReallocatableArray<PtclComm> ptcl_send_;
    PS::ReallocatableArray<PS::S32> rank_send_ptcl_;
    PS::ReallocatableArray<PS::S32> n_ptcl_send_;
    PS::ReallocatableArray<PS::S32> n_ptcl_disp_send_;
    PS::ReallocatableArray<PS::S32> rank_recv_ptcl_;
    PS::ReallocatableArray<PS::S32> n_ptcl_recv_;
    PS::ReallocatableArray<PS::S32> n_ptcl_disp_recv_;
    template<class T>
    void packDataToThread0(T * data){
        const PS::S32 n_thread = PS::Comm::getNumberOfThread();
        PS::S32 size_data_th0 = data[0].size();
        PS::S32 size_data = size_data_th0;
        for(PS::S32 i=1; i<n_thread; i++){
            size_data += data[i].size();
        }
#ifdef HARD_DEBUG
        assert(size_data<ARRAY_ALLOW_LIMIT);
#endif        
        
        data[0].resizeNoInitialize(size_data);
#pragma omp parallel
        {
            const PS::S32 ith = PS::Comm::getThreadNum();
            if(ith > 0){
                PS::S32 offset = size_data_th0;
                for(PS::S32 i=1; i<ith; i++) offset += data[i].size();
                for(PS::S32 i=0; i<data[ith].size(); i++) data[0][i+offset] = data[ith][i];
            }
        }
    }

    void searchClusterImpl(PtclCluster * target,
                           PtclCluster * p_first,
                           PtclCluster *& p_top,
                           PS::S32 & n_ptcl_in_cluster,
                           PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > & adr_ngb_array,
                           bool & fg_isolated){
        if(target->adr_sys_ < 0) fg_isolated = false;
        target->flag_searched_ = true;
        target->next_ = NULL;
        n_ptcl_in_cluster++;
        PS::S32 target_adr = target->adr_ngb_head_;
        const PS::S32 n_ngb = target->n_ngb_;
        for(PS::S32 i=0; i<n_ngb; i++){
            PS::S32 adr_ngb = adr_ngb_array[target_adr+i].second;
            if( (p_first+adr_ngb)->flag_searched_ == false){
                p_top->next_ = p_first+adr_ngb;
                p_top = p_first+adr_ngb;
                searchClusterImpl(p_first+adr_ngb, p_first, p_top, n_ptcl_in_cluster,
                                  adr_ngb_array, fg_isolated);
            }
        }
    }

    void setNgbAdrHead(PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > id_ngb_multi_cluster[]){
        const PS::S32 size_ngb = id_ngb_multi_cluster[0].size();
        PS::S32 n_cnt = 0;
        for(PS::S32 i=0; i<size_ngb; i++){
            if(id_ngb_multi_cluster[0][i].first == ptcl_cluster_[0][n_cnt].id_){
                ptcl_cluster_[0][n_cnt].adr_ngb_head_ = i;
                n_cnt++;
            }
        }
    }

    void mergePtclCluster(const PS::ReallocatableArray<PtclOuter> ptcl_outer[],
                          PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > id_ngb_multi_cluster[]){
        if(ptcl_outer[0].size() == 0) return;
        const PS::S32 size = ptcl_outer[0].size();
        PS::S32 id_tmp = ptcl_outer[0][0].id_;
        PS::S32 adr_ngb_head = id_ngb_multi_cluster[0].size();
        PS::S32 rank_tmp = ptcl_outer[0][0].rank_org_;
        PS::S32 n_cnt = 0;
        for(PS::S32 i=0; i<size; i++){
            if(id_tmp != ptcl_outer[0][i].id_){
                ptcl_cluster_[0].push_back( PtclCluster(id_tmp, -1, adr_ngb_head, n_cnt, false, NULL, rank_tmp) );
                id_tmp = ptcl_outer[0][i].id_;
                rank_tmp = ptcl_outer[0][i].rank_org_;
                adr_ngb_head += n_cnt;
                n_cnt = 0;
            }
            id_ngb_multi_cluster[0].push_back( std::pair<PS::S32, PS::S32>
                                               (ptcl_outer[0][i].id_,
                                                ptcl_outer[0][i].id_ngb_));
            n_cnt++;
        }
        ptcl_cluster_[0].push_back( PtclCluster(id_tmp, -1, adr_ngb_head, n_cnt, false, NULL, rank_tmp) );
    }

    struct OPEqualID{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_ == right.id_;
        }
    };
    struct OPEqualSecond{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.second == right.second;
        }
    };
    struct OPLessID{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_ < right.id_;
        }
    };
    struct OPLessIDCluster{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.id_cluster_ < right.id_cluster_;
        }
    };
    struct OPLessRankOrg{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.rank_org_ < right.rank_org_;
        }
    };
    struct OPLessRank{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.rank_ < right.rank_;
        }
    };
    struct OPLessFirst{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.first < right.first;
        }
    };
    struct OPLessSecond{
        template<class T> bool operator() (const T & left, const T & right) const {
            return left.second < right.second;
        }
    };

public:
    void initialize(){
        const PS::S32 n_thread = PS::Comm::getNumberOfThread();
        adr_sys_one_cluster_  = new PS::ReallocatableArray<PS::S32>[n_thread];
        ptcl_cluster_ = new PS::ReallocatableArray<PtclCluster>[n_thread];
    }


    //! identify whether the neighbor satisfy velocity criterion
    /*! 
      @param[out] _index: neighbor index pass check
      @param[in] _pi: particle i
      @param[in] _pb: neighbor particles (including self)
      @param[in] _nb: number of neighbors
      @param[in] _G:  gravitational constant
      @param[in] _radius_factor: radius_factor to check neighbors (peri*radius_factor)
      \return new neighbor number
     */
    template<class Tpsoft, class Tepj>
    PS::S32 checkNeighborWithVelocity(PS::S32* _index, Tpsoft& _pi, Tepj *_pb, const PS::S32 _nb, const PS::F64 _G, const PS::F64 _radius_factor) {
        
        PS::S32 n_nb_new = 0;
        PS::F64 r_crit_i = _pi.changeover.getRout();

        ParticleBase pi;
        pi.DataCopy(_pi);
        auto& pi_cm = _pi.group_data.cm;

#ifdef CLUSTER_DEBUG
        assert(pi_cm.mass>=0.0);
#endif
        // if i particle is member
        if (pi_cm.mass>0.0) {
            pi.mass   = pi_cm.mass;
            pi.vel[0] = pi_cm.vel.x;
            pi.vel[1] = pi_cm.vel.y;
            pi.vel[2] = pi_cm.vel.z;
#ifdef CLUSTER_DEBUG
            if (_pi.group_data.artificial.status==0.0 && _pi.group_data.artificial.mass_backup!=0.0) 
                std::cout<<"Warning: may not be pcm data! idi "<<_pi.id<<" status="<<_pi.group_data.artificial.status<<" mass_bk="<<_pi.group_data.artificial.mass_backup<<std::endl;
#endif
        }

        for (PS::S32 j=0; j<_nb; j++) {
            if (_pi.id==_pb[j].id) continue;
            auto& pbj_cm = _pb[j].group_data.cm;
#ifdef CLUSTER_DEBUG
            assert(pbj_cm.mass>=0.0);
#endif
            ParticleBase pj;
            pj.pos = _pb[j].pos;
            // particle member case
            if (pbj_cm.mass>0.0) {
                pj.vel[0] = pbj_cm.vel.x;
                pj.vel[1] = pbj_cm.vel.y;
                pj.vel[2] = pbj_cm.vel.z;
                pj.mass   = pbj_cm.mass; 
#ifdef CLUSTER_DEBUG
                if (_pb[j].group_data.artificial.status==0.0 && _pb[j].group_data.artificial.mass_backup!=0.0) 
                    std::cout<<"Warning: may not be pcm data! idj "<<_pb[j].id<<" status="<<_pb[j].group_data.artificial.status<<" mass_bk="<<_pb[j].group_data.artificial.mass_backup<<std::endl;
#endif
            }
            else {
                pj.mass= _pb[j].mass;
                pj.vel = _pb[j].vel;
            }
            
            PS::F64 semi,ecc, r,rv;
            COMM::Binary::particleToSemiEcc(semi, ecc, r, rv, pi, pj, _G);
            PS::F64 peri = semi*(1-ecc);
            PS::F64 r_crit_j = _pb[j].r_out;
            PS::F64 r_crit_max = _radius_factor*std::max(r_crit_i, r_crit_j);
            
            if (peri < r_crit_max && !(rv<0 && r>r_crit_max) ) {
                _index[n_nb_new++] = j;
            }
#ifdef CLUSTER_DEBUG_PRINT
            else 
                std::cerr<<"Reject id.i "<<_pi.id<<" id.j "<<_pb[j].id<<" peri "
                         <<peri<<" r_out.i "<<r_crit_i<<" r_out.j "<<r_crit_j<<" dr "<<r<<" drdv "<<rv
                         <<" cm.vel.i "<<pi.vel[0]<<" "<<pi.vel[1]<<" "<<pi.vel[2]<<" m.i "<<pi.mass<<" "
                         <<" status.i "<<_pi.group_data.artificial.status<<" mass_bk.i "<<_pi.group_data.artificial.mass_backup
                         <<" cm.vel.j "<<pj.vel[0]<<" "<<pj.vel[1]<<" "<<pj.vel[2]<<" m.j "<<pj.mass<<" "
                         <<" status.j "<<_pb[j].group_data.artificial.status<<" mass_bk.j "<<_pb[j].group_data.artificial.mass_backup
                         <<std::endl;
#endif            
        }
        
        return n_nb_new;
    }

    //! search neighbors and separate isolated and multiple clusters
    template<class Tsys, class Ttree, class Tepj>
    void searchNeighborOMP(Tsys & sys,
                           Ttree & tree,
                           const PS::F64ort pos_domain[],
                           const PS::F64 _G,
                           const PS::F64 _radius_factor){
        static PS::ReallocatableArray<PtclOuter> * ptcl_outer = NULL;
        static PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > * id_ngb_multi_cluster = NULL;
        const PS::S32 n_thread = PS::Comm::getNumberOfThread();
        if(ptcl_outer==NULL) ptcl_outer = new PS::ReallocatableArray<PtclOuter>[n_thread];
        if(id_ngb_multi_cluster==NULL) id_ngb_multi_cluster = new PS::ReallocatableArray< std::pair<PS::S32, PS::S32> >[n_thread];
        const PS::S32 my_rank = PS::Comm::getRank();
        //        const PS::S32 n_proc_tot = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc = sys.getNumberOfParticleLocal();
#pragma omp parallel
        {
            const PS::S32 ith = PS::Comm::getThreadNum();
//            const PS::S32 ith = 0;
            adr_sys_one_cluster_[ith].clearSize();
            id_ngb_multi_cluster[ith].clearSize();
            for(PS::S32 i=0; i<ptcl_cluster_[ith].size(); i++) ptcl_cluster_[ith][i].clear();
            ptcl_cluster_[ith].clearSize();
            for(PS::S32 i=0; i<ptcl_outer[ith].size(); i++) ptcl_outer[ith][i].clear();
            ptcl_outer[ith].clearSize();
#pragma omp for
            for(PS::S32 i=0; i<n_loc; i++){
                if(sys[i].n_ngb == 1){
                    // no neighbor
                    adr_sys_one_cluster_[ith].push_back(i);
#ifdef CLUSTER_DEBUG
                    assert(sys[i].group_data.artificial.isSingle());

                    Tepj * nbl = NULL;
                    PS::S32 n_ngb_tree_i = tree.getNeighborListOneParticle(sys[i], nbl) ;
                    if(n_ngb_tree_i!=sys[i].n_ngb) {
                        std::cerr<<"Error: particle "<<i<<" Tree neighbor search number ("<<n_ngb_tree_i<<") is inconsistent with force kernel neighbor number ("<<sys[i].n_ngb<<")!\n Neighbor ID from tree: ";
                        for(int j=0; j<n_ngb_tree_i; j++) {
                            std::cerr<<(nbl+j)->id<<" ";
                        }
                        std::cerr<<std::endl;
                        abort();
                    }
#endif
                }
                else{
                    // has neighbor
                    Tepj * nbl = NULL;
#ifdef CLUSTER_DEBUG
                    PS::S32 n_ngb_force_i = sys[i].n_ngb;
#endif
                    sys[i].n_ngb = tree.getNeighborListOneParticle(sys[i], nbl) - 1;
#ifdef CLUSTER_DEBUG
                    if(n_ngb_force_i<sys[i].n_ngb+1) {
                        std::cerr<<"Error: particle "<<i<<" Tree neighbor search number ("<<sys[i].n_ngb+1<<") is inconsistent with force kernel neighbor number ("<<n_ngb_force_i<<")!\n Neighbor ID from tree: ";
                        for(int j=0; j<sys[i].n_ngb+1; j++) {
                            std::cerr<<(nbl+j)->id<<" ";
                        }
                        std::cerr<<std::endl;
                        abort();
                    }
#endif

                    // no neighbor
                    if(sys[i].n_ngb == 0){
#ifdef CLUSTER_DEBUG
                        assert(sys[i].group_data.artificial.isSingle());
#endif
                        adr_sys_one_cluster_[ith].push_back(i);
                        continue;
                    }
                    
                    PS::S32 adr_ngb_head_i = id_ngb_multi_cluster[ith].size();
                    PS::S32 n_ngb_i = sys[i].n_ngb;
                    
#ifdef CLUSTER_VELOCITY
                    // Use velocity criterion to select neighbors
                    PS::S32 neighbor_index[n_ngb_i];
                    PS::S32 n_ngb_check = checkNeighborWithVelocity(neighbor_index, sys[i], nbl, n_ngb_i+1, _G, _radius_factor);
#ifdef CLUSTER_DEBUG
                    assert(n_ngb_check>=0);
#endif
                    sys[i].n_ngb = n_ngb_check;

                    // no neighbor
                    if (n_ngb_check==0) {
                        adr_sys_one_cluster_[ith].push_back(i);
                        continue;
                    }
                    else {
                        for(PS::S32 j=0; j<n_ngb_check; j++) {
                            PS::S32 index = neighbor_index[j];
                            auto* nbj = nbl + index;

                            id_ngb_multi_cluster[ith].push_back( std::pair<PS::S32, PS::S32>(sys[i].id, nbj->id) );
                            if( nbj->rank_org != my_rank ){
                                ptcl_outer[ith].push_back(PtclOuter(nbj->id, sys[i].id, nbj->rank_org));
                            }
                        }
                        ptcl_cluster_[ith].push_back( PtclCluster(sys[i].id, i, adr_ngb_head_i, n_ngb_check, false, NULL, my_rank) );
                    }
#else
                    // use neighbor lists from tree neighbor search
                    for (int j=0; j<n_ngb_i+1; j++) {
                        auto* nbj = nbl + j;
                        if (sys[i].id==nbj->id) continue;
                        
                        id_ngb_multi_cluster[ith].push_back( std::pair<PS::S32, PS::S32>(sys[i].id, nbj->id) );
                        if( nbj->rank_org != my_rank ){
                            ptcl_outer[ith].push_back(PtclOuter(nbj->id, sys[i].id, nbj->rank_org));
                        }
                    }
                    ptcl_cluster_[ith].push_back( PtclCluster(sys[i].id, i, adr_ngb_head_i, n_ngb_i, false, NULL, my_rank) );
#endif
                }
            }
        } // end of OMP parallel 
        packDataToThread0(adr_sys_one_cluster_);
        packDataToThread0(ptcl_cluster_);
        packDataToThread0(id_ngb_multi_cluster);
        packDataToThread0(ptcl_outer);
        setNgbAdrHead(id_ngb_multi_cluster);
        n_pcluster_self_node_ = ptcl_cluster_[0].size();
        std::sort(ptcl_outer[0].getPointer(), ptcl_outer[0].getPointer(ptcl_outer[0].size()), OPLessID());

/*        if(PS::Comm::getRank() == 0 && ptcl_outer[0].size() > 0){
            PS::S32 n_tmp = 1;
            PS::S32 n_tmp_tmp = 0;
            //std::cerr<<"ptcl_outer[0].size()="<<ptcl_outer[0].size()<<std::endl;
            PS::S32 ref_tmp = ptcl_outer[0][0].id_;
            for(PS::S32 i=0; i<ptcl_outer[0].size(); i++){
                if( ref_tmp != ptcl_outer[0][i].id_){
                    ref_tmp = ptcl_outer[0][i].id_;
                    n_tmp++;
                    n_tmp_tmp = 0;
                }
                n_tmp_tmp++;
            }
            }*/
        mergePtclCluster(ptcl_outer, id_ngb_multi_cluster); // ptcl_cluster_ is complited
        id_to_adr_pcluster_.clear(); // temporally add 
        for(PS::S32 i=0; i<ptcl_cluster_[0].size(); i++){
            id_to_adr_pcluster_.insert(std::pair<PS::S32, PS::S32>(ptcl_cluster_[0][i].id_, i));
        }
        adr_ngb_multi_cluster_.clearSize();
#ifdef HARD_DEBUG
        assert(id_ngb_multi_cluster[0].size()<ARRAY_ALLOW_LIMIT);
#endif        
        adr_ngb_multi_cluster_.resizeNoInitialize(id_ngb_multi_cluster[0].size());
        for(PS::S32 i=0; i<id_ngb_multi_cluster[0].size(); i++){
            const PS::S32 adr_self = id_to_adr_pcluster_[id_ngb_multi_cluster[0][i].first];
            const PS::S32 adr_ngb  = id_to_adr_pcluster_[id_ngb_multi_cluster[0][i].second];
            adr_ngb_multi_cluster_[i] = std::pair<PS::S32, PS::S32>(adr_self, adr_ngb);
        }
    }

    template<class Tsys>
    void checkPtclCluster(const Tsys & sys){
        for(PS::S32 i=0; i<n_pcluster_self_node_; i++){
            if(ptcl_cluster_[0][i].id_ != sys[ptcl_cluster_[0][i].adr_sys_].id){
                std::cerr<<"i="<<i<<" ptcl_cluster_[0][i].adr_sys_="<<ptcl_cluster_[0][i].adr_sys_
                         <<" ptcl_cluster_[0][i].id_="<<ptcl_cluster_[0][i].id_
                         <<" sys[ptcl_cluster_[0][i].adr_sys_].id="<<sys[ptcl_cluster_[0][i].adr_sys_].id<<std::endl;
            }
            assert(ptcl_cluster_[0][i].id_ == sys[ptcl_cluster_[0][i].adr_sys_].id);
        }
        std::cerr<<"PASS: checkPtclCluster"<<std::endl;
    }


    void searchClusterLocal(){
        const PS::S32 my_rank = PS::Comm::getRank();
        const PS::S32 n_loc = ptcl_cluster_[0].size();

        adr_sys_multi_cluster_isolated_.clearSize();
        for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
            mediator_sorted_id_cluster_[i].clear();
        }
        mediator_sorted_id_cluster_.clearSize();
        static PS::ReallocatableArray<Cluster> cluster_isolated;
        for(PS::S32 i=0; i<cluster_isolated.size(); i++) cluster_isolated[i].clear();
        cluster_isolated.clearSize();
        n_ptcl_in_multi_cluster_isolated_.clearSize();
        n_ptcl_in_multi_cluster_isolated_offset_.clearSize();
        n_ptcl_in_multi_cluster_isolated_offset_.push_back(0);

        for(PS::S32 i=0; i<n_loc; i++){
            bool flag_isolated = true;
            if(ptcl_cluster_[0][i].flag_searched_ == false){
                PS::S32 n_ptcl_in_cluster = 0;
                PtclCluster * p_top = ptcl_cluster_[0].getPointer(i);
                searchClusterImpl(ptcl_cluster_[0].getPointer(i), ptcl_cluster_[0].getPointer(),
                                  p_top, n_ptcl_in_cluster,
                                  adr_ngb_multi_cluster_, flag_isolated);
                if(flag_isolated){
                    PtclCluster * p_tmp = ptcl_cluster_[0].getPointer(i);
                    PS::S32 id_cluster = p_tmp->id_;
                    PS::S32 adr_sys_multi_cluster_isolated_head = adr_sys_multi_cluster_isolated_.size();
                    while(p_tmp != NULL){
                        if(p_tmp->id_ < id_cluster) id_cluster = p_tmp->id_;
                        adr_sys_multi_cluster_isolated_.push_back(p_tmp->adr_sys_);
                        p_tmp = p_tmp->next_;
                    }
                    cluster_isolated.push_back( Cluster(id_cluster, n_ptcl_in_cluster, n_ptcl_in_cluster, adr_sys_multi_cluster_isolated_head, my_rank) );
                    n_ptcl_in_multi_cluster_isolated_.push_back(n_ptcl_in_cluster);
                    n_ptcl_in_multi_cluster_isolated_offset_.push_back(n_ptcl_in_cluster+n_ptcl_in_multi_cluster_isolated_offset_.back());
                }
                else{
                    PtclCluster * p_tmp = ptcl_cluster_[0].getPointer(i);
                    PS::S32 id_cluster = p_tmp->id_;
                    while(p_tmp != NULL){
                        if(p_tmp->id_ < id_cluster) id_cluster = p_tmp->id_;
                        PS::S32 adr_pcluster = p_tmp-ptcl_cluster_[0].getPointer();
                        mediator_sorted_id_cluster_.push_back(Mediator(p_tmp->id_, p_tmp->adr_sys_, adr_pcluster, -1, PS::Comm::getRank()) ); // temporarily send_rank_ is my_rank
                        p_tmp = p_tmp->next_;
                    }
                }
            }
        }
    }

    template<class Tsys>
    void checkMediator(const Tsys & sys){
        for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
            assert(mediator_sorted_id_cluster_[i].id_ == ptcl_cluster_[0][mediator_sorted_id_cluster_[i].adr_pcluster_].id_);
        }
        std::cerr<<"PASS: checkMediator (1st)"<<std::endl;
        for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
            if( mediator_sorted_id_cluster_[i].adr_sys_ == -1) continue;
            assert(mediator_sorted_id_cluster_[i].id_ == sys[mediator_sorted_id_cluster_[i].adr_sys_].id);
        }
        std::cerr<<"PASS: checkMediator (2nd)"<<std::endl;
    }

    void setIdClusterLocal(){
        for(PS::S32 i=0; i<ptcl_cluster_[0].size(); i++){
            ptcl_cluster_[0][i].setIdClusterImpl(ptcl_cluster_[0], adr_ngb_multi_cluster_);
        }
    }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL    
    template <class Ttree>
    void connectNodes(PS::F64ort pos_domain[], Ttree & tree){
        const PS::S32 n_proc_tot = PS::Comm::getNumberOfProc();
        const PS::S32 my_rank = PS::Comm::getRank();
        static PS::ReallocatableArray<PS::S32> rank_neighbor;
        static PS::ReallocatableArray<PS::S32> n_send;
        static PS::ReallocatableArray<PS::S32> n_recv;
        static PS::ReallocatableArray<PS::S32> n_send_disp;
        static PS::ReallocatableArray<PS::S32> n_recv_disp;
        static PS::ReallocatableArray<PS::S32> id_send;
        static PS::ReallocatableArray<PS::S32> id_recv;
        static PS::ReallocatableArray<MPI_Request> req_send;
        static PS::ReallocatableArray<MPI_Request> req_recv;
        static PS::ReallocatableArray<MPI_Status> stat_send;
        static PS::ReallocatableArray<MPI_Status> stat_recv;
        rank_neighbor.clearSize();
        n_send.resizeNoInitialize(n_proc_tot);
        n_recv.resizeNoInitialize(n_proc_tot);
        n_send_disp.resizeNoInitialize(n_proc_tot+1);
        n_recv_disp.resizeNoInitialize(n_proc_tot+1);
        req_send.resizeNoInitialize(n_proc_tot);
        req_recv.resizeNoInitialize(n_proc_tot);
        stat_send.resizeNoInitialize(n_proc_tot);
        stat_recv.resizeNoInitialize(n_proc_tot);
        SetRankComm(pos_domain, tree, rank_neighbor);
        for(PS::S32 i=0; i<n_proc_tot; i++){
            n_send[i] = n_recv[i] = 0;
        }
        for(PS::S32 i=n_pcluster_self_node_; i<ptcl_cluster_[0].size(); i++){
            assert(ptcl_cluster_[0][i].rank_org_ != my_rank);
            n_send[ptcl_cluster_[0][i].rank_org_]++;
        }
        for(PS::S32 i=0; i<rank_neighbor.size(); i++){
            PS::S32 rank = rank_neighbor[i];
            MPI_Isend(n_send.getPointer(rank),  1, PS::GetDataType<PS::S32>(),
                      rank, 1837, MPI_COMM_WORLD, req_send.getPointer(i));
            MPI_Irecv(n_recv.getPointer(rank),  1, PS::GetDataType<PS::S32>(),
                      rank, 1837, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        MPI_Waitall(rank_neighbor.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_neighbor.size(), req_recv.getPointer(), stat_recv.getPointer());
        n_send_disp[0] = n_recv_disp[0] = 0;
        for(PS::S32 i=0; i<n_proc_tot; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            n_send_disp[i+1] = n_send_disp[i] + n_send[i];
        }
        id_send.resizeNoInitialize(n_send_disp[n_proc_tot]);
        id_recv.resizeNoInitialize(n_recv_disp[n_proc_tot]);
        static PS::ReallocatableArray<PS::S32> n_cnt_send;
        n_cnt_send.resizeNoInitialize(n_proc_tot);
        for(PS::S32 i=0; i<n_proc_tot; i++) n_cnt_send[i] = 0;
        for(PS::S32 i=n_pcluster_self_node_; i<ptcl_cluster_[0].size(); i++){
            PS::S32 rank = ptcl_cluster_[0][i].rank_org_;
            PS::S32 adr = n_send_disp[rank] + n_cnt_send[rank];
            id_send[adr] = ptcl_cluster_[0][i].id_;
            n_cnt_send[rank]++;
        }
        rank_send_cluster_.resizeNoInitialize(0);
        rank_recv_cluster_.resizeNoInitialize(0);
        n_cluster_send_.resizeNoInitialize(0);
        n_cluster_recv_.resizeNoInitialize(0);
        PS::S32 n_proc_send = 0;
        PS::S32 n_proc_recv = 0;
        for(PS::S32 i=0; i<rank_neighbor.size(); i++){
            PS::S32 rank = rank_neighbor[i];
            if(n_send[rank] > 0){
                rank_recv_cluster_.push_back(rank); // NOTE: recv not send
                n_cluster_recv_.push_back(n_send[rank]); // NOTE: recv not send
                MPI_Isend(id_send.getPointer(n_send_disp[rank]), n_send[rank], PS::GetDataType<PS::S32>(),
                          rank, 1871, MPI_COMM_WORLD, req_send.getPointer(n_proc_recv));
                n_proc_recv++;
            }
            if(n_recv[rank] > 0){
                rank_send_cluster_.push_back(rank); // NOTE: send
                n_cluster_send_.push_back(n_recv[rank]); // NOTE: send
                MPI_Irecv(id_recv.getPointer(n_recv_disp[rank]), n_recv[rank], PS::GetDataType<PS::S32>(),
                          rank, 1871, MPI_COMM_WORLD, req_recv.getPointer(n_proc_send));
                n_proc_send++;
            }
        }
        MPI_Waitall(n_proc_recv, req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(n_proc_send, req_recv.getPointer(), stat_recv.getPointer());
        adr_pcluster_send_.clearSize();
        adr_pcluster_recv_.clearSize();
        for(PS::S32 i=0; i<n_recv_disp.back(); i++){
#ifdef CLUSTER_DEBUG
            auto itr = id_to_adr_pcluster_.find(id_recv[i]);
            assert(itr != id_to_adr_pcluster_.end());
            assert(id_to_adr_pcluster_[id_recv[i]] < n_pcluster_self_node_);
#endif
            adr_pcluster_send_.push_back(id_to_adr_pcluster_[id_recv[i]]);
        }
        for(PS::S32 i=0; i<n_send_disp.back(); i++){
#ifdef CLUSTER_DEBUG
            auto itr = id_to_adr_pcluster_.find(id_send[i]);
            assert(itr != id_to_adr_pcluster_.end());
            /*
              if(id_to_adr_pcluster_[id_send[i]] < n_pcluster_self_node_ && my_rank == 0){
              std::cout<<"my_rank= "<<my_rank
              <<" id_send[i]= "<<id_send[i]
              <<" id_to_adr_pcluster_[id_send[i]]= "<<id_to_adr_pcluster_[id_send[i]]
              <<" n_pcluster_self_node_= "<<n_pcluster_self_node_
              <<std::endl;
              }
            */
            assert(id_to_adr_pcluster_[id_send[i]] >= n_pcluster_self_node_);
#endif
            adr_pcluster_recv_.push_back(id_to_adr_pcluster_[id_send[i]]);
        }
        n_cluster_disp_send_.resizeNoInitialize(n_proc_send+1);
        n_cluster_disp_recv_.resizeNoInitialize(n_proc_recv+1);
        n_cluster_disp_send_[0] = n_cluster_disp_recv_[0] = 0;
        for(PS::S32 i=0; i<n_proc_send; i++) n_cluster_disp_send_[i+1] = n_cluster_disp_send_[i] + n_cluster_send_[i];
        for(PS::S32 i=0; i<n_proc_recv; i++) n_cluster_disp_recv_[i+1] = n_cluster_disp_recv_[i] + n_cluster_recv_[i];
    }

    void setIdClusterGlobalIteration(){
        //        const PS::S32 my_rank = PS::Comm::getRank();
        const PS::S32 n_proc_tot = PS::Comm::getNumberOfProc();
#ifdef HARD_DEBUG
        //assert(n_cluster_disp_send_.back()<ARRAY_ALLOW_LIMIT);
        if(n_cluster_disp_send_.back()>ARRAY_ALLOW_LIMIT) {
            std::cerr<<"Error: size overflow: rank: "<<PS::Comm::getRank()<<" n_cluster_disp_send_.back()="<<n_cluster_disp_send_.back()<<" size="<<n_cluster_disp_send_.size()<<std::endl;
        }
        if(n_cluster_disp_recv_.back()>ARRAY_ALLOW_LIMIT) {
            std::cerr<<"Error: size overflow: rank: "<<PS::Comm::getRank()<<" n_cluster_disp_recv_.back()="<<n_cluster_disp_recv_.back()<<" size="<<n_cluster_disp_recv_.size()<<std::endl;
        }
#endif        
        id_cluster_send_.resizeNoInitialize(n_cluster_disp_send_.back());
        id_cluster_recv_.resizeNoInitialize(n_cluster_disp_recv_.back());
        static PS::ReallocatableArray<MPI_Request> req_send;
        static PS::ReallocatableArray<MPI_Request> req_recv;
        static PS::ReallocatableArray<MPI_Status> stat_send;
        static PS::ReallocatableArray<MPI_Status> stat_recv;
        req_send.reserve(n_proc_tot);
        req_recv.reserve(n_proc_tot);
        stat_send.reserve(n_proc_tot);
        stat_recv.reserve(n_proc_tot);
        bool flag_itr_glb = true;
        PS::S32 n_loop = 0;
        while(flag_itr_glb){
            flag_itr_glb = false;
            //if(PS::Comm::getRank() == 0) std::cout<<"n_loop="<<n_loop<<std::endl;
            //id_cluster_send_.resizeNoInitialize(0);
            PS::S32 n_proc_send = 0;
            PS::S32 n_proc_recv = 0;
            for(PS::S32 i=0; i<adr_pcluster_send_.size(); i++){
                PS::S32 adr = adr_pcluster_send_[i];
                assert(adr <= ptcl_cluster_[0].size());
                id_cluster_send_[i] = ptcl_cluster_[0][adr].id_cluster_;
            }
            for(PS::S32 i=0; i<rank_send_cluster_.size(); i++){
                PS::S32 rank = rank_send_cluster_[i];
                MPI_Isend(id_cluster_send_.getPointer(n_cluster_disp_send_[i]), n_cluster_send_[i], PS::GetDataType<PS::S32>(),
                          rank, 1947, MPI_COMM_WORLD, req_send.getPointer(n_proc_send));
                n_proc_send++;
            }
            for(PS::S32 i=0; i<rank_recv_cluster_.size(); i++){
                PS::S32 rank = rank_recv_cluster_[i];
                MPI_Irecv(id_cluster_recv_.getPointer(n_cluster_disp_recv_[i]), n_cluster_recv_[i], PS::GetDataType<PS::S32>(),
                          rank, 1947, MPI_COMM_WORLD, req_recv.getPointer(n_proc_recv));
                n_proc_recv++;
            }
            MPI_Waitall(n_proc_send, req_send.getPointer(), stat_send.getPointer());
            MPI_Waitall(n_proc_recv, req_recv.getPointer(), stat_recv.getPointer());

            bool flag_itr_loc = false;
            for(PS::S32 i=0; i<id_cluster_recv_.size(); i++){
                PS::S32 adr = adr_pcluster_recv_[i];
                if(id_cluster_recv_[i] < ptcl_cluster_[0][adr].id_cluster_){
                    ptcl_cluster_[0][adr].id_cluster_ = id_cluster_recv_[i];
                    flag_itr_loc |= true;
                }
            }
            for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
                PS::S32 adr = mediator_sorted_id_cluster_[i].adr_pcluster_;
                PS::S32 id_clucster_prev = ptcl_cluster_[0][adr].id_cluster_;
                ptcl_cluster_[0][adr].setIdClusterImpl(ptcl_cluster_[0], adr_ngb_multi_cluster_);
                if( id_clucster_prev != ptcl_cluster_[0][adr].id_cluster_){
                    flag_itr_loc |= true;
                }
            }
            flag_itr_glb = PS::Comm::synchronizeConditionalBranchOR(flag_itr_loc);
            n_loop++;
        }
        for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
            PS::S32 adr = mediator_sorted_id_cluster_[i].adr_pcluster_;
            mediator_sorted_id_cluster_[i].id_cluster_ = ptcl_cluster_[0][adr].id_cluster_;
        }
        std::sort(mediator_sorted_id_cluster_.getPointer(0), mediator_sorted_id_cluster_.getPointer(mediator_sorted_id_cluster_.size()), OPLessIDCluster());
    }

    template<class Tsys>
    void sendAndRecvCluster(const Tsys & sys){
        PS::S32 my_rank = PS::Comm::getRank();
        PS::S32 n_proc = PS::Comm::getNumberOfProc();
        static PS::ReallocatableArray<Cluster> cluster_loc;
        cluster_loc.clearSize();
        if(mediator_sorted_id_cluster_.size() > 0){
            PS::S32 id_cluster_ref = mediator_sorted_id_cluster_[0].id_cluster_;
            cluster_loc.push_back( Cluster(id_cluster_ref, 0, 0, 0, my_rank) );
            for(PS::S32 i=0; i<mediator_sorted_id_cluster_.size(); i++){
                if( id_cluster_ref != mediator_sorted_id_cluster_[i].id_cluster_){
                    id_cluster_ref = mediator_sorted_id_cluster_[i].id_cluster_;
                    cluster_loc.push_back( Cluster(id_cluster_ref, 0, 0, i, my_rank) );
                }
                if(mediator_sorted_id_cluster_[i].adr_sys_>=0) cluster_loc.back().n_ptcl_stored_++;
                cluster_loc.back().n_ptcl_++;
            }
        }
        //std::cerr<<"sendAndRecvCluster 0: "<<my_rank<<std::endl;

        static PS::ReallocatableArray<PS::S32> n_cluster_recv;
        static PS::ReallocatableArray<PS::S32> n_cluster_recv_disp;
        n_cluster_recv.resizeNoInitialize(n_proc);
        n_cluster_recv_disp.resizeNoInitialize(n_proc+1);
        PS::S32 n_cluster_tot_loc = cluster_loc.size();
        PS::Comm::allGather(&n_cluster_tot_loc, 1, n_cluster_recv.getPointer());
        n_cluster_recv_disp[0] = 0;
        for(PS::S32 i=0; i<n_proc; i++){
            n_cluster_recv_disp[i+1] = n_cluster_recv[i] + n_cluster_recv_disp[i];
        }

        static PS::ReallocatableArray<Cluster> cluster_recv;
        cluster_recv.resizeNoInitialize(n_cluster_recv_disp[n_proc]);
        if(n_cluster_tot_loc > 0){
            PS::Comm::allGatherV(cluster_loc.getPointer(), n_cluster_tot_loc,
                                 cluster_recv.getPointer(),
                                 n_cluster_recv.getPointer(), n_cluster_recv_disp.getPointer());
        }
        else{
            Cluster tmp;
            PS::Comm::allGatherV(&tmp, n_cluster_tot_loc,
                                 cluster_recv.getPointer(),
                                 n_cluster_recv.getPointer(), n_cluster_recv_disp.getPointer());
        }

        //std::cerr<<"sendAndRecvCluster 1: "<<my_rank<<std::endl;

        // exchange cluster info
        //////////////////////////

        for(PS::S32 i0=0; i0<n_cluster_tot_loc; i0++){
            const PS::S32 id_cluster = cluster_loc[i0].id_;
            const PS::S32 n_ptcl_in_cluster = cluster_loc[i0].n_ptcl_stored_;
            PS::S32 rank_send_tmp = my_rank;
            PS::S32 n_ptcl_max = n_ptcl_in_cluster;
            for(PS::S32 i1=0; i1<n_proc; i1++){
                if(i1 == my_rank) continue;
                /*
                  if(my_rank == 0){
                  std::cout<<"i1="<<i1<<std::endl;
                  std::cout<<"n_cluster_recv_disp[i1]="<<n_cluster_recv_disp[i1]
                  <<" n_cluster_recv_disp[i1+1]="<<n_cluster_recv_disp[i1+1]<<std::endl;
                  }
                */
                for(PS::S32 i2=n_cluster_recv_disp[i1]; i2<n_cluster_recv_disp[i1+1]; i2++){
                    if(id_cluster == cluster_recv[i2].id_){
                        /*
                          if(my_rank == 0){
                          std::cout<<"cluster_recv[i2].id_="<<cluster_recv[i2].id_<<std::endl;
                          std::cout<<"n_ptcl_in_cluster="<<n_ptcl_in_cluster<<std::endl;
                          std::cout<<"cluster_recv[i2].n_ptcl_stored_="<<cluster_recv[i2].n_ptcl_stored_<<std::endl;
                          }
                        */
                        if( (n_ptcl_max < cluster_recv[i2].n_ptcl_stored_) ||
                            (n_ptcl_max == cluster_recv[i2].n_ptcl_stored_ && rank_send_tmp > i1) ){
                            rank_send_tmp = i1;
                            n_ptcl_max = cluster_recv[i2].n_ptcl_stored_;
                        }
                    }
                }
            }
            cluster_loc[i0].rank_ = rank_send_tmp;
        }
        std::sort(cluster_loc.getPointer(), cluster_loc.getPointer(cluster_loc.size()), OPLessRank());
        //std::cerr<<"sendAndRecvCluster 2: "<<my_rank<<std::endl;

        ////////////
        // pack and send particles
        ptcl_send_.clearSize();
        rank_send_ptcl_.clearSize();
        n_ptcl_send_.clearSize();
        adr_sys_ptcl_send_.clearSize();
        if(cluster_loc.size() > 0){
            PS::S32 rank_send_ref = -999999;
            for(PS::S32 i=0; i<cluster_loc.size(); i++){
                if(cluster_loc[i].rank_ == my_rank) continue;
                if( rank_send_ref != cluster_loc[i].rank_){
                    rank_send_ref = cluster_loc[i].rank_;
                    rank_send_ptcl_.push_back(rank_send_ref);
                    n_ptcl_send_.push_back(0);
                }

                PS::S32 n_tmp = cluster_loc[i].n_ptcl_;
                PS::S32 n_cnt = 0;
                for(PS::S32 ii=0; ii<n_tmp; ii++){
                    PS::S32 adr_sys = mediator_sorted_id_cluster_[cluster_loc[i].adr_head_+ii].adr_sys_;
                    if(adr_sys < 0){
                        continue;
                    }
                    mediator_sorted_id_cluster_[cluster_loc[i].adr_head_+ii].rank_send_ = rank_send_ref; // 2006.09.06
                    const auto &p = sys[adr_sys];
                    adr_sys_ptcl_send_.push_back(adr_sys);
                    ptcl_send_.push_back(PtclComm(p));
                    ptcl_send_.back().id_cluster = cluster_loc[i].id_;
#ifdef HARD_DEBUG
                    if(ptcl_send_.back().id==0&&ptcl_send_.back().group_data.artificial.isUnused()) {
                        std::cerr<<"Error! sending particle is unused! adr="<<adr_sys<<std::endl;
                        abort();
                    }
#endif
                    n_ptcl_send_.back()++;
                    n_cnt++;
                }
                /*
                  if(cluster_loc[i].n_ptcl_stored_ != n_cnt){
                  std::cerr<<"cluster_loc[i].n_ptcl_stored_="<<cluster_loc[i].n_ptcl_stored_
                  <<" n_cnt="<<n_cnt
                  <<std::endl;
                  }
                */
                assert(cluster_loc[i].n_ptcl_stored_ == n_cnt);
            }
        }

        //std::cerr<<"sendAndRecvCluster 3: "<<my_rank<<std::endl;
	
        //static PS::ReallocatableArray<PS::S32> n_ptcl_disp_send;
        n_ptcl_disp_send_.resizeNoInitialize(rank_send_ptcl_.size()+1);
        n_ptcl_disp_send_[0] = 0;
        for(PS::S32 i=0; i<rank_send_ptcl_.size(); i++){
            n_ptcl_disp_send_[i+1] = n_ptcl_disp_send_[i] + n_ptcl_send_[i];
        }

        static PS::ReallocatableArray<MPI_Request> req_send;
        static PS::ReallocatableArray<MPI_Status> stat_send;
        req_send.resizeNoInitialize(n_proc);
        stat_send.resizeNoInitialize(n_proc);
        for(PS::S32 i=0; i<rank_send_ptcl_.size(); i++){
            PS::S32 rank = rank_send_ptcl_[i];
            MPI_Isend(ptcl_send_.getPointer(n_ptcl_disp_send_[i]),  n_ptcl_send_[i], 
                      PS::GetDataType<PtclComm>(),
                      rank, 2136, MPI_COMM_WORLD, req_send.getPointer(i));
        }

        //std::cerr<<"sendAndRecvCluster 4: "<<my_rank<<std::endl;
        // pack and send particles
        ////////////
	
        ////////////
        // make and recv particles
        rank_recv_ptcl_.clearSize();
        n_ptcl_recv_.clearSize();
        for(PS::S32 i0=0; i0<n_proc; i0++){
            if(i0 == my_rank) continue;
            bool flag_recv = false;
            n_ptcl_recv_.push_back(0);
            rank_recv_ptcl_.push_back(i0);
            for(PS::S32 i1=n_cluster_recv_disp[i0]; i1<n_cluster_recv_disp[i0+1]; i1++){
                PS::S32 id_cluster = cluster_recv[i1].id_;
                for(PS::S32 i2=0; i2<cluster_loc.size(); i2++){
                    if(id_cluster == cluster_loc[i2].id_ && cluster_loc[i2].rank_ == my_rank){
                        n_ptcl_recv_.back() += cluster_recv[i1].n_ptcl_stored_;
                        flag_recv = true;
                    }
                }
            }
            if(!flag_recv){
                n_ptcl_recv_.resizeNoInitialize(n_ptcl_recv_.size()-1);
                rank_recv_ptcl_.resizeNoInitialize(rank_recv_ptcl_.size()-1);
            }
        }

        //static PS::ReallocatableArray<PS::S32> n_ptcl_disp_recv_;
        n_ptcl_disp_recv_.resizeNoInitialize(n_ptcl_recv_.size()+1);
        n_ptcl_disp_recv_[0] = 0;
        for(PS::S32 i=0; i<n_ptcl_recv_.size(); i++){
            n_ptcl_disp_recv_[i+1] = n_ptcl_disp_recv_[i] + n_ptcl_recv_[i];
        }
        //std::cerr<<"sendAndRecvCluster 4.5: "<<my_rank<<std::endl;
	
        static PS::ReallocatableArray<MPI_Request> req_recv;
        static PS::ReallocatableArray<MPI_Status> stat_recv;
        ptcl_recv_.resizeNoInitialize(n_ptcl_disp_recv_[n_ptcl_recv_.size()]);
        req_recv.resizeNoInitialize(n_proc);
        stat_recv.resizeNoInitialize(n_proc);
        for(PS::S32 i=0; i<rank_recv_ptcl_.size(); i++){
            PS::S32 rank = rank_recv_ptcl_[i];
            MPI_Irecv(ptcl_recv_.getPointer(n_ptcl_disp_recv_[i]), n_ptcl_recv_[i], 
                      PS::GetDataType<PtclComm>(),
                      rank, 2136, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        //std::cerr<<"sendAndRecvCluster 4.7: "<<my_rank<<std::endl;
        MPI_Waitall(rank_send_ptcl_.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_recv_ptcl_.size(), req_recv.getPointer(), stat_recv.getPointer());

        //std::cerr<<"sendAndRecvCluster 5: "<<my_rank<<std::endl;

        // make and recv particles
        ////////////
    }

    //! Send and receive the remote particles due to the change of kick
    /* Send local single particles to remote nodes and receive remote single particles.
       Notice _ptcl_hard are not overlap with ptcl_send
       @param[in,out] _sys: particle system. Notice the local particles of non-group members are updated.
       @param[in,out] _ptcl_hard: local partical in system_hard_connected
     */
    template<class Tsys, class Tphard>
    void SendSinglePtcl(Tsys & _sys,
                        PS::ReallocatableArray<Tphard> & _ptcl_hard){
        //// First, write back local kicked single particles to sys.
        //for(PS::S32 i=0; i<n; i++) {
        //    const PS::S32 adr = _ptcl_hard[i].adr_org;
        //    if(adr >=0){
//#ifdef HARD_DEBUG
        //        assert( _sys[adr].id == _ptcl_hard[i].id );
//#endif  // 
        //        _sys[adr].DataCopy(_ptcl_hard[i]);
        //    }
        //}
        // write kicked single particle to sending buffer
        for(PS::S32 i=0; i<ptcl_send_.size(); i++){
            PS::S32 adr = adr_sys_ptcl_send_[i];
#ifdef HARD_DEBUG
            assert(_sys[adr].id == ptcl_send_[i].id);
#endif
            // write local single particle
            if(ptcl_send_[i].group_data.artificial.isSingle()) {
                ptcl_send_[i].DataCopy(_sys[adr]);
            }
        }
        static PS::ReallocatableArray<MPI_Request> req_send;
        static PS::ReallocatableArray<MPI_Status> stat_send;
        static PS::ReallocatableArray<MPI_Request> req_recv;
        static PS::ReallocatableArray<MPI_Status> stat_recv;
        req_send.resizeNoInitialize(rank_send_ptcl_.size());
        stat_send.resizeNoInitialize(rank_send_ptcl_.size());
        req_recv.resizeNoInitialize(rank_recv_ptcl_.size());
        stat_recv.resizeNoInitialize(rank_recv_ptcl_.size());

        for(PS::S32 i=0; i<rank_send_ptcl_.size(); i++){
            PS::S32 rank = rank_send_ptcl_[i];
            MPI_Isend(ptcl_send_.getPointer(n_ptcl_disp_send_[i]),  n_ptcl_send_[i],
                      PS::GetDataType<PtclComm>(),
                      rank, 2239, MPI_COMM_WORLD, req_send.getPointer(i));
        }
        for(PS::S32 i=0; i<rank_recv_ptcl_.size(); i++){
            PS::S32 rank = rank_recv_ptcl_[i];
            MPI_Irecv(ptcl_recv_.getPointer(n_ptcl_disp_recv_[i]), n_ptcl_recv_[i],
                      PS::GetDataType<PtclComm>(),
                      rank, 2239, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        MPI_Waitall(rank_send_ptcl_.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_recv_ptcl_.size(), req_recv.getPointer(), stat_recv.getPointer());

        // Receive remote single particle data
        const PS::S32 n = _ptcl_hard.size();
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 adr = _ptcl_hard[i].adr_org;
            if(adr <0 && _ptcl_hard[i].group_data.artificial.isSingle()){
#ifdef HARD_DEBUG
                assert( ptcl_recv_[-(adr+1)].id == _ptcl_hard[i].id );
                assert( ptcl_recv_[-(adr+1)].group_data.artificial.isSingle());
#endif
                _ptcl_hard[i].DataCopy(ptcl_recv_[-(adr+1)]);
            }
        }

    }
    
    //! send and receive particles on remote
    /* @param[in,out] _sys: particle system
       @param[in] _ptcl_hard: local partical in system_hard_connected
       @param[out] _removelist: address on _sys of particles that are need to be removed are stored.
     */
    template<class Tsys, class Tphard>
    void writeAndSendBackPtcl(Tsys & _sys,
                              const PS::ReallocatableArray<Tphard> & _ptcl_hard,
                              PS::ReallocatableArray<PS::S32> &_removelist){
        //  const PS::S32 my_rank = PS::Comm::getRank();
        //  const PS::S32 n_proc  = PS::Comm::getNumberOfProc();
        const PS::S32 n = _ptcl_hard.size();
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 adr = _ptcl_hard[i].adr_org;
            if( adr >= 0){
#ifdef HARD_DEBUG
                assert( _sys[adr].id == _ptcl_hard[i].id);
#endif
                _sys[adr].DataCopy(_ptcl_hard[i]);
                if(_sys[adr].id==0.0&&_sys[adr].group_data.artificial.isUnused()) _removelist.push_back(adr);
            }
            else{
                //assert( ptcl_recv_[-(adr+1)].id == _ptcl_hard[i].id );
                ptcl_recv_[-(adr+1)].DataCopy(_ptcl_hard[i]);
            }
        }
        static PS::ReallocatableArray<MPI_Request> req_recv;
        static PS::ReallocatableArray<MPI_Status> stat_recv;
        static PS::ReallocatableArray<MPI_Request> req_send;
        static PS::ReallocatableArray<MPI_Status> stat_send;
        req_recv.resizeNoInitialize(rank_recv_ptcl_.size());
        stat_recv.resizeNoInitialize(rank_recv_ptcl_.size());
        req_send.resizeNoInitialize(rank_send_ptcl_.size());
        stat_send.resizeNoInitialize(rank_send_ptcl_.size());

        for(PS::S32 i=0; i<rank_recv_ptcl_.size(); i++){
            PS::S32 rank = rank_recv_ptcl_[i];
            MPI_Isend(ptcl_recv_.getPointer(n_ptcl_disp_recv_[i]), n_ptcl_recv_[i],
                      PS::GetDataType<PtclComm>(),
                      rank, 2303, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        for(PS::S32 i=0; i<rank_send_ptcl_.size(); i++){
            PS::S32 rank = rank_send_ptcl_[i];
            MPI_Irecv(ptcl_send_.getPointer(n_ptcl_disp_send_[i]),  n_ptcl_send_[i],
                      PS::GetDataType<PtclComm>(),
                      rank, 2303, MPI_COMM_WORLD, req_send.getPointer(i));
        }
        MPI_Waitall(rank_send_ptcl_.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_recv_ptcl_.size(), req_recv.getPointer(), stat_recv.getPointer());

        for(PS::S32 i=0; i<ptcl_send_.size(); i++){
            PS::S32 adr = adr_sys_ptcl_send_[i];
#ifdef HARD_DEBUG
            assert(_sys[adr].id == ptcl_send_[i].id);
#endif
            _sys[adr].DataCopy(ptcl_send_[i]);
            if(_sys[adr].id==0&&_sys[adr].group_data.artificial.isUnused()) _removelist.push_back(adr);
        }
        // remove empty particles cannot do here
        // _sys.removeParticle(_removelist.getPointer(), _removelist.size());
    }


#endif


    const PS::ReallocatableArray<PS::S32> & getAdrSysOneCluster(){
        return adr_sys_one_cluster_[0];
    }
    
    //! get the address list of ptcl need to send in connected clusters
    const PS::ReallocatableArray<PS::S32> & getAdrSysConnectClusterSend(){
        return adr_sys_ptcl_send_;
    }

};


