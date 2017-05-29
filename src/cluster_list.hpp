#pragma once
#include<particle_simulator.hpp>
//#include"hard_force.hpp"

template <class Ttree>
void SetRankComm(const PS::F64ort pos_domain[], Ttree &tree,
                 PS::ReallocatableArray<PS::S32> & rank_neighbor,
                 const PS::F64 eps_sq = 0.0){
    rank_neighbor.clearSize();
//    static const PS::F64 SAFTY_FACTOR_FOR_SEARCH_SQ = 1.01;
    const PS::S32 my_rank = PS::Comm::getRank();
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
//    const PS::F64ort pos_my_domain = pos_domain[my_rank];
    const PS::F64ort pos_my_domain = tree.getOuterBoundaryOfLocalTree();
    const PS::F64ort pos_my_domain_in = tree.getInnerBoundaryOfLocalTree();
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
}

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
    PtclCluster(): id_(-1), adr_sys_(-1), adr_ngb_head_(-1), n_ngb_(0), 
                   flag_searched_(false), next_(NULL), rank_org_(-1), id_cluster_(id_){}
    PtclCluster(const PS::S32 _id,    const PS::S32 _adr_sys,    const PS::S32 _adr_ngb_head, 
                const PS::S32 _n_ngb, const bool _flag_searched, PtclCluster * const _next, 
                PS::S32 _rank_org):
        id_(_id), adr_sys_(_adr_sys), adr_ngb_head_(_adr_ngb_head),  
        n_ngb_(_n_ngb), flag_searched_(_flag_searched), next_(_next), rank_org_(_rank_org), 
        id_cluster_(_id){}
    void clear(){
        id_ = -1; adr_sys_ = -1; adr_ngb_head_ = -1; n_ngb_ = 0; 
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
    PtclOuter(): id_(-1), id_ngb_(-1), rank_org_(-1){}
    PtclOuter(const PS::S32 _id, const PS::S32 _id_ngb, const PS::S32 _rank_org): 
        id_(_id), id_ngb_(_id_ngb), rank_org_(_rank_org){}
    void clear(){
        id_ = -1; id_ngb_ = -1; rank_org_ = -1;
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
    Mediator():id_(-1), adr_sys_(-1), adr_pcluster_(-1), id_cluster_(-1), rank_send_(-1){} 
    Mediator(const PS::S32 _id, const PS::S32 _adr_sys, 
             const PS::S32 _adr_pcluster, const PS::S32 _id_cluster,
             const PS::S32 _rank_send):
        id_(_id), adr_sys_(_adr_sys), adr_pcluster_(_adr_pcluster), id_cluster_(_id_cluster),
        rank_send_(_rank_send){}
    void clear(){
        id_ = -1; adr_sys_ = -1; adr_pcluster_ = -1; id_cluster_ = -1; rank_send_ = -1;
    }
    void dump(){
        std::cout<<" id="<<id_<<" adr_sys="<<adr_sys_
                 <<" adr_pcluster_="<<adr_pcluster_<<" id_cluster="<<id_cluster_
                 <<" rank_send_="<<rank_send_<<std::endl;
    }
};

class PtclComm{
public:
    PS::S32 id_;
    PS::F64 mass_;
    PS::F64vec pos_;
    PS::F64vec vel_;
    PS::F64 r_out_;
    PS::S32 n_ngb_;
    PS::S32 id_cluster_;
    PtclComm(): id_(-1), mass_(0.0), pos_(0.0), vel_(0.0), r_out_(0.0), n_ngb_(0) {}
    PtclComm(const PS::S32 _id, const PS::F64 _mass, const PS::F64vec _pos, const PS::F64vec _vel, const PS::F64 _r_out, const PS::S32 _n_ngb):
        id_(_id), mass_(_mass), pos_(_pos), vel_(_vel), r_out_(_r_out), n_ngb_(_n_ngb) {}
    void dump(){
        std::cout<<" id="<<id_<<" mass="<<mass_
                 <<" pos="<<pos_<<" vel="<<vel_
                 <<" id_cluster_="<<id_cluster_
                 <<" r_out_="<<r_out_
                 <<" n_ngb_="<<n_ngb_
                 <<std::endl;
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


    // * adr_sys_one_cluster_
    // * ptcl_cluster_ 
    // * adr_ngb_multi_cluster_ 
    // * id_to_adr_pcluster_
    template<class Tsys, class Ttree, class Tepj>
    void searchNeighborAndCalcHardForceOMP(Tsys & sys,
                                           Ttree & tree,
                                           const PS::F64 r_out,
                                           const PS::F64 r_in,
                                           const PS::F64ort pos_domain[],
                                           const PS::F64 eps_sq=0.0){
        static PS::ReallocatableArray<PtclOuter> * ptcl_outer;
        static PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > * id_ngb_multi_cluster;
        const PS::S32 n_thread = PS::Comm::getNumberOfThread();
        ptcl_outer = new PS::ReallocatableArray<PtclOuter>[n_thread];
        id_ngb_multi_cluster = new PS::ReallocatableArray< std::pair<PS::S32, PS::S32> >[n_thread];
        const PS::S32 my_rank = PS::Comm::getRank();
        //        const PS::S32 n_proc_tot = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc = sys.getNumberOfParticleLocal();
#pragma omp parallel
        {
            const PS::S32 ith = PS::Comm::getThreadNum();
            adr_sys_one_cluster_[ith].clearSize();
            id_ngb_multi_cluster[ith].clearSize();
            for(PS::S32 i=0; i<ptcl_cluster_[ith].size(); i++) ptcl_cluster_[ith][i].clear();
            ptcl_cluster_[ith].clearSize();
            for(PS::S32 i=0; i<ptcl_outer[ith].size(); i++) ptcl_outer[ith][i].clear();
            ptcl_outer[ith].clearSize();
#pragma omp for
            for(PS::S32 i=0; i<n_loc; i++){
                Tepj * nbl = NULL;
                sys[i].n_ngb = tree.getNeighborListOneParticle(sys[i], nbl) - 1;
                assert(sys[i].n_ngb >= 0);
                if(sys[i].n_ngb == 0){
                    // no neighbor
                    adr_sys_one_cluster_[ith].push_back(i);
                }
                else{
                    // has neighbor
                    ptcl_cluster_[ith].push_back( PtclCluster(sys[i].id, i, id_ngb_multi_cluster[ith].size(), sys[i].n_ngb, false, NULL, my_rank) );
                    //PS::S32 n_tmp2 = 0;
                    //PS::S32 n_tmp3 = 0;
                    //PS::S32 n_tmp4 = 0;
                    for(PS::S32 ii=0; ii<sys[i].n_ngb+1; ii++){
                        if( (nbl+ii)->id == sys[i].id ){
                            //n_tmp2++;
                            continue;
                        }
                        id_ngb_multi_cluster[ith].push_back( std::pair<PS::S32, PS::S32>(sys[i].id, (nbl+ii)->id) );
                        CalcAccPotShortWithLinearCutoff
                            (sys[i].pos,    sys[i].acc,      sys[i].pot_tot,
                             (nbl+ii)->pos, (nbl+ii)->mass,  eps_sq,
                             r_out,         r_in);
                        if( (nbl+ii)->rank_org != my_rank ){
                            ptcl_outer[ith].push_back(PtclOuter((nbl+ii)->id, sys[i].id, (nbl+ii)->rank_org));
                            //n_tmp3++;
                        }
                    }
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

        if(PS::Comm::getRank() == 0 && ptcl_outer[0].size() > 0){
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
        }
        mergePtclCluster(ptcl_outer, id_ngb_multi_cluster); // ptcl_cluster_ is complited
        id_to_adr_pcluster_.clear(); // temporally add 
        for(PS::S32 i=0; i<ptcl_cluster_[0].size(); i++){
            id_to_adr_pcluster_.insert(std::pair<PS::S32, PS::S32>(ptcl_cluster_[0][i].id_, i));
        }
        adr_ngb_multi_cluster_.clearSize();
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

    template <class Ttree>
    void conectNodes(PS::F64ort pos_domain[], Ttree & tree){
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
                      rank, 0, MPI_COMM_WORLD, req_send.getPointer(i));
            MPI_Irecv(n_recv.getPointer(rank),  1, PS::GetDataType<PS::S32>(),
                      rank, 0, MPI_COMM_WORLD, req_recv.getPointer(i));
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
                          rank, 0, MPI_COMM_WORLD, req_send.getPointer(n_proc_send));
                n_proc_send++;
            }
            if(n_recv[rank] > 0){
                rank_send_cluster_.push_back(rank); // NOTE: send
                n_cluster_send_.push_back(n_recv[rank]); // NOTE: send
                MPI_Irecv(id_recv.getPointer(n_recv_disp[rank]), n_recv[rank], PS::GetDataType<PS::S32>(),
                          rank, 0, MPI_COMM_WORLD, req_recv.getPointer(n_proc_recv));
                n_proc_recv++;
            }
        }
        MPI_Waitall(n_proc_send, req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(n_proc_recv, req_recv.getPointer(), stat_recv.getPointer());
        adr_pcluster_send_.clearSize();
        adr_pcluster_recv_.clearSize();
        for(PS::S32 i=0; i<n_recv_disp.back(); i++){
            auto itr = id_to_adr_pcluster_.find(id_recv[i]);
            assert(itr != id_to_adr_pcluster_.end());
            assert(id_to_adr_pcluster_[id_recv[i]] < n_pcluster_self_node_);
            adr_pcluster_send_.push_back(id_to_adr_pcluster_[id_recv[i]]);
        }
        for(PS::S32 i=0; i<n_send_disp.back(); i++){
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
                          rank, 0, MPI_COMM_WORLD, req_send.getPointer(n_proc_send));
                n_proc_send++;
            }
            for(PS::S32 i=0; i<rank_recv_cluster_.size(); i++){
                PS::S32 rank = rank_recv_cluster_[i];
                MPI_Irecv(id_cluster_recv_.getPointer(n_cluster_disp_recv_[i]), n_cluster_recv_[i], PS::GetDataType<PS::S32>(),
                          rank, 0, MPI_COMM_WORLD, req_recv.getPointer(n_proc_recv));
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
                    ptcl_send_.push_back(PtclComm(p.id, p.mass, p.pos, p.vel, p.r_out, p.n_ngb));
                    ptcl_send_.back().id_cluster_ = cluster_loc[i].id_;
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
                      rank, 0, MPI_COMM_WORLD, req_send.getPointer(i));
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
                      rank, 0, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        //std::cerr<<"sendAndRecvCluster 4.7: "<<my_rank<<std::endl;
        MPI_Waitall(rank_send_ptcl_.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_recv_ptcl_.size(), req_recv.getPointer(), stat_recv.getPointer());

        //std::cerr<<"sendAndRecvCluster 5: "<<my_rank<<std::endl;

        // make and recv particles
        ////////////
    }

    template<class Tsys, class Tphard>
    void writeAndSendBackPtcl(Tsys & sys,
                              const PS::ReallocatableArray<Tphard> & ptcl_hard){
        //  const PS::S32 my_rank = PS::Comm::getRank();
        //  const PS::S32 n_proc  = PS::Comm::getNumberOfProc();
        const PS::S32 n = ptcl_hard.size();
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 adr = ptcl_hard[i].adr_org;
            if( adr >= 0){
                assert( sys[adr].id == ptcl_hard[i].id);
                sys[adr].id       = ptcl_hard[i].id;
                sys[adr].mass     = ptcl_hard[i].mass;
                sys[adr].pos      = ptcl_hard[i].pos;
                sys[adr].vel      = ptcl_hard[i].vel;
                sys[adr].r_out    = ptcl_hard[i].r_out;
            }
            else{
                assert( ptcl_recv_[-(adr+1)].id_ == ptcl_hard[i].id );
                ptcl_recv_[-(adr+1)].id_       = ptcl_hard[i].id;
                ptcl_recv_[-(adr+1)].mass_     = ptcl_hard[i].mass;
                ptcl_recv_[-(adr+1)].pos_      = ptcl_hard[i].pos;
                ptcl_recv_[-(adr+1)].vel_      = ptcl_hard[i].vel;
                ptcl_recv_[-(adr+1)].r_out_    = ptcl_hard[i].r_out;
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
                      rank, 0, MPI_COMM_WORLD, req_recv.getPointer(i));
        }
        for(PS::S32 i=0; i<rank_send_ptcl_.size(); i++){
            PS::S32 rank = rank_send_ptcl_[i];
            MPI_Irecv(ptcl_send_.getPointer(n_ptcl_disp_send_[i]),  n_ptcl_send_[i],
                      PS::GetDataType<PtclComm>(),
                      rank, 0, MPI_COMM_WORLD, req_send.getPointer(i));
        }
        MPI_Waitall(rank_send_ptcl_.size(), req_send.getPointer(), stat_send.getPointer());
        MPI_Waitall(rank_recv_ptcl_.size(), req_recv.getPointer(), stat_recv.getPointer());

        for(PS::S32 i=0; i<ptcl_send_.size(); i++){
            PS::S32 adr = adr_sys_ptcl_send_[i];
            assert(sys[adr].id == ptcl_send_[i].id_);
            sys[adr].mass     = ptcl_send_[i].mass_;
            sys[adr].pos      = ptcl_send_[i].pos_;
            sys[adr].vel      = ptcl_send_[i].vel_;
            sys[adr].r_out    = ptcl_send_[i].r_out_;
        }
    }


    const PS::ReallocatableArray<PS::S32> & getAdrSysOneCluster(){
        return adr_sys_one_cluster_[0];
    }


};


template<class Tptcl>
class SearchClusterHard{
private:
    ReallocatableArray<Tptcl> ptcl_;               ///new ptcl with c.m.
    ReallocatableArray<PS::S32> ptcl_map_;      ///map from new ptcl_ to original ptcl
    ReallocatableArray<ReallocatableArray<Tptcl>> group_ptcl_;        ///partner group ptcl
    ReallocatableArray<ReallocatableArray<PS::S32>> pert_list_;       ///perturber list
    ReallocatableArray<ReallocatableArray<PS::S32>> group_list_;   /// partner_group member list
    

///    void searchPerturber() {
///        // perturber list
///        n_ptr.resizeNoInitialize(n_ptcl);
///        ptr_list.resizeNoInitialize(n_ptcl);
///        ptr_order.resizeNoInitialize(n_ptcl);
/// 
///        // partner list
///        n_part.resizeNoInitialize(n_ptcl);
///        part_list.resizeNoInitialize(n_ptcl);
/// 
///        // initialization array
///        for(int i=0; i<n_ptcl; i++) {
///            int nbi = ptcl_org[i].n_ngb;
///#ifdef HARD_DEBUG
///            assert(nbi>0);
///            assert(nbi<=n_ptcl);
///            assert(ptr_list[i].capacity()==0);
///            assert(ptr_order[i].capacity()==0);
///            assert(parr_list[i].capacity()==0);
///#endif
///            ptr_list[ i].resizeNoInitialize(nbi);
///            ptr_order[i].resizeNoInitialize(nbi);
///            n_ptr[i] = 0;
/// 
///            part_list[i].resizeNoInitialize(nbi);
///            n_part[i]= 0;
///        }
///        
///        // find perturber and partner
///        for(int i=0; i<n_ptcl; i++) {
///            Tptcl* pi = &ptcl_org[i];
///            int nbi = pi->n_ngb;
///            
///            for(int j=i+1; j<n_ptcl; j++) {
///                PS::F64vec dr = ptcl_org[i].pos-ptcl_org[j].pos;
///                PS::F64 r2 = dr*dr;
///                PS::F64 rout = std::max(ptcl_org[i].r_out,ptcl_org[j].r_out);
///                PS::F64 rout2 = rout*rout+SAFTY_OFFSET_FOR_SEARCH;
///                if (r2<rout2) {
///                    ptr_list [n_ptr[i]] = j;
///                    ptr_order[n_ptr[i]] = n_ptr[j];
/// 
///                    ptr_list [n_ptr[j]] = i;
///                    ptr_order[n_ptr[j]] = n_ptr[i];
/// 
///                    n_ptr[i]++;
///                    n_ptr[j]++;
///#ifdef HARD_DEBUG
///                    assert(n_ptr[i]>nbi);
///                    assert(n_ptr[j]>ptcl_org[j].n_ngb);
///                    assert(n_ptr[i]!=ptr_list[i].size());
///                    assert(n_ptr[j]!=ptr_list[j].size());
///#endif
///                }
///                if (r2<rin2) {
///                    part_list[n_part[i]] = j;
///                    part_list[n_part[j]] = i;
///                    n_part[i]++;
///                    n_part[j]++;
///                }
///            }
///        }
///    }

    void searchPartner(ReallocatableArray<ReallocatableArray<PS::S32>> & ptr_list,
                       Tptcl *ptcl,
                       PS::F64 r_crit2) {
        ptr_list.resizeNoInitialize(n_ptcl);
        
#ifdef HARD_DEBUG
        for(int i=0; i<n_ptcl; i++) {
            assert(ptr_list[i].size()==0);
        }
#endif
        
        // find partner
        for(int i=0; i<n_ptcl; i++) {
            Tptcl* pi = &ptcl[i];
            
            for(int j=i+1; j<n_ptcl; j++) {
                PS::F64vec dr = ptcl[i].pos-ptcl[j].pos;
                PS::F64 r2 = dr*dr;
                if (r2<r_crit2) {
                    ptr_list[i].push_back(j);
                    ptr_list[j].push_back(i);
                }
            }
        }
    }

    void searchPerturber(ReallocatableArray<ReallocatableArray<PS::S32>> & ptr_list,
                         Tptcl *ptcl) {
        // perturber list
        ptr_list.resizeNoInitialize(n_ptcl);
        
#ifdef HARD_DEBUG
        for(int i=0; i<n_ptcl; i++) {
            assert(ptr_list[i].size()==0);
        }
#endif
        
        // find perturber
        for(int i=0; i<n_ptcl; i++) {
            Tptcl* pi = &ptcl[i];
            
            for(int j=i+1; j<n_ptcl; j++) {
                PS::F64vec dr = ptcl[i].pos-ptcl[j].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 rout = std::max(ptcl[i].r_out,ptcl[j].r_out);
                PS::F64 rout2 = rout*rout+SAFTY_OFFSET_FOR_SEARCH;
                if (r2<rout2) {
                    ptr_list[i].push_back(j);
                    ptr_list[j].push_back(i);
                }
            }
        }
    }
    
    void mergeCluster(ReallocatableArray<ReallocatableArray<PS::S32>> & group_list,
                      ReallocatableArray<ReallocatableArray<PS::S32>> & part_list) {
        // partner index with marker
        ReallocatableArray<std::pair<PS::S32,PS::S32>> partner_index; 
        // map index from ptcl_org to partner_index
        ReallocatableArray<PS::S32> reverse_list; 
        reverse_list.resizeNoInitialize(n_ptcl);
        
#ifdef HARD_DEBUG
        assert(group_list.size()==0);
        assert(partner_index.size()==0);
#endif

        for(int i=0; i<n_ptcl; i++) {
            if(part_list[i].size()>0) {
                reverse_list[i] = partner_index.size();
                std:::pair<PS::S32,PS::S32> ipart;
                ipart.first = i;
                ipart.second = -1;
                partner_index.push_back(ipart);
            }
#ifdef HARD_DEBUG
            else {
                reverse_list[i] = -1;
            }
#endif
        }
        PS::S32 n_tot = partner_index.size();

        for(int i=0; i<n_tot; i++) {
            PS::S32 k = partner_index[i].first;
            if(partner_index[i].second>=0) continue;

            PS::S32 npart = connectGroups(i,i,part_list,partner_index,reverse_list);
            group_list[n_groups].resizeNoInitialize(npart);
            PS::S32 inext=partner_index[i].second;
            PS::S32 k=0;
            while (inext!=i) {
                group_list[n_groups][k] = partner_index[inext].first;
                inext=partner_index[inext].second;
                k++;
#ifdef HARD_DEBUG
                assert(k<=npart);
#endif
            }
            n_group.push_back(npart);
            n_groups++;
        }
    }

    PS::S32 connectGroups(const PS::S32 ip,
                          const PS::S32 iend,
                          ReallocatableArray<ReallocatableArray<PS::S32>> & part_list,
                          ReallocatableArray<std::pair<PS::S32,PS::S32>> & partner_index,
                          ReallocatableArray<PS::S32> & reverse_list) {
        PS::S32 n_connected = 0;
        PS::S32 jst,inow;
        PS::S32 kp = partner_index[ip].first;
        for(jst=0; jst<part_list[kp].size(); jst++) {
            PS::S32 inext = reverse_list[part_list[kp][jst]];
            if(partner_index[inext].second<0) {
                partner_index[ip].second = inext;
                inow = inext;
                n_connected++;
                break;
            }
        }
        
        for(int j=jst+1; j<part_list[kp].size(); j++) {
            PS::S32 inext = reverse_list[part_list[kp][j]];
            if(partner_index[inext].second<0) {
                partner_index[inow].second = inext;
                n_connected += connectGroups(inow,inext,part_list,partner_index,reverse_list);
                inow = inext;
                n_connected++;
            }
        }
        return n_connected;
    }

    void generateNewPtcl(ReallocatableArray<Tptcl> ptcl,
                         ReallocatableArray<PS::S32> ptcl_map,
                         Tptcl* ptcl_org,
                         const PS::S32 n_ptcl,
                         ReallocatableArray<ReallocatableArray<Tptcl>> group_ptcl,
                         ReallocatableArray<ReallocatableArray<PS::S32>> group_list) {
#ifdef HARD_DEBUG
        assert(ptcl.size()==0);
        assert(group_ptcl.size()==0);
        assert(ptcl_map.size()==0);
#endif
        PS::S32 n_group_tot = 0;
        PS::S32 masks[n_ptcl]={0};
        for (int i=0; i<group_list.size(); i++) {
#ifdef HARD_DEBUG
            assert(group_ptcl[i].size()==0);
#endif
            PS::S32 n_group = group_list[i].size();
            group_ptcl[i].resizeNoInitialize(n_group);
            for (int j=0; j<n_group; j++) {
                PS::S32 k = group_list[i][j]
                group_ptcl[i][j] = ptcl_org[k];
                n_group_tot++;
                masks[k] = 1;
            }
        }
        PS::S32 n_ptcl_new = n_ptcl-n_group_tot+group_list.size()
#ifdef HARD_DEBUG
        assert(n_ptcl_new>0);
#endif
        ptcl.resizeNoInitialize(n_ptcl_new);
        ptcl_map.resizeNoInitialize(n_ptcl_new);
        for (int i=0; i<group_list.size(); i++) {
            getCenterOfMass(ptcl[i],group_ptcl[i].getPointer(),group_ptcl[i].size());
            ptcl_map[i] = -i;
        }

        PS::S32 ik = group_list.size();
        for (int i=0; i<n_ptcl; i++) {
            if(masks[i]>0) continue;
            ptcl[ik]     = ptcl_org[i];
            ptcl_map[ik] = i;
            ik++;
        }
#ifdef HARD_DEBUG
        assert(n_ptcl_new!=ik);
#endif
    };

    void getCenterOfMass(Tptcl &cm, Tptcl* ptcl, const PS::S32 n_ptcl) {
        PS::F64vec cmr=0;
        PS::F64vec cmv=0;
        PS::F64 cmm = 0;
        cm.r_out = p[0].r_out;
        
        for (int i=0;i<num;i++) {
            cmr += p[i].pos * p[i].mass;
            cmv += p[i].vel * p[i].mass;
            cmm += p[i].mass;
#ifdef HARD_DEBUG
            assert(cm.r_out==p[i].r_out);
#endif
        }
        cmr /= cmm; 
        cmv /= cmm; 
      
        cm.mass = cmm;
        cm.pos = cmr;
        cm.vel = cmv;

    }

//    void UpdatePertList() {
//        ReallocatableArray<PS::S32> markers;
//        markers.ReallocatableArray(n_ptcl);
//        for(int i=0; i<n_groups; i++) {
//            for(int j=1; j<n_ptcl; j++) markers[j]=0;
//            for(int j=1; j<n_group[i]; i++) {
//                PS::S32 inow = group_list[i][j];
//                for(int k=0; k<n_ptr[inow]; k++) {
//                    PS::S32 iptr  = ptr_list[inow][k];
//                    if(markers[iptr]!=0) {
//                        
//                    }
//                }
//            }
//        }
//    }

public:
    void generatelist(Tptcl *ptcl_org, const PS::S32 n_ptcl, const PS::F64 rin){
        ReallocatableArray<ReallocatableArray<PS::S32>> part_list;      ///partner list
        
        searchPartner(part_list, ptcl_org, rin_*rin_);
        mergeCluster(group_list_, part_list);
        generateNewPtcl(ptcl_, ptcl_map_, ptcl_org, n_ptcl, group_ptcl_, group_list_);
        searchPerturber(pert_list_, ptcl_);
    }

    void reverseCopy(Tptcl *ptcl_org, const PS::S32 n_ptcl) {
#ifdef HARD_DEBUG
        int checker[n_ptcl]={0};
#endif
        for (int i=group_list_.size(); i<ptcl_.size(); i++) {
            ptcl_org[ptcl_map_[i]] = ptcl_[i];
#ifdef HARD_DEBUG
            assert(ptcl_map_[i]<n_ptcl);
            checker[ptcl_map_[i]]++;
#endif
        }
        for (int i=0; i<group_list_.size(); i++) 
            for (int j=0; j<group_list_[i].size(); j++) {
#ifdef HARD_DEBUG
                assert(group_list_[i][j]<n_ptcl);
                checker[group_list_[i][j]]++;
#endif
                ptcl_org[group_list_[i][j]] = group_ptcl_[i][j];
            }

#ifdef HARD_DEBUG
        for (int i=0; i<n_ptcl; i++) assert(checker[i]==1);
#endif
    }

};
