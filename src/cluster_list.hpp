#pragma once
#include<particle_simulator.hpp>
#include<unordered_map>
#include<map>
#include"ptcl.hpp"
#include"kepler.hpp"

//extern const PS::F64 SAFTY_OFFSET_FOR_SEARCH;

#define ID_PHASE_SHIFT 4
#define ID_PHASE_MASKER 0xF

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
    PtclComm(const Tp &p): Ptcl(p) {}
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
                                           const PS::F64 r_out,  // 1/(r_out-r_in)
                                           const PS::F64 r_in,
                                           const PS::F64ort pos_domain[],
                                           const PS::F64 eps_sq=0.0){
        static PS::ReallocatableArray<PtclOuter> * ptcl_outer = NULL;
        static PS::ReallocatableArray< std::pair<PS::S32, PS::S32> > * id_ngb_multi_cluster = NULL;
        const PS::S32 n_thread = PS::Comm::getNumberOfThread();
        if(ptcl_outer==NULL) ptcl_outer = new PS::ReallocatableArray<PtclOuter>[n_thread];
        if(id_ngb_multi_cluster==NULL) id_ngb_multi_cluster = new PS::ReallocatableArray< std::pair<PS::S32, PS::S32> >[n_thread];
        const PS::S32 my_rank = PS::Comm::getRank();
        //        const PS::S32 n_proc_tot = PS::Comm::getNumberOfProc();
        const PS::S32 n_loc = sys.getNumberOfParticleLocal();
        const PS::F64 r_oi_inv = 1.0/(r_out-r_in);
        const PS::F64 r_A = (r_out-r_in)/(r_out+r_in);
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
                // self-potential correction 
                if (sys[i].status==0) sys[i].pot_tot += sys[i].mass/r_out;
                else if (sys[i].status<0) sys[i].pot_tot += sys[i].mass_bk/r_out;

                if(sys[i].n_ngb == 0){
                    // no neighbor
                    adr_sys_one_cluster_[ith].push_back(i);
                }
                else{
                    // has neighbor
                    ptcl_cluster_[ith].push_back( PtclCluster(sys[i].id, i, id_ngb_multi_cluster[ith].size(), sys[i].n_ngb, false, NULL, my_rank) );
#ifdef HARD_DEBUG
                    if(sys[i].id<0&&sys[i].status<0) {
                        std::cerr<<"Error: ghost particle exist! i="<<i<<std::endl;
                        abort();
                    }
//                    if(sys[i].id==10477) {
//                        sys[i].print(std::cerr);
//                        std::cerr<<std::endl;
//                    }
                    assert(i<n_loc);
#endif
                    //PS::S32 n_tmp2 = 0;
                    //PS::S32 n_tmp3 = 0;
                    //PS::S32 n_tmp4 = 0;
                    for(PS::S32 ii=0; ii<sys[i].n_ngb+1; ii++){
                        if( (nbl+ii)->id == sys[i].id ){
                            //n_tmp2++;
                            continue;
                        }
                        id_ngb_multi_cluster[ith].push_back( std::pair<PS::S32, PS::S32>(sys[i].id, (nbl+ii)->id) );
                        PS::S32 pot_control_flag;
                        Tepj* nj = nbl+ii;
                        if (nj->id>0) {
                            if (nj->status==0) pot_control_flag = 0; // single
                            else if (nj->status<0) pot_control_flag = 1; // member
                            else pot_control_flag = 2; // fake

                            CalcAccPotShortWithLinearCutoff
                                (sys[i].pos,    sys[i].acc,      sys[i].pot_tot,
                                 (nbl+ii)->pos, (nbl+ii)->mass,  (nbl+ii)->mass_bk,  pot_control_flag, eps_sq,
                                 r_oi_inv,    r_A,   r_out,  r_in);
                        }
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
                    ptcl_send_.push_back(PtclComm(p));
                    ptcl_send_.back().id_cluster = cluster_loc[i].id_;
#ifdef HARD_DEBUG
                    if(ptcl_send_.back().id<0&&ptcl_send_.back().status<0) {
                        std::cerr<<"Error! sending particle is ghost! adr="<<adr_sys<<std::endl;
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
                              const PS::ReallocatableArray<Tphard> & ptcl_hard,
                              PS::ReallocatableArray<PS::S32> &removelist){
        //  const PS::S32 my_rank = PS::Comm::getRank();
        //  const PS::S32 n_proc  = PS::Comm::getNumberOfProc();
        const PS::S32 n = ptcl_hard.size();
        for(PS::S32 i=0; i<n; i++){
            const PS::S32 adr = ptcl_hard[i].adr_org;
            if( adr >= 0){
                //assert( sys[adr].id == ptcl_hard[i].id);
                sys[adr].DataCopy(ptcl_hard[i]);
                if(sys[adr].id<0&&sys[adr].status<0) removelist.push_back(adr);
            }
            else{
                //assert( ptcl_recv_[-(adr+1)].id == ptcl_hard[i].id );
                ptcl_recv_[-(adr+1)].DataCopy(ptcl_hard[i]);
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
            //assert(sys[adr].id == ptcl_send_[i].id);
            sys[adr].DataCopy(ptcl_send_[i]);
            if(sys[adr].id<0&&sys[adr].status<0) removelist.push_back(adr);
        }
        // remove empty particles cannot do here
        // sys.removeParticle(removelist.getPointer(), removelist.size());
    }
#endif


    const PS::ReallocatableArray<PS::S32> & getAdrSysOneCluster(){
        return adr_sys_one_cluster_[0];
    }


};


template<class Tptcl>
class SearchGroup{
private:
    typedef std::pair<PS::S32, PS::S32> PLinker;
//    typedef std::pair<PS::S32, PS::F64> RCList;

    PS::ReallocatableArray<PS::S32> p_list_;
    // group member list
    PS::ReallocatableArray<PS::S32> group_list_;
    PS::ReallocatableArray<PS::S32> group_list_disp_;
    PS::ReallocatableArray<PS::S32> group_list_n_;
    //PS::ReallocatableArray<PS::S32> group_list_pert_adr_; // soft pert adr in soft_pert_list_
    // fake particle list
    PS::ReallocatableArray<PS::S32> soft_pert_list_;
    // PS::ReallocatableArray<PS::S32> soft_pert_list_disp_;

//    PS::ReallocatableArray<Tptcl> ptcl_;               ///new ptcl with c.m.
//    PS::ReallocatableArray<PS::S32> ptcl_map_;      ///map from new ptcl_ to original ptcl
//    PS::ReallocatableArray<PS::ReallocatableArray<Tptcl>> group_ptcl_;        ///partner group ptcl
//    PS::ReallocatableArray<RCList> Rout_change_list_;  /// Rout change list

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
///            assert(part_list[i].capacity()==0);
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

    void searchPartner(PS::ReallocatableArray<PS::S32> & part_list,
                       PS::ReallocatableArray<PS::S32> & part_list_disp,
                       PS::ReallocatableArray<PS::S32> & part_list_n,
                       PS::ReallocatableArray<PS::S32> & p_list,
                       Tptcl *ptcl,
//                       const PS::S32 n_ptcl,
                       const PS::F64 r_crit2) {
        PS::S32 n = p_list.size();
        part_list.clearSize();
        part_list_disp.reserve(n);
        part_list_disp.resizeNoInitialize(n);
        part_list_n.reserve(n);
        part_list_n.resizeNoInitialize(n);
        
        // find partner
        PS::S32 offset = 0;
        for(int i=0; i<n; i++) {
            PS::S32 ip = p_list[i];
            ptcl[ip].mass_bk = 0.0;
            part_list_n[i] = 0;
            part_list_disp[i] = offset;
            for(int j=0; j<n; j++) {
                PS::S32 jp = p_list[j];
                if(ip==jp) continue;
                PS::F64vec dr = ptcl[ip].pos-ptcl[jp].pos;
                PS::F64 r2 = dr*dr;
                PS::F64 mass_factor=std::max(ptcl[ip].mass,ptcl[jp].mass)/std::min(ptcl[ip].mass,ptcl[jp].mass);
                //PS::F64 mass_factor=std::max(ptcl[ip].mass,ptcl[jp].mass)*Ptcl::mean_mass_inv;
                if (r2<r_crit2*mass_factor) {
                    part_list.push_back(j);
                    part_list_n[i]++;
                    offset++;
                }
                else {
                    //store perturbation force
                    ptcl[ip].mass_bk += ptcl[jp].mass/r2; 
                }
            }
        }
    }

    void mergeCluster(PS::ReallocatableArray<PS::S32> & group_list,
                      PS::ReallocatableArray<PS::S32> & group_list_disp,
                      PS::ReallocatableArray<PS::S32> & group_list_n,
                      PS::ReallocatableArray<PS::S32> & p_list,
                      PS::S32 part_list[],
                      PS::S32 part_list_disp[],
                      PS::S32 part_list_n[]) {

        const PS::S32 n_ptcl = p_list.size();
        // partner index with marker
        PS::ReallocatableArray<PLinker> partner_index; 
        partner_index.reserve(n_ptcl);

        // map index from ptcl_org to partner_index
        PS::ReallocatableArray<PS::S32> reverse_list; 
        reverse_list.reserve(n_ptcl);
        reverse_list.resizeNoInitialize(n_ptcl);

#ifdef HARD_DEBUG
        assert(group_list.size()==0);
        assert(group_list_disp.size()==0);
        assert(group_list_n.size()==0);
#endif

        for(int i=0; i<n_ptcl; i++) {
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
            group_list.push_back(p_list[partner_index[i].first]);
            PS::S32 inext=partner_index[i].second;
            PS::S32 k=1;
            while (inext!=i) {
                group_list.push_back(p_list[partner_index[inext].first]);
                inext=partner_index[inext].second;
                k++;
#ifdef HARD_DEBUG
                assert(k<=npart);
#endif
            }
#ifdef HARD_DEBUG
            assert(k==npart);
            if(npart>n_ptcl) {
                std::cerr<<"Error: connect group particle number mismatch: npart ="<<npart<<" ; n_ptcl = "<<n_ptcl<<std::endl;
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


    template<class Tptree>
    PS::S32 setGroupMemberPars(Tptree &bin, const PS::S32 n_split, const PS::S64 id_offset, const PS::S32 bid, const bool is_top) {
        PS::S32 nloc = 0;
        for(int i=0; i<2; i++) {
            if(bin.member[i]->status!=0) {
                if(is_top) nloc += setGroupMemberPars(*(Tptree*)bin.member[i], n_split, id_offset, bin.member[i]->id, false);
                else nloc += setGroupMemberPars(*(Tptree*)bin.member[i], n_split, id_offset, bid, false);
            }
            else {
                if(is_top) bin.member[i]->status = -std::abs(id_offset+bin.member[i]->id*n_split);
                else bin.member[i]->status = -std::abs(id_offset+bid*n_split);
                bin.member[i]->mass_bk = bin.member[i]->mass;
                bin.member[i]->r_search = bin.r_search;
#ifdef SPLIT_MASS
                bin.member[i]->mass    = 0.0;
#endif
                nloc += 1;
            }
        }
        return nloc;
    }

    //! generate kepler sampling artificial particles
    /* also store the binary parameters in mass_bk of first 4 pairs
       acc, ecc
       peri, tstep (integrator step estimation ),
       inc, OMG,
       omg, ecca
       
     */
    template<class Tptree>
    void keplerOrbitGenerator(Tptcl* ptcl_org,
                              PS::ReallocatableArray<Tptcl> & ptcl_extra,
                              PS::ReallocatableArray<PS::S32> & empty_list,
                              Tptree &bin,
                              const PS::S64 id_offset,
                              const PS::S32 n_split) {
#ifdef TIDAL_TENSOR
        const PS::F64 dE = 8.0*atan(1.0)/(n_split-4);
#else
        const PS::F64 dE = 8.0*atan(1.0)/n_split;
#endif
        //const PS::F64 dt = bin.peri/n_split;
        if (n_split<8) {
            std::cerr<<"N_SPLIT to small to save binary parameters, should be >= 8!";
            abort();
        }
        PS::F64 bindata[5][2];
        //Tptcl* plist[n_split][2];
        /*
          acc, ecc
          peri, tstep,
          inc, OMG,
          omg, ecca
          tperi, stable_factor
         */
        bindata[0][0] = bin.ax;
        bindata[0][1] = bin.ecc;
        bindata[1][0] = bin.peri;
        bindata[1][1] = bin.tstep;
        bindata[2][0] = bin.inc; 
        bindata[2][1] = bin.OMG; 
        bindata[3][0] = bin.omg; 
        bindata[3][1] = bin.ecca;
        bindata[4][0] = bin.tperi;
        bindata[4][1] = bin.stable_factor;

#ifdef SPLIT_MASS        
        PS::F64* pm[n_split][2];
        PS::F64 mfactor;
        PS::F64 mnormal=0.0;
#endif

        // Make sure the ptcl_extra will not make new array due to none enough capacity during the following loop, otherwise the p[j] pointer will point to wrong position
        ptcl_extra.reserveEmptyAreaAtLeast(2*n_split-empty_list.size());
        // First 4 is used for tidal tensor points
        // remaining is used for sample points
        for (int i=0; i<n_split; i++) {
            Tptcl* p[2];
            for(int j=0; j<2; j++) {
                if(empty_list.size()>0) {
                    PS::S32 k = empty_list.back();
                    empty_list.decreaseSize(1);
                    p[j] = &ptcl_org[k];
                }
                else {
                    ptcl_extra.push_back(Tptcl());
                    p[j] = &ptcl_extra.back();
                }
                p[j]->mass = bin.member[j]->mass;
                p[j]->id = id_offset + (bin.member[j]->id)*n_split + i;
                //p[j]->r_search = bin.member[j]->r_search;
                p[j]->r_search = bin.r_search;
                p[j]->status = (bin.id<<ID_PHASE_SHIFT)|i;
#ifdef SPLIT_MASS
#ifdef TIDAL_TENSOR
                if(i>=4) 
#endif
                    pm[i][j] = &(p[j]->mass);
#endif
            }
#ifdef TIDAL_TENSOR
            if (i>=4) {
                PS::S32 iph = i-4;
#else 
                PS::S32 iph = i;
#endif
                // center_of_mass_shift(*(Tptcl*)&bin,p,2);
                // generate particles at different orbitial phase
                OrbParam2PosVel(p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass,
                                bin.ax, bin.ecc, bin.inc, bin.OMG, bin.omg, dE*iph);
                //DriveKeplerOrbParam(p[0]->pos, p[1]->pos, p[0]->vel, p[1]->vel, p[0]->mass, p[1]->mass, (i+1)*dt, bin.ax, bin.ecc, bin.inc, bin.OMG, bin.omg, bin.peri, bin.ecca);
#ifdef TIDAL_TENSOR
            }
            else {
                /* Assume apo-center distance is the maximum length inside box
                   Then the lscale=apo/(2*sqrt(2))
                 */
                
                PS::F64 lscale = bin.ax*(1+bin.ecc)*0.35;
                // 8 points box 
                switch(i) {
                case 0:
                    p[0]->pos = PS::F64vec(lscale, 0,      -lscale) + bin.pos;
                    p[1]->pos = PS::F64vec(0,      lscale, -lscale) + bin.pos;
                    break;
                case 1:
                    p[0]->pos = PS::F64vec(-lscale, 0,      -lscale) + bin.pos;
                    p[1]->pos = PS::F64vec(0,      -lscale, -lscale) + bin.pos;
                    break;
                case 2:
                    p[0]->pos = PS::F64vec(lscale, 0,      lscale) + bin.pos;
                    p[1]->pos = PS::F64vec(0,      lscale, lscale) + bin.pos;
                    break;
                case 3:
                    p[0]->pos = PS::F64vec(-lscale, 0,      lscale) + bin.pos;
                    p[1]->pos = PS::F64vec(0,      -lscale, lscale) + bin.pos;
                    break;
                default:
                    std::cerr<<"Error: index >= 4!\n";
                    abort();
                }
                p[0]->vel = bin.vel;
                p[1]->vel = bin.vel;
                p[0]->mass = p[1]->mass = 0.0;
            }
#endif
            // binary parameters
            if(i<5) {
                p[0]->mass_bk = bindata[i][0];
                p[1]->mass_bk = bindata[i][1];
            }
            else if(i==5) {
                p[0]->mass_bk = p[0]->mass;
                p[1]->mass_bk = p[1]->mass;
            }
            else if(i==6) {
                p[0]->mass_bk = 0.0; // indicate the order
                p[1]->mass_bk = 1.0;
            }
#ifdef TIDAL_TENSOR
            if(i>=4) {
#endif
#ifdef SPLIT_MASS
                PS::F64vec dvvec= p[0]->vel - p[1]->vel;
                PS::F64 odv = 1.0/std::sqrt(dvvec*dvvec);
                for(int j=0; j<2; j++) p[j]->mass *= odv;
                //if(i==0) p[j]->mass /= 2.0*n_split;
                //else p[j]->mass /= n_split;
                mnormal += odv;
#else
                for(int j=0; j<2; j++) p[j]->mass = 0;
#endif
                center_of_mass_correction(*(Tptcl*)&bin,p,2);
#ifdef HARD_DEBUG
                //check rsearch consistence:
                PS::F64 rsearch_bin = bin.r_search+bin.ax*(1+bin.ecc);
                for(int j=0; j<2; j++) {
                    PS::F64vec dp = p[j]->pos-bin.pos;
                    PS::F64 dr = dp*dp;
                    assert(dr<=rsearch_bin*rsearch_bin);
//                    if(p[j]->id==10477) {
//                        std::cerr<<"i="<<i<<" dr="<<sqrt(dr)<<std::endl;
//                        p[j]->print(std::cerr);
//                        std::cerr<<std::endl;
//                    }
                }
#endif
#ifdef TIDAL_TENSOR
            }
#endif
        }

#ifdef SPLIT_MASS
        mfactor = 1.0/mnormal;
#ifdef TIDAL_TENSOR
        for (int i=4; i<n_split; i++) 
#else
        for (int i=0; i<n_split; i++) 
#endif
            for (int j=0; j<2; j++) 
                *pm[i][j] *= mfactor;
#endif
        //mfactor *= 0.5;
        //*pm[0][0] *= mfactor;
        //*pm[0][1] *= mfactor;
        //mfactor /= n_split;
        PS::S32 nbin = setGroupMemberPars(bin, n_split, id_offset, bin.id, true);
        Tptcl* pcm;
        if(empty_list.size()>0) {
            pcm = &ptcl_org[empty_list.back()];
            empty_list.decreaseSize(1);
        }
        else {
            ptcl_extra.push_back(Tptcl());
            pcm = &ptcl_extra.back();
        }
        pcm->mass_bk = bin.mass;
        pcm->mass = 0.0;
        pcm->pos = bin.pos;
        pcm->vel = bin.vel;
        pcm->id  = - std::abs(bin.id);
        pcm->r_search = bin.r_search;
#ifdef TIDAL_TENSOR
        pcm->r_search += bin.ax*(1+bin.ecc);  // depend on the mass ratio, the upper limit distance to c.m. from all members and artifical particles is apo-center distance
#endif
        pcm->status = nbin;

    }

    template<class Tptree>
    void generateNewPtcl(Tptcl* ptcl_org,
//                         const PS::S32 n_ptcl,
                         PS::ReallocatableArray<PS::S32> & p_list,
                         PS::ReallocatableArray<Tptcl> & ptcl_extra,
                         PS::ReallocatableArray<PS::S32> & group_list,
                         PS::ReallocatableArray<PS::S32> & group_list_disp,
                         PS::ReallocatableArray<PS::S32> & group_list_n,
                         PS::ReallocatableArray<PS::S32> & empty_list,
                         const PS::F64 rbin,
                         const PS::F64 rin,
                         const PS::F64 rout,
                         const PS::F64 dt_tree,
                         const PS::S64 id_offset,
                         const PS::S32 n_split){
#ifdef HARD_DEBUG
//        assert(ptcl.size()==0);
//        assert(group_ptcl.size()==0);
//        assert(n_ptcl<=n_ptcl_org);
//        assert(ptcl_map.size()==0);
#endif
        for (int i=0; i<p_list.size(); i++) ptcl_org[p_list[i]].status = 0;

        const PS::S32 n_groups = group_list_n.size();
        for (int i=0; i<n_groups; i++) {
            PS::ReallocatableArray<Tptree> bins;
            PS::ReallocatableArray<Tptree*> stab_bins;
            bins.reserve(n_groups);
            stab_bins.reserve(n_groups);
            bins.resizeNoInitialize(group_list_n[i]-1);
            keplerTreeGenerator(bins.getPointer(), &group_list[group_list_disp[i]], group_list_n[i], ptcl_org, dt_tree);
            PS::S32 fstab = stabilityCheck<Tptcl>(stab_bins,bins.back(),rbin,rin,rout);
            if (fstab) {
                keplerOrbitGenerator(ptcl_org, ptcl_extra, empty_list, bins.back(), id_offset, n_split);
            }
            else {
                for (int i=0; i<stab_bins.size(); i++) {
                    keplerOrbitGenerator(ptcl_org, ptcl_extra, empty_list, *stab_bins[i], id_offset, n_split);
                }
            }
        }

        for (int i=0; i<empty_list.size(); i++) {
            PS::S32 ik = empty_list[i];
            ptcl_org[ik].mass = 0.0;
            ptcl_org[ik].id = -1;
            ptcl_org[ik].status = -1;
        }

//#ifdef HARD_DEBUG
//        PS::S32 n_ptcl_new = n_ptcl-n_group_tot+group_list.size();
//        if(n_ptcl_new<=0) {
//            std::cerr<<"Error: data size unmatch! n_ptcl = "<<n_ptcl<<" n_group_tot = "<<n_group_tot<<" group_list.size = "<<group_list.size()<<std::endl;
//            abort();
//        }
//#endif
        // ptcl.resizeNoInitialize(n_ptcl_new);
        // ptcl_map.resizeNoInitialize(n_ptcl_new);
        //PS::S32 nstart = 0;
        //PS::S32 nnext  = group_list.size();
        //for (int i=0; i<group_list.size(); i++) {
        //    if(masks[i]==0) {
        //        while(masks[nnext]==0) {
        //            nnext++;
        //            assert(nnext>=n_ptcl);
        //        }
        //        ptcl_org[nnext] = ptcl_org[i];
        //        ptcl_org[nnext].id_group = -1;
        //        nnext++;
        //    }
        //}
        //PS::S32 nend = n_ptcl-1;
        //for (int i=nnext, i<n_ptcl_new; i++) {
        //    if(masks[nnext]>0) {
        //        while(masks[nend]>0) nend--;
        //        if(nnext<nend) {
        //            ptcl_org[nnext] = ptcl_org[nend];
        //            ptcl_org[nnext].id_group = -1;
        //        }
        //    }
        //    else ptcl_org[i].id_group = -1;
        //}
//#ifdef HARD_DEBUG
//        for (int i=group_list.size(); i<n_ptcl_new; i++) {
//            assert(masks[i]==0);
//            assert(ptcl_org[i].id_group==-1);
//        }
//#endif
        //nnext = n_ptcl_new;
        //for (int i=0; i<group_list.size(); i++) {
        //    if(masks[i]>=0) {
        //    
        //        if (i<n_groups_glb) {
        //            gnow = &group_ptcl_glb[id_group_glb[i]];
        //        ptcl_org[i].id_group = id_group_glb[i];
        //    }
        //    else {
        //        gnow = &group_ptcl_loc[i-n_groups_glb];
        //    }
        //    getCenterOfMass(ptcl_org[i],gnow->getPointer(),gnow->size());
        //    // ptcl_map[i] = -i;
        //    for (int j=0; j<group_ptcl[i].size(); j++) {
        //        ptcl_org[nnext] = ptcl_org[i];
        //        ptcl_org[nnext].mass = 0.0;
        //        nnext++;
        //    }
        //}
        //assert(nnext==n_ptcl);

        // PS::S32 ik = group_list.size();
        // for (int i=0; i<n_ptcl; i++) {
        //     if(masks[i]>0) continue;
        //     ptcl[ik]     = ptcl_org[i];
        //     ptcl_map[ik] = i;
        //     ik++;
        // }
        //#ifdef HARD_DEBUG
        //        assert(n_ptcl_new!=ik);
        //#endif
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

//    void UpdatePertList() {
//        PS::ReallocatableArray<PS::S32> markers;
//        markers.PS::ReallocatableArray(n_ptcl);
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

    //! Find groups
    /* Algorithm:
       1. Find c.m. and create cm_adr (map<c.m.id, p_list.index>)
       2. Find artificial members and create: 
               fake_cm   (map<fake.id, positive c.m.id>)
               fake_order(map<fake.id, order in soft_pert_list_ (0,1)>)
               fake_adr  (map<fake.id, position index in soft_pert_list_)
               soft_pert_list_ (fake order depend on which is found first)
       3. Find single and group members
               group member order depend on which is found first
               group member status is updated to fake_order
       4. correct mass of c.m. and members
     */
    void findGroups(Tptcl* ptcl_org,
                    const PS::S32 n_ptcl,
                    const PS::S32 n_split) {
        p_list_.reserve(n_ptcl);
        p_list_.clearSize();

        group_list_disp_.reserve(n_ptcl);
        group_list_disp_.clearSize();
        group_list_n_.reserve(n_ptcl);
        group_list_n_.clearSize();

        std::map<PS::S32,PS::S32> cm_adr;
        std::map<PS::S32,PS::S32> fake_cm;
        std::map<PS::S32,PS::S32> fake_adr;
        std::map<PS::S32,PS::S32> fake_order;

#ifdef HARD_DEBUG
        PS::ReallocatableArray<PS::S32> icount;
        icount.reserve(n_ptcl);
        for (int i=0; i<n_ptcl; i++) icount.pushBackNoCheck(0);
#endif

        // find cm
        PS::S32 n_members = 0;
        PS::S32 n_unused = 0;
        for (int i=0; i<n_ptcl; i++) {
            if(ptcl_org[i].id<0) {
                if(ptcl_org[i].status>0) {
                    cm_adr[ptcl_org[i].id] = p_list_.size();
                    p_list_.pushBackNoCheck(i);
                    group_list_disp_.pushBackNoCheck(n_members);
                    group_list_n_.pushBackNoCheck(0);
                    n_members += ptcl_org[i].status;
#ifdef HARD_DEBUG
                    icount[i]++;
#endif
                }
                else n_unused++;
            }
        }

        PS::S32 n_groups = p_list_.size();  // group number

        group_list_.reserve(n_members);
        group_list_.resizeNoInitialize(n_members);
        //group_list_pert_adr_.reserve(n_members);
        //group_list_pert_adr_.resizeNoInitialize(n_members);

        //std::map<PS::S32,std::pair<PS::S32,PS::S32>> mem_order;
        //gloffset.reserve(n_members);

        PS::S32 soft_size = n_groups*2*n_split+n_unused;
        PS::S32 i_unused = n_groups*2*n_split;
        soft_pert_list_.reserve(soft_size);
        soft_pert_list_.resizeNoInitialize(soft_size);

        PS::ReallocatableArray<PS::S32> soft_pert_list_disp;
        soft_pert_list_disp.reserve(n_groups);
        soft_pert_list_disp.resizeNoInitialize(n_groups);

        PS::ReallocatableArray<PS::S32> soft_pert_list_n;
        soft_pert_list_n.reserve(soft_size);
        soft_pert_list_.resizeNoInitialize(soft_size);

        for (int i=0; i<n_groups; i++) {
            soft_pert_list_disp[i] = i*n_split*2;
            soft_pert_list_n[i] = 0;
        }

#ifdef HARD_DEBUG
        PS::ReallocatableArray<PS::S32> scount;
        scount.reserve(soft_size);
        for (int i=0; i<soft_size; i++) scount.pushBackNoCheck(0);
        PS::ReallocatableArray<PS::S32> gcount;
        gcount.reserve(n_members);
        for (int i=0; i<n_members; i++) gcount.pushBackNoCheck(0);
#endif

        for (int i=0; i<n_ptcl; i++) {
            if(ptcl_org[i].id>0&&ptcl_org[i].status>0) { // fake members
                PS::S32 idcm = (ptcl_org[i].status>>ID_PHASE_SHIFT);
                PS::S32 phase = ptcl_org[i].status&ID_PHASE_MASKER;
                PS::S32 icm = cm_adr.at(-idcm);
                PS::S32 id0 = ptcl_org[i].id-phase;
                auto itr = fake_cm.find(id0);
                if(itr==fake_cm.end()) {
                    fake_cm[id0] = idcm;
                    fake_order[id0] = soft_pert_list_n[icm];
                    PS::S32 ipos = soft_pert_list_disp[icm] + soft_pert_list_n[icm]++;
                    fake_adr[id0] = ipos;
                    soft_pert_list_[ipos+phase*2] = i;
#ifdef HARD_DEBUG
                    icount[i]++;
                    scount[ipos+phase*2]++;
#endif
                }
                else {
                    soft_pert_list_[fake_adr[id0]+phase*2] = i;
#ifdef HARD_DEBUG
                    icount[i]++;
                    scount[fake_adr[id0]+phase*2]++;
#endif
                }
                //std::pair<PS::S32,PS::S32>* iadr = &mem_order[ptcl_org[i].id];
                // PS::S32 icm = iadr->first;
                // PS::S32 iorder = iadr->second;
                // PS::S32 iphase = ptcl_org[i].status;
            }
        }

        for (int i=0; i<n_ptcl; i++) {
            if(ptcl_org[i].id>=0) {
                if(ptcl_org[i].status==0) {
                    p_list_.pushBackNoCheck(i); //single
#ifdef HARD_DEBUG
                    icount[i]++;
#endif
                }
                else if(ptcl_org[i].status<0) {  // members
                    PS::S32 id_fake = -ptcl_org[i].status;
                    PS::S32 iadr_cm   = cm_adr.at(-fake_cm.at(id_fake));
                    PS::S32 ilst    = group_list_disp_[iadr_cm]+group_list_n_[iadr_cm]++;
                    group_list_[ilst] = i;
                    //group_list_pert_adr_[ilst] = fake_adr.at(id_fake);

                    // give status the order of fake member order saved in soft_pert_list
                    ptcl_org[i].status = fake_order.at(id_fake); 
                    
                    //mem_order[ptcl_org[i].id] = std::pair<PS::S32, PS::S32>(iadr_cm, group_list_n_[iadr_cm]);  // record cm index and member id order
#ifdef HARD_DEBUG
                    assert(ilst<group_list_.size());
                    assert(group_list_n_[iadr_cm]<=ptcl_org[p_list_[iadr_cm]].status);
                    icount[i]++;
                    gcount[group_list_disp_[iadr_cm]+group_list_n_[iadr_cm]-1]++;
#endif
                }
            }
            else if(ptcl_org[i].status<0){
                soft_pert_list_[i_unused++] = i; // unused 
#ifdef HARD_DEBUG
                icount[i]++;
                scount[i_unused-1]++;
#endif
            }
        }


#ifdef HARD_DEBUG
        assert(i_unused-n_groups*2*n_split==n_unused);
        if(p_list_.size()+n_members+soft_size!=n_ptcl) {
            std::cerr<<"Rank: "<<PS::Comm::getRank()<<" n_members="<<n_members<<" p_list_.size="<<p_list_.size()<<" soft_size="<<soft_size<<" n_ptcl="<<n_ptcl<<" n_groups="<<n_groups<<std::endl;
        }
        for(int i=0; i<icount.size(); i++) assert(icount[i]==1);
        for(int i=0; i<gcount.size(); i++) assert(gcount[i]==1);
        for(int i=0; i<scount.size(); i++) assert(scount[i]==1);
        assert(p_list_.size()+n_members+soft_size==n_ptcl);
#endif
        
        for (int i=0; i<n_groups; i++) { // get correct mass
            //PS::F64 masscm = 0.0;
            for (int j=0; j<group_list_n_[i]; j++) {
                PS::S32 i_mem  = group_list_[group_list_disp_[i]+j];
                //PS::S32 i_fake = soft_pert_list_[soft_pert_list_disp_[i]+j];
                ptcl_org[i_mem].mass = ptcl_org[i_mem].mass_bk;
                //masscm += ptcl_org[i_mem].mass;
            }
            const PS::S32 pid = p_list_[i];
            ptcl_org[pid].mass = ptcl_org[pid].mass_bk;
        }
        
        // get 
//#ifdef HARD_DEBUG
//        assert(group_list.size()==0);
//        assert(adr_cm.size()==0);
//#endif
//        n_groups = 0;
//        id_group.resizeNoInitialize(n_ptcl);
//        id_group_map.resizeNoInitialize(n_ptcl);
//        for (int i=0; i<n_ptcl; i++) {
//            id_group[i] = ptcl_org[i].id_group;
//            if(id_group[i]>=0) {
//                id_group_map[i] = n_groups;
//                group_list.push_back(id_group[i]);
//                adr_cm.push_back(i);
//                n_groups++;
//            }
//            else {
//                id_group_map[i] = -1;
//            }
//        }
    }

    //void searchPerturber(Tptcl *ptcl_org, const PS::S32 n_ptcl) {
    //    searchPerturber(pert_list_, ptcl_org, n_ptcl);
    //}

    void searchAndMerge(Tptcl *ptcl_org, const PS::F64 rcrit){
        PS::ReallocatableArray<PS::S32> part_list;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_disp;      ///partner list
        PS::ReallocatableArray<PS::S32> part_list_n;      ///partner list
        
        searchPartner(part_list, part_list_disp, part_list_n, p_list_, ptcl_org, rcrit*rcrit);
        mergeCluster(group_list_, group_list_disp_, group_list_n_, p_list_, part_list.getPointer(), part_list_disp.getPointer(), part_list_n.getPointer());
        // #ifdef HARD_DEBUG
        // assert(Rout_change_list_.size()==0);
        // #endif
        // checkRoutChange(Rout_change_list_, group_list, ptcl_org);

    }

    /* rmax: maximum binary radius
       dt_tree: tree time step
    */
    void generateList(Tptcl *ptcl_org, 
                      PS::ReallocatableArray<Tptcl> & ptcl_extra,
                      const PS::F64 rbin,
                      const PS::F64 rin,
                      const PS::F64 rout,
                      const PS::F64 dt_tree,
                      const PS::S64 id_offset,
                      const PS::S32 n_split) {
        if (n_split>(1<<ID_PHASE_SHIFT)) {
            std::cerr<<"Error! ID_PHASE_SHIFT is too small for phase split! shift bit: "<<ID_PHASE_SHIFT<<" n_split: "<<n_split<<std::endl;
            abort();
        }
        generateNewPtcl<PtclTree<Tptcl>>(ptcl_org, p_list_, ptcl_extra, group_list_, group_list_disp_, group_list_n_, soft_pert_list_, rbin, rin, rout, dt_tree, id_offset,n_split);
        //searchPerturber(pert_list_, ptcl_org, n_ptcl);
    }

    void resolveGroups() {
        const PS::S32 ng = group_list_n_.size();
        for (int i=0; i<ng; i++) {
            soft_pert_list_.push_back(p_list_[i]);
            p_list_[i] = group_list_[i];
        }
        PS::S32 n = group_list_.size();
#ifdef HARD_DEBUG
        assert(n>=0);
        assert(n+p_list_.size()-ng<=p_list_.capacity());
#endif
        for (int i=ng; i<n; i++) {
#ifdef HARD_DEBUG
            assert(p_list_.size()<p_list_.capacity());
#endif
            p_list_.pushBackNoCheck(group_list_[i]);
        }
        group_list_.clearSize();
        group_list_n_.clearSize();
        group_list_disp_.clearSize();
    }

//    Tptcl* getPtcl() {
//        return ptcl_.getPointer();
//    };

    //PS::ReallocatableArray<PS::S32> *getPerts() {
    //    return pert_list_.getPointer();
    //}

    const PS::S32 getPtclN() const {
        return p_list_.size();
    }

    PS::S32* getPtclList() const {
        return p_list_.getPointer();
    }

    const PS::S32 getPtclIndex(const std::size_t i) const {
        return p_list_.getPointer()[i];
    }

//    PS::ReallocatableArray<Tptcl> *getPGroups() {
//        return group_ptcl_.getPointer();
//    }

    const PS::S32 getNumOfGroups() const {
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

    PS::S32 getGroupPertAdr(const std::size_t igroup, const std::size_t imember, const std::size_t iphase, const PS::S32 n_split) const {
        return soft_pert_list_[igroup*2*n_split+imember+2*iphase];
    }
    // 
    PS::S32* getGroupPertList(const std::size_t igroup, const PS::S32 n_split) {
        return &soft_pert_list_[igroup*2*n_split];
    }

    // assume the binary information stored in artificial star mass_bk
    void getBinPars(Binary &bin, const Tptcl ptcl_org[], const std::size_t igroup, const PS::S32 n_split) {
        const PS::S32 ioff = igroup*2*n_split;
        if (ptcl_org[soft_pert_list_[ioff+12]].mass_bk==0.0) {
            bin.ax   = ptcl_org[soft_pert_list_[ioff  ]].mass_bk;
            bin.ecc  = ptcl_org[soft_pert_list_[ioff+1]].mass_bk;
            bin.peri = ptcl_org[soft_pert_list_[ioff+2]].mass_bk;
            bin.tstep= ptcl_org[soft_pert_list_[ioff+3]].mass_bk;
            bin.inc  = ptcl_org[soft_pert_list_[ioff+4]].mass_bk;
            bin.OMG  = ptcl_org[soft_pert_list_[ioff+5]].mass_bk;
            bin.omg  = ptcl_org[soft_pert_list_[ioff+6]].mass_bk;
            bin.ecca = ptcl_org[soft_pert_list_[ioff+7]].mass_bk;
            bin.tperi= ptcl_org[soft_pert_list_[ioff+8]].mass_bk;
            bin.stable_factor= ptcl_org[soft_pert_list_[ioff+9]].mass_bk;
            bin.m1   = ptcl_org[soft_pert_list_[ioff+10]].mass_bk;
            bin.m2   = ptcl_org[soft_pert_list_[ioff+11]].mass_bk;
        }
        else {
            bin.ax   = ptcl_org[soft_pert_list_[ioff+1]].mass_bk;
            bin.ecc  = ptcl_org[soft_pert_list_[ioff  ]].mass_bk;
            bin.peri = ptcl_org[soft_pert_list_[ioff+3]].mass_bk;
            bin.tstep= ptcl_org[soft_pert_list_[ioff+2]].mass_bk;
            bin.inc  = ptcl_org[soft_pert_list_[ioff+5]].mass_bk;
            bin.OMG  = ptcl_org[soft_pert_list_[ioff+4]].mass_bk;
            bin.omg  = ptcl_org[soft_pert_list_[ioff+7]].mass_bk;
            bin.ecca = ptcl_org[soft_pert_list_[ioff+6]].mass_bk;
            bin.tperi= ptcl_org[soft_pert_list_[ioff+9]].mass_bk;
            bin.stable_factor= ptcl_org[soft_pert_list_[ioff+8]].mass_bk;
            bin.m1   = ptcl_org[soft_pert_list_[ioff+11]].mass_bk;
            bin.m2   = ptcl_org[soft_pert_list_[ioff+10]].mass_bk;
        }
    }
    
//    RCList* getRoutChangeList() {
//        return Rout_change_list_.getPointer();
//    }

//    void reverseCopy(Tptcl *ptcl_org, const PS::S32 n_ptcl) {
//#ifdef HARD_DEBUG
//        int checker[n_ptcl]={0};
//#endif
//        for (int i=group_list_.size(); i<ptcl_.size(); i++) {
//            ptcl_org[ptcl_map_[i]] = ptcl_[i];
//#ifdef HARD_DEBUG
//            assert(ptcl_map_[i]<n_ptcl);
//            checker[ptcl_map_[i]]++;
//#endif
//        }
//        for (int i=0; i<group_list_.size(); i++) 
//            for (int j=0; j<group_list_[i].size(); j++) {
//#ifdef HARD_DEBUG
//                assert(group_list_[i][j]<n_ptcl);
//                checker[group_list_[i][j]]++;
//#endif
//                ptcl_org[group_list_[i][j]] = group_ptcl_[i][j];
//            }
// 
//#ifdef HARD_DEBUG
//        for (int i=0; i<n_ptcl; i++) assert(checker[i]==1);
//#endif
//    }

};
