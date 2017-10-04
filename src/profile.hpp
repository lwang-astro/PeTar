#include<particle_simulator.hpp>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<map>

/*
template<class Tdinfo, class Tsystem, class Ttree>
class Profile{

private:
    Tdinfo * dinfo_;
    Tsystem * system_;
    Ttree * tree_; 
    PS::F64 n_op_ep_ep_;
    PS::F64 n_op_ep_sp_;
    PS::F64 flops_per_core_;
public:

    static void dumpTimeProfile(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
        PS::S32 id = 0;
        fout<<"collect_sample_particle= "<<tp.collect_sample_particle<<", max= "<<tp_max.collect_sample_particle<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"decompose_domain= "<<tp.decompose_domain<<", max= "<<tp_max.decompose_domain<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"exchange_particle= "<<tp.exchange_particle<<", max= "<<tp_max.exchange_particle<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree<<", max= "<<tp_max.set_particle_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"set_root_cell= "<<tp.set_root_cell<<", max= "<<tp_max.set_root_cell<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_local_tree= "<<tp.make_local_tree<<", max= "<<tp_max.make_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    morton_sort_local_tree= "<<tp.morton_sort_local_tree<<", max= "<<tp_max.morton_sort_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    link_cell_local_tree= "<<tp.link_cell_local_tree<<", max= "<<tp_max.link_cell_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree<<", max= "<<tp_max.calc_moment_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_LET_1st= "<<tp.make_LET_1st<<", max= "<<tp_max.make_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st<<", max= "<<tp_max.exchange_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_LET_2nd= "<<tp.make_LET_2nd<<", max= "<<tp_max.make_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd<<", max= "<<tp_max.exchange_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree<<", max= "<<tp_max.set_particle_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    morton_sort_global_tree= "<<tp.morton_sort_global_tree<<", max= "<<tp_max.morton_sort_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    link_cell_global_tree= "<<tp.link_cell_global_tree<<", max= "<<tp_max.link_cell_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_global_tree= "<<tp.make_global_tree<<", max= "<<tp_max.make_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree<<", max= "<<tp_max.calc_moment_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_force = "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<std::endl;
    } 

    static void dumpTimeProfile0(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
        PS::S32 id = 0;
        fout<<"collect_sample_particle= "<<tp.collect_sample_particle<<", max= "<<tp_max.collect_sample_particle<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"decompose_domain= "<<tp.decompose_domain<<", max= "<<tp_max.decompose_domain<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<std::endl;
    } 

    static void dumpTimeProfile1(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
        PS::S32 id = 2;
        fout<<"exchange_particle= "<<tp.exchange_particle<<", max= "<<tp_max.exchange_particle<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<std::endl;
    } 

    static void dumpTimeProfile2(const PS::TimeProfile & tp, const PS::TimeProfile & tp_max, const PS::S32 rank_max[], std::ostream & fout){
        PS::S32 id = 3;
        fout<<"set_particle_local_tree= "<<tp.set_particle_local_tree<<", max= "<<tp_max.set_particle_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"set_root_cell= "<<tp.set_root_cell<<", max= "<<tp_max.set_root_cell<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_local_tree= "<<tp.make_local_tree<<", max= "<<tp_max.make_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    morton_sort_local_tree= "<<tp.morton_sort_local_tree<<", max= "<<tp_max.morton_sort_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    link_cell_local_tree= "<<tp.link_cell_local_tree<<", max= "<<tp_max.link_cell_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_moment_local_tree= "<<tp.calc_moment_local_tree<<", max= "<<tp_max.calc_moment_local_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_LET_1st= "<<tp.make_LET_1st<<", max= "<<tp_max.make_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"exchange_LET_1st= "<<tp.exchange_LET_1st<<", max= "<<tp_max.exchange_LET_1st<<", rank= "<<rank_max[id++]<<std::endl;

        fout<<"    exchange_LET_1st__a2a_n= "<<tp.exchange_LET_1st__a2a_n<<", max= "<<tp_max.exchange_LET_1st__a2a_n<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    exchange_LET_1st__a2a_ep= "<<tp.exchange_LET_1st__a2a_ep<<", max= "<<tp_max.exchange_LET_1st__a2a_ep<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    exchange_LET_1st__a2a_sp= "<<tp.exchange_LET_1st__a2a_sp<<", max= "<<tp_max.exchange_LET_1st__a2a_sp<<", rank= "<<rank_max[id++]<<std::endl;

        fout<<"make_LET_2nd= "<<tp.make_LET_2nd<<", max= "<<tp_max.make_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"exchange_LET_2nd= "<<tp.exchange_LET_2nd<<", max= "<<tp_max.exchange_LET_2nd<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"set_particle_global_tree= "<<tp.set_particle_global_tree<<", max= "<<tp_max.set_particle_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"make_global_tree= "<<tp.make_global_tree<<", max= "<<tp_max.make_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    morton_sort_global_tree= "<<tp.morton_sort_global_tree<<", max= "<<tp_max.morton_sort_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"    link_cell_global_tree= "<<tp.link_cell_global_tree<<", max= "<<tp_max.link_cell_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_moment_global_tree= "<<tp.calc_moment_global_tree<<", max= "<<tp_max.calc_moment_global_tree<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<"calc_force = "<<tp.calc_force<<", max= "<<tp_max.calc_force<<", rank= "<<rank_max[id++]<<std::endl;
        fout<<std::endl;
    }

    static void getTimeProfileMax(const PS::TimeProfile & tp, const PS::S32 rank, PS::TimeProfile & tp_max, PS::S32 rank_max[]){
        PS::S32 id = 0;
        PS::Comm::getMaxValue(tp.collect_sample_particle, rank, tp_max.collect_sample_particle, rank_max[id++]);
        PS::Comm::getMaxValue(tp.decompose_domain, rank, tp_max.decompose_domain, rank_max[id++]);

        PS::Comm::getMaxValue(tp.exchange_particle, rank, tp_max.exchange_particle, rank_max[id++]);

        PS::Comm::getMaxValue(tp.set_particle_local_tree, rank, tp_max.set_particle_local_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.set_root_cell, rank, tp_max.set_root_cell, rank_max[id++]);
        PS::Comm::getMaxValue(tp.make_local_tree, rank, tp_max.make_local_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.morton_sort_local_tree, rank, tp_max.morton_sort_local_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.link_cell_local_tree, rank, tp_max.link_cell_local_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.calc_moment_local_tree, rank, tp_max.calc_moment_local_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.make_LET_1st, rank, tp_max.make_LET_1st, rank_max[id++]);
        PS::Comm::getMaxValue(tp.exchange_LET_1st, rank, tp_max.exchange_LET_1st, rank_max[id++]);

        PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_n, rank, tp_max.exchange_LET_1st__a2a_n, rank_max[id++]);
        PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_ep, rank, tp_max.exchange_LET_1st__a2a_ep, rank_max[id++]);
        PS::Comm::getMaxValue(tp.exchange_LET_1st__a2a_sp, rank, tp_max.exchange_LET_1st__a2a_sp, rank_max[id++]);

        PS::Comm::getMaxValue(tp.make_LET_2nd, rank, tp_max.make_LET_2nd, rank_max[id++]);
        PS::Comm::getMaxValue(tp.exchange_LET_2nd, rank, tp_max.exchange_LET_2nd, rank_max[id++]);
        PS::Comm::getMaxValue(tp.set_particle_global_tree, rank, tp_max.set_particle_global_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.make_global_tree, rank, tp_max.make_global_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.morton_sort_global_tree, rank, tp_max.morton_sort_global_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.link_cell_global_tree, rank, tp_max.link_cell_global_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.calc_moment_global_tree, rank, tp_max.calc_moment_global_tree, rank_max[id++]);
        PS::Comm::getMaxValue(tp.calc_force, rank, tp_max.calc_force, rank_max[id++]);
    }
    
    Profile(Tdinfo  * _dinfo, 
            Tsystem  * _system, 
            Ttree * _tree, 
            const PS::F64 _n_op_ep_ep, 
            const PS::F64 _n_op_ep_sp, 
            const PS::F64 _flops_per_core){
        dinfo_ = _dinfo;
        system_ = _system;
        tree_ = _tree;
        n_op_ep_ep_ = _n_op_ep_ep;
        n_op_ep_sp_ = _n_op_ep_sp;
        flops_per_core_ = _flops_per_core;
    }

    void dump(std::ofstream & fout, const PS::F64 time_sys, const PS::S64 n_loop, const PS::F64 wtime_tot){
        //static PS::S64 n_loop_old = 0;
        //const PS::S64 dn_loop = n_loop - n_loop_old;
        const PS::S64 dn_loop = n_loop;
        PS::TimeProfile tp_dinfo = dinfo_->getTimeProfile();
        PS::TimeProfile tp_system = system_->getTimeProfile();
        PS::TimeProfile tp_tree = tree_->getTimeProfile();
        PS::TimeProfile tp_dinfo_max, tp_system_max, tp_tree_max;
        const PS::S32 n_profile_max = 100; 
        PS::S32 rank_dinfo_max[n_profile_max], rank_system_max[n_profile_max], rank_tree_max[n_profile_max];
        getTimeProfileMax(tp_dinfo, PS::Comm::getRank(), tp_dinfo_max, rank_dinfo_max);
        getTimeProfileMax(tp_system, PS::Comm::getRank(), tp_system_max, rank_system_max);
        getTimeProfileMax(tp_tree, PS::Comm::getRank(), tp_tree_max, rank_tree_max);
        PS::CountT n_int_ep_ep = tree_->getNumberOfInteractionEPEPGlobal();
        PS::CountT n_int_ep_sp = tree_->getNumberOfInteractionEPSPGlobal();
        PS::CountT n_op_tot = n_int_ep_ep * n_op_ep_ep_ + n_int_ep_sp * n_op_ep_sp_;
        fout<<"soft part break down"<<std::endl;
        fout<<"time_sys= "<<time_sys<<" n_loop= "<<n_loop<<std::endl;
	fout<<"n_loc= "<<system_->getNumberOfParticleLocal()<<" n_glb= "<<system_->getNumberOfParticleGlobal()<<std::endl;
        fout<<"speed= "<<(PS::F64)(n_op_tot)/(wtime_tot)*1e-12<<"[Tflops]"<<std::endl;
        fout<<"PS::Comm::getNumberOfThread()= "<<PS::Comm::getNumberOfThread()<<std::endl;
        fout<<"efficiency= "<<(PS::F64)(n_op_tot)/(wtime_tot)/(flops_per_core_*PS::Comm::getNumberOfProc()*PS::Comm::getNumberOfThread())<<std::endl;
        fout<<"wtime_tot= "<<wtime_tot<<std::endl;
        fout<<"n_op_tot= "<<n_op_tot<<std::endl;
        //timer.dump(fout);
        dumpTimeProfile0(tp_dinfo, tp_dinfo_max, rank_dinfo_max, fout);
        dumpTimeProfile1(tp_system, tp_system_max, rank_system_max, fout);
        fout<<"n_int_ep_ep= "<<n_int_ep_ep<<" n_int_ep_sp= "<<n_int_ep_sp<<std::endl;
        fout<<"ni_ave= "<<(PS::F64)(system_->getNumberOfParticleGlobal() * dn_loop) / tree_->getNumberOfWalkGlobal()
            <<" nj_ave(EP-EP)= "<<(PS::F64)(n_int_ep_ep) / (system_->getNumberOfParticleGlobal() * dn_loop)
            <<" nj_ave(EP-SP)= "<<(PS::F64)(n_int_ep_sp) / (system_->getNumberOfParticleGlobal() * dn_loop)<<std::endl;
        dumpTimeProfile2(tp_tree, tp_tree_max, rank_tree_max, fout);
        //n_loop_old = n_loop;
    }

    void clear(){
        dinfo_->clearTimeProfile();
        system_->clearTimeProfile();
	tree_->clearCounterAll();
        //tree_->clearTimeProfile();
        //tree_->clearNumberOfInteraction();
    }

};
*/

struct Tprofile{
    PS::F64 time;
    const char* name;
    
    Tprofile(const char* name_): time(0.0), name(name_) {}
    
    void start(){
        time -= PS::GetWtime();
    }

    void end(){
        time += PS::GetWtime();
    }
    
    void print(std::ostream & fout, const PS::S32 divider=1){
        fout<<name<<": "<<time/divider<<std::endl;
    }

    void dump(std::ostream & fout, const PS::S32 width=20, const PS::S32 divider=1){
        fout<<std::setw(width)<<time/divider;
    }

    void reset() {
        time = 0.0;
    }

};

struct NumCounter{
    PS::S64 n_;
    const char* name_;
    
    NumCounter(const char* name): n_(0), name_(name) {}

    NumCounter &operator++() {
        (this->n_)++;
        return *this;
    }

    NumCounter &operator+=(const PS::S64 n) {
        this->n_ += n;
        return *this;
    }

    NumCounter &operator=(const PS::S64 n) {
        this->n_ = n;
        return *this;
    }

    void print(std::ostream & fout, const PS::S32 divider=1){
        fout<<name_<<": "<<((divider==1)?n_:(PS::F64)n_/divider)<<std::endl;
    }

    void dump(std::ostream & fout, const PS::S32 width=20, const PS::S32 divider=1){
        fout<<std::setw(width)<<((divider==1)?n_:(PS::F64)n_/divider);
    }

};

class SysProfile{
public:
	Tprofile tot;		   
	Tprofile hard_tot;	   
	Tprofile hard_single;	   
	Tprofile hard_isolated;
	Tprofile hard_connected;
	Tprofile soft_tot;	   
	Tprofile search_cluster;
    const PS::S32 n_profile;
    
    SysProfile(): tot           (Tprofile("Total         ")),
                  hard_tot      (Tprofile("Hard total    ")),
                  hard_single   (Tprofile("Hard single   ")),
                  hard_isolated (Tprofile("Hard isolated ")),
                  hard_connected(Tprofile("Hard connected")),
                  soft_tot      (Tprofile("Soft total    ")),
                  search_cluster(Tprofile("Search cluster")),
                  n_profile(7) {}

	void print(std::ostream & fout, const PS::F64 time_sys, const PS::S64 n_loop=1){
        fout<<"Time: "<<time_sys<<std::endl;
        
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->print(fout, n_loop);
        }
    }

    void dump(std::ofstream & fout, const PS::S32 width=20, const PS::S64 n_loop=1){
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->dump(fout, width, n_loop);
        }
    }

    void clear(){
        for(PS::S32 i=0; i<n_profile; i++) {
            Tprofile* iptr = (Tprofile*)this+i;
            iptr->reset();
        }
    }

};

class SysCounts{
public:
    NumCounter hard_single;
    NumCounter hard_isolated;
    NumCounter hard_connected;
    NumCounter cluster_isolated;
    NumCounter cluster_connected;
    NumCounter ARC_substep_sum;
    const PS::S32 n_counter;
    std::map<PS::S32,PS::S32> n_cluster; ///<Histogram of number of particles in clusters

    SysCounts(): hard_single     (NumCounter("Hard single   ")),
                 hard_isolated   (NumCounter("Hard isolated ")),
                 hard_connected  (NumCounter("Hard connected")),
                 cluster_isolated (NumCounter("Cluster isolated ")),
                 cluster_connected(NumCounter("Cluster connected")),
                 ARC_substep_sum  (NumCounter("ARC sub-steps sum")),
                 n_counter(6) {}

    void cluster_count(const PS::S32 n, const PS::S32 ntimes=1) {
        if (n_cluster.count(n)) n_cluster[n] += ntimes;
        else n_cluster[n]=ntimes;
    }

    void print(std::ostream & fout, const PS::S32 width=20, const PS::S64 n_loop=1) {
        for(PS::S32 i=0; i<n_counter; i++) {
            NumCounter* iptr = (NumCounter*)this+i;
            iptr->print(fout, n_loop);
        }
        fout<<"Number of members in clusters:\n";
        for(auto i=n_cluster.begin(); i!=n_cluster.end(); ++i) fout<<std::setw(width)<<i->first;
        fout<<std::endl;
        for(auto i=n_cluster.begin(); i!=n_cluster.end(); ++i) fout<<std::setw(width)<<i->second/((n_loop==1)?1:(PS::F64)n_loop);
        fout<<std::endl;
    }

    void dump(std::ofstream & fout, const PS::S32 width=20, const PS::S64 n_loop=1){
        for(PS::S32 i=0; i<n_counter; i++) {
            NumCounter* iptr = (NumCounter*)this+i;
            iptr->dump(fout, width, n_loop);
        }
        for(auto i=n_cluster.begin(); i!=n_cluster.end(); ++i)
            fout<<std::setw(width)<<i->first<<std::setw(width)<<i->second/((n_loop==1)?1:(PS::F64)n_loop);
    }

    void clear() {
        for(PS::S32 i=0; i<n_counter; i++) {
            NumCounter* iptr = (NumCounter*)this+i;
            *iptr = 0;
        }
        n_cluster.clear();
    }
};
