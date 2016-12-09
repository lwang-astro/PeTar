
#include<particle_simulator.hpp>
#include<map>

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

class Wtime{
public:
    static PS::F64 cum;
    static PS::F64 interval;
    static PS::F64 take_snp;
    static PS::F64 soft;
    static PS::F64 soft_force;
    static PS::F64 hard;
    static PS::F64 hard_1body;
    static PS::F64 hard_2body;
    static PS::F64 hard_multi;
    static PS::F64 hard_gather_data;
    static PS::F64 hard_scatter_data;
    static PS::F64 hard_copy_h2s;
    static PS::F64 hard_merge;

    static PS::F64 cum_offset;
    static PS::F64 interval_offset;
    static PS::F64 take_snp_offset;
    static PS::F64 soft_offset;
    static PS::F64 soft_force_offset;
    static PS::F64 hard_offset;
    static PS::F64 hard_1body_offset;
    static PS::F64 hard_2body_offset;
    static PS::F64 hard_multi_offset;
    static PS::F64 hard_gather_data_offset;
    static PS::F64 hard_scatter_data_offset;
    static PS::F64 hard_copy_h2s_offset;
    static PS::F64 hard_merge_offset;

    static PS::F64 hard_2body_select_system;
    static PS::F64 hard_2body_select_system_offset;

    static void dump(std::ofstream & fout, const PS::F64 time_sys, const PS::S64 n_loop){
	fout<<"time_sys= "<<time_sys
	    <<" wtime_cum= "<<cum
	    <<" wtime_interval= "<<interval<<" "<<interval / n_loop
	    <<" wtime_hard= "<<hard<<" "<<hard / n_loop
	    <<" wtime_soft= "<<soft<<" "<<soft / n_loop
	    <<" wtime_soft_force= "<<soft_force<<" "<<soft_force / n_loop
	    <<" wtime_hard_1body= "<<hard_1body<<" "<<hard_1body / n_loop
	    <<" wtime_hard_2body_select_system= "<<hard_2body_select_system<<" "<<hard_2body_select_system / n_loop
	    <<" wtime_hard_2body= "<<hard_2body<<" "<<hard_2body / n_loop
	    <<" wtime_hard_multi= "<<hard_multi<<" "<<hard_multi / n_loop
	    <<" wtime_hard_gather_data= "<<hard_gather_data<<" "<<hard_gather_data / n_loop
	    <<" wtime_copy_h2s= "<<hard_copy_h2s<<" "<<hard_copy_h2s / n_loop
	    <<" wtime_hard_merge= "<<hard_merge<<" "<<hard_merge / n_loop
	    <<" wtime_hard_scatter_data= "<<hard_scatter_data<<" "<<hard_scatter_data/n_loop
	    <<" take_snp= "<<take_snp<<std::endl;
    }

    static void clear(){
	cum = interval = take_snp 
	    = soft = soft_force = hard = hard_1body = hard_2body = hard_multi 
	    = hard_gather_data = hard_scatter_data = hard_copy_h2s = hard_merge 
	    = hard_2body_select_system = 0.0;
    }

};
PS::F64 Wtime::cum;
PS::F64 Wtime::interval;
PS::F64 Wtime::take_snp;
PS::F64 Wtime::soft;
PS::F64 Wtime::soft_force;
PS::F64 Wtime::hard;
PS::F64 Wtime::hard_1body;
PS::F64 Wtime::hard_2body;
PS::F64 Wtime::hard_multi;
PS::F64 Wtime::hard_gather_data;
PS::F64 Wtime::hard_scatter_data;
PS::F64 Wtime::hard_copy_h2s;
PS::F64 Wtime::hard_merge;

PS::F64 Wtime::cum_offset;
PS::F64 Wtime::interval_offset;
PS::F64 Wtime::take_snp_offset;
PS::F64 Wtime::soft_offset;
PS::F64 Wtime::soft_force_offset;
PS::F64 Wtime::hard_offset;
PS::F64 Wtime::hard_1body_offset;
PS::F64 Wtime::hard_2body_offset;
PS::F64 Wtime::hard_multi_offset;
PS::F64 Wtime::hard_gather_data_offset;
PS::F64 Wtime::hard_scatter_data_offset;
PS::F64 Wtime::hard_copy_h2s_offset;
PS::F64 Wtime::hard_merge_offset;


PS::F64 Wtime::hard_2body_select_system;
PS::F64 Wtime::hard_2body_select_system_offset;

//! for event counts (L.Wang)
class Counts{
public:
  std::map<PS::S32,PS::S32> Ncluster; ///<Histogram of number of particles in clusters

  void Ncluster_count(const PS::S32 n) {
    if (Ncluster.count(n)) Ncluster[n]++;
    else Ncluster[n]=1;
  }
};
