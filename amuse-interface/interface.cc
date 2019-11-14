#include "petar.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Interface code
 */

    static PeTar* ptr=NULL;
    static double time_start = 0.0;
 
    int initialize_code() {
        ptr = new PeTar;

        int argc = 0;
        char **argv=NULL;

        //No second MPI init
        //ptr->initialFDPS(argc,argv);
        ptr->initial_fdps_flag = true;

        // default input
        int flag= ptr->readParameters(argc,argv);

        // set writing flag to false
        ptr->input_parameters.write_flag = false;
        ptr->write_flag = false;

        // set restart flat to false
        ptr->file_header.nfile = 0; 
        
        ptr->system_soft.initialize();

        return flag;
    }

    int cleanup_code() {
        delete ptr;
        ptr=NULL;
        return 0;
    }

    int commit_parameters() {
        if (!ptr->read_parameters_flag) return -1;
        return 0;
    }

    int recommit_parameters() {
        // not allown
        return -1;
    }

    int new_particle(int* index_of_the_particle,  double mass, double x, double y, double z, double vx, double vy, double vz, double radius) {
        // if not yet initial the system
        ptr->read_data_flag = true;
        
        FPSoft p;
        long long int id_offset = ptr->hard_manager.ap_manager.id_offset;
        long long int n_glb = ptr->stat.n_real_glb;
        long long int n_loc = ptr->system_soft.getNumberOfParticleLocal();
        if (n_loc!=ptr->stat.n_real_loc) {std::cerr<<"Error: number of particle not match !"; return -1;}
        if (n_glb<id_offset) return -1;
        p.mass = mass;
        p.pos.x = x;
        p.pos.y = y;
        p.pos.z = z;
        p.vel.x = vx;
        p.vel.y = vy;
        p.vel.z = vz;
        p.status.d = 0.0;
        p.mass_bk.d = 0.0;
        if (ptr->initial_parameters_flag) p.calcRSearch(ptr->input_parameters.dt_soft.value);
        p.id = n_glb+1;
        p.adr = n_loc;
        p.rank_org = ptr->my_rank;
        ptr->system_soft.addOneParticle(p);
        ptr->add_particle_index_map(p.adr, p.rank_org, p.id);
        ptr->stat.n_real_loc++;
        ptr->stat.n_real_glb++;
        *index_of_the_particle = p.id;
        return 0;
    }

    int delete_particle(int index_of_the_particle) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if(ptr->my_rank==rank) {
            ptr->system_soft.removeParticle(&index, 1);
            ptr->stat.n_real_loc--;
        }
        ptr->stat.n_real_glb--;
        return 0;
    }

    int get_state(int index_of_the_particle,
                  double * mass, 
                  double * x, double * y, double * z,
                  double * vx, double * vy, double * vz, double * radius){
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (rank<0) return -1;
        if (ptr->my_rank == rank) {
            if (index<0||index>ptr->stat.n_real_loc) return -1;
            FPSoft* p = &(ptr->system_soft[index]);
            *mass = p->mass;
            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
        }
        return 0;
    }

    int set_state(int index_of_the_particle,
                  double mass, 
                  double x, double y, double z,
                  double vx, double vy, double vz, double radius) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->mass = mass;
            p->pos.x = x;
            p->pos.y = y;
            p->pos.z = z;
            p->vel.x = vx;
            p->vel.y = vy;
            p->vel.z = vz;
        }
        return 0;
    }

    int get_mass(int index_of_the_particle, double * mass) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            *mass = p->mass;
        }    
        return 0;
    }

    int set_mass(int index_of_the_particle, double mass) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->mass = mass;
        }    
        return 0;
    }

    int get_radius(int index_of_the_particle, double * radius) {
        // not allown
        return -1;
    }

    int set_radius(int index_of_the_particle, double radius) {
        // not allown
        return -1;
    }

    int set_position(int index_of_the_particle,
                     double x, double y, double z) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->pos.x = x;
            p->pos.y = y;
            p->pos.z = z;
        }    
        return 0;
    }

    int get_position(int index_of_the_particle,
                     double * x, double * y, double * z) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            *x = p->pos.x;
            *y = p->pos.y;
            *z = p->pos.z;
        }    
        return 0;
    }

    int set_velocity(int index_of_the_particle,
                     double vx, double vy, double vz) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->vel.x = vx;
            p->vel.y = vy;
            p->vel.z = vz;
        }    
        return 0;
    }

    int get_velocity(int index_of_the_particle,
                     double * vx, double * vy, double * vz) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            *vx = p->vel.x;
            *vy = p->vel.y;
            *vz = p->vel.z;
        }    
        return 0;
    }

    int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            *ax = p->acc.x;
            *ay = p->acc.y;
            *az = p->acc.z;
        }    
        return 0;
    }

    int set_acceleration(int index_of_the_particle, double ax, double ay, double az) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            p->acc.x = ax; 
            p->acc.y = ay; 
            p->acc.z = az; 
        }    
        return 0;
    }

    int get_potential(int index_of_the_particle, double * potential) {
        if (index_of_the_particle>ptr->stat.n_real_glb||index_of_the_particle<0) return -1;
        int index, rank;
        ptr->get_particle_index_from_id(index, rank, index_of_the_particle);
        if (ptr->my_rank == rank) {
            FPSoft* p = &(ptr->system_soft[index]);
            *potential = p->pot_tot;
        }    
        return 0;
    }

    int evolve_model(double time_next) {
        if (!ptr->initial_step_flag) return -1;
        ptr->input_parameters.time_end.value = time_next*2;
        ptr->evolveToTime(time_next);
        return 0;
    }

    int commit_particles() {
        if (!ptr->read_parameters_flag) return -1;
        if (!ptr->read_data_flag) return -1;
        ptr->initialParameters();
        ptr->initialStep();
        return 0;
    }

    int synchronize_model() {
        return 0;
    }

    int recommit_particles() {
        ptr->initial_step_flag = false;
        return 0;
    }

    int get_eps2(double * epsilon_squared) {
        *epsilon_squared = ptr->input_parameters.eps.value ;
        return 0;
    }

    int set_eps2(double epsilon_squared) {
        ptr->input_parameters.eps.value = epsilon_squared;
        return 0;
    }

    int get_kinetic_energy(double * kinetic_energy) {
        *kinetic_energy = ptr->stat.energy.ekin;
        return 0;
    }

    int get_potential_energy(double * potential_energy) {
        *potential_energy = ptr->stat.energy.epot;
        return 0;
    }

    int get_time(double * time) {
        * time = ptr->stat.time;
        return 0;
    }

    int get_begin_time(double * time) {
        * time = time_start;
        return 0;
    }

    int set_begin_time(double time) {
        time_start = time;
        ptr->stat.time = time_start;
        return 0;
    }

    int get_time_step(double * time_step) {
        * time_step = ptr->input_parameters.dt_soft.value;
        return 0;
    }

    int get_total_mass(double * mass) {
        * mass = ptr->stat.pcm.mass;
        return 0;
    }

    int get_center_of_mass_position(double * x, double * y, double * z) {
        * x = ptr->stat.pcm.pos.x;
        * y = ptr->stat.pcm.pos.y;
        * z = ptr->stat.pcm.pos.z;
        return 0;
    }

    int get_center_of_mass_velocity(double * x, double * y, double * z) {
        * x = ptr->stat.pcm.pos.x;
        * y = ptr->stat.pcm.pos.y;
        * z = ptr->stat.pcm.pos.z;
        return 0;
    }

    int get_total_radius(double * radius) {
        * radius = ptr->stat.half_mass_radius;
        return 0;
    }

    int get_number_of_particles(int * number_of_particles) {
        * number_of_particles = ptr->stat.n_real_glb;
        return 0;
    }

    int get_index_of_first_particle(int * index_of_the_particle) {
        * index_of_the_particle = 1;
        return 0;
    }

    int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle) {
        * index_of_the_next_particle = index_of_the_particle + 1;
        return 0;
    }

#ifdef __cplusplus
}
#endif
