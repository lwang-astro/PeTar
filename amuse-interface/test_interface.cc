#include "interface.h"
#include <cstdio>
#include <cassert>
#include "mpi.h"

int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);
    initialize_code();
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    commit_parameters();

    int index[4];
    new_particle(&index[0], 1, 2, 3, 4, 5, 6, 7, 8);
    new_particle(&index[1], 11, 12, 13, 14, 15, 16, 17, 18);
    new_particle(&index[2], 21, 22, 23, 24, 25, 26, 27, 28);

    commit_particles();

    MPI_Bcast(index, 4, MPI_INT, 0, MPI_COMM_WORLD);

    double m,x,y,z,vx,vy,vz,r;
    int error = get_state(index[0],&m,&x,&y,&z,&vx,&vy,&vz,&r);
    if (my_rank==0) {
        if (error<0) printf("get state error\n");
        assert(m==1);
        assert(x==2);
        assert(y==3);
        assert(z==4);
        assert(vx==5);
        assert(vy==6);
        assert(vz==7);
    }
    error = get_state(index[2],&m,&x,&y,&z,&vx,&vy,&vz,&r);
    if (my_rank==0) {
        if (error<0) printf("get state error\n");
        assert(m==21);
        assert(x==22);
        assert(y==23);
        assert(z==24);
        assert(vx==25);
        assert(vy==26);
        assert(vz==27);
    }
    //printf("I%d m:%f x:%f y:%f z:%f vx:%f vy:%f vz:%f r:%f\n", index[0],m,x,y,z,vx,vy,vz,r);
    
    error = get_mass(index[1],&m);
    if (my_rank==0) {
        if (error<0) printf("get mass error\n");
        assert(m==11);
    }

    error = get_position(index[1],&x, &y, &z);
    if (my_rank==0) {
        if (error<0) printf("get position error\n");
        assert(x==12);
        assert(y==13);
        assert(z==14);
    }
    error = get_velocity(index[1], &vx, &vy, &vz);
    if (my_rank==0) {
        if (error<0) printf("get velocity error\n");
        assert(vx==15);
        assert(vy==16);
        assert(vz==17);
    }

    m=41;
    x=42;
    y=43;
    z=44;
    vx=45;
    vy=46;
    vz=47;
    error = set_state(index[1], m, x, y, z, vx, vy, vz, r);
    if (my_rank==0) {
        if (error<0) printf("set state error\n");
    }
    error = get_state(index[1],&m,&x,&y,&z,&vx,&vy,&vz,&r);
    if (my_rank==0) {
        if (error<0) printf("get state error\n");
        assert(m==41);
        assert(x==42);
        assert(y==43);
        assert(z==44);
        assert(vx==45);
        assert(vy==46);
        assert(vz==47);
    }

    m=51;
    x=52;
    y=53;
    z=54;
    vx=55;
    vy=56;
    vz=57;
    error = set_mass(index[1], m);
    if (my_rank==0) {
        if (error<0) printf("set mass error\n");
    }
    error = get_mass(index[1],&m);
    if (my_rank==0) {
        if (error<0) printf("get mass error\n");
        assert(m==51);
    }

    error = set_position(index[1],x, y, z);
    if (my_rank==0) {
        if (error<0) printf("set position error\n");
    }
    error = get_position(index[1],&x, &y, &z);
    if (my_rank==0) {
        if (error<0) printf("get position error\n");
        assert(x==52);
        assert(y==53);
        assert(z==54);
    }

    error = set_velocity(index[1], vx, vy, vz);
    if (my_rank==0) {
        if (error<0) printf("set velocity error\n");
    }
    error = get_velocity(index[1], &vx, &vy, &vz);
    if (my_rank==0) {
        if (error<0) printf("get velocity error\n");
        assert(vx==55);
        assert(vy==56);
        assert(vz==57);
    }
    recommit_particles();

    for (int k=0; k<10; k++) {
        error = evolve_model(0.25*k);
        if (my_rank==0) {
            if (error<0) printf("evolve model error\n");
        }
    }
    
    double time, ekin, epot;
    get_kinetic_energy(&ekin);
    get_potential_energy(&epot);
    get_time(&time);
    if (my_rank==0) {
        printf("T=%f Ekin=%f Epot=%f \n",time,ekin,epot);
    }

    cleanup_code();

    return 0;
}
