#include <iostream>
#include <fstream>
#include <getopt.h>
#include "galpy_interface.h"
#include "../src/io.hpp"

struct Particle{
    double mass;
    double pos[3];
    double vel[3];
    double acc[3];
    double pot;

    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf",
                          &mass, &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
        if(rcount<7) {
            std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass"
             <<std::setw(_width)<<"pos.x"
             <<std::setw(_width)<<"pos.y"
             <<std::setw(_width)<<"pos.z"
             <<std::setw(_width)<<"vel.x"
             <<std::setw(_width)<<"vel.y"
             <<std::setw(_width)<<"vel.z"
             <<std::setw(_width)<<"acc.x"
             <<std::setw(_width)<<"acc.y"
             <<std::setw(_width)<<"acc.z"
             <<std::setw(_width)<<"pot";
    }


    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20) const{
        _fout<<std::setw(_width)<<mass
             <<std::setw(_width)<<pos[0]
             <<std::setw(_width)<<pos[1]
             <<std::setw(_width)<<pos[2]
             <<std::setw(_width)<<vel[0]
             <<std::setw(_width)<<vel[1]
             <<std::setw(_width)<<vel[2]
             <<std::setw(_width)<<acc[0]
             <<std::setw(_width)<<acc[1]
             <<std::setw(_width)<<acc[2]
             <<std::setw(_width)<<pot;
    }
};

int main(int argc, char** argv){
    
    int arg_label;

    IOParamsGalpy galpy_io;
    opterr = 0;

    // reset optind
    static int long_flag=-1;
    static struct option long_options[] = {
        {0,0,0,0}
    };
    
    int opt_used = 0;
    optind=0;
    int option_index;
    bool help_flag=false;
    bool unit_astro_flag = false;
    bool measure_flag = false;

    while ((arg_label = getopt_long(argc, argv, "-muh", long_options, &option_index)) != -1)
        switch (arg_label) {
        case 0:
            switch (long_flag) {
            default:
                break;
            }
            break;
        case 'm':
            measure_flag = true;
            std::cout<<"Create mesh points to measure the potential evolution in x-y and x-z planes\n";
            opt_used ++;
            break;
        case 'u':
            unit_astro_flag = true;
            std::cout<<"Use the astronomical unit set (Myr, pc, Msun)\n";
            opt_used ++;
            break;
        case 'h':
            std::cout<<"The tool to calculate acceleration and potential for a given particle list \n"
                     <<"Usage: petar.galpy [options] [data file]\n"
                     <<"       data file format: if -m, file contains the mesh parameters\n"
                     <<"                             one line: time, dt_evolve, n_step_evolve, dt_output, x_min, x_max, n_x, y_min, y_max, n_y, z_min, z_max, n_z\n"
                     <<"                         else file contains a particle list\n"
                     <<"                             first line: number of particles, time, cm pos offset (3), cm vel offset (30)\n"
                     <<"                             following lines: mass, pos(3), vel(3)\n"
                     <<"Options:\n"
                     <<"    -m    : instead of reading particle list, generate a mesh of points in the x-y plane and x-z plane to create the acceleration and potential map.\n"
                     <<"            Time-dependent potential is also supported.\n"
                     <<"    -u    : input data use astronomical unit set (Myr, pc, Msun) and set unit scaling factor for Galpy automatically.\n"
                     <<"    -h    : help\n";
            help_flag=true;
            break;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }
    
    galpy_io.print_flag = true;
    galpy_io.read(argc,argv);
    
    if (help_flag) return 0;

    // argc is 1 no input is given
    opt_used++;
    std::string filename;
    if (opt_used<argc) {
        filename = argv[opt_used++];
    }

    FILE* fp;
    if( (fp = fopen(filename.c_str(),"r")) == NULL) {
        fprintf(stderr,"Error: Cannot open file %s.\n", filename.c_str());
        abort();
    }

    GalpyManager galpy_manager;
    if (unit_astro_flag) galpy_io.setStdUnit();

    if (measure_flag) {
        double time, time_out, dt, dt_out, xmin, xmax, ymin, ymax, zmin, zmax;
        int n_step, nx, ny, nz;
        int rcount = fscanf(fp, "%lf %lf %d %lf %lf %lf %d %lf %lf %d %lf %lf %d", 
                            &time, &dt, &n_step, &dt_out, &xmin, &xmax, &nx, &ymin, &ymax, &ny, &zmin, &zmax, &nz);
        if (rcount<12) {
            std::cerr<<"Error: Data reading fails! requiring data number is 12, only obtain "<<rcount<<".\n";
            abort();
        }
        std::ofstream fxy,fxz;
        galpy_manager.initial(galpy_io, time, std::string(), false, true);
        time_out = time;

        for (int i=0; i<=n_step; i++) {
            bool out_flag = (time>=time_out);
            galpy_manager.updatePotential(time, out_flag);

            if (out_flag) {
                fxy.open(("xy"+std::to_string(i)).c_str(), std::ifstream::out);
                fxz.open(("xz"+std::to_string(i)).c_str(), std::ifstream::out);
                fxy<<time<<" "<<nx<<" "<<ny<<std::endl;
                fxz<<time<<" "<<nx<<" "<<nz<<std::endl;
                Particle particle_xy[nx][ny];
                Particle particle_xz[nx][ny];
                for (int j=0; j<nx; j++) {
                    double x = xmin + (xmax-xmin)/(nx-1)*j;
                    for (int k=0; k<ny; k++) {
                        auto& pjk = particle_xy[j][k];
                        pjk.mass = 0;
                        pjk.pos[0] = x;
                        pjk.pos[1] = ymin + (ymax-ymin)/(ny-1)*k;
                        pjk.pos[2] = 0;
                        pjk.vel[0] = pjk.vel[1] = pjk.vel[2] = 0;

                        galpy_manager.calcAccPot(pjk.acc, pjk.pot, time, pjk.pos, pjk.pos); 
                        pjk.printColumn(fxy);
                        fxy<<std::endl;
                    }
                
                    for (int k=0; k<nz; k++) {
                        auto& pjk = particle_xz[j][k];
                        pjk.mass = 0;
                        pjk.pos[0] = x;
                        pjk.pos[1] = 0;
                        pjk.pos[2] = zmin + (zmax-zmin)/(nz-1)*k;
                        pjk.vel[0] = pjk.vel[1] = pjk.vel[2] = 0;

                        galpy_manager.calcAccPot(pjk.acc, pjk.pot, time, pjk.pos, pjk.pos); 
                        pjk.printColumn(fxz);
                        fxz<<std::endl;
                    }
                }
                fxy.close();
                fxz.close();

                time_out += dt_out;
            }
            time += dt;
        }
    }
    else {
        int n = 0;
        double time = 0.0;
        double pos_offset[3], vel_offset[3];
        int rcount = fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf", &n, &time, &pos_offset[0], &pos_offset[1], &pos_offset[2], &vel_offset[0], &vel_offset[1], &vel_offset[2]);
        if(rcount<8) {
            std::cerr<<"Error: Data reading fails! requiring data number is 8, only obtain "<<rcount<<".\n";
            abort();
        }

        assert(n>0);

        galpy_manager.initial(galpy_io, time, std::string(), false, true);

        Particle particles[n];

        Particle::printColumnTitle(std::cout);
        std::cout<<std::endl;

        for (int i=0; i<n; i++) {
            particles[i].readAscii(fp);
            double pos[3] = {particles[i].pos[0] + pos_offset[0],
                             particles[i].pos[1] + pos_offset[1],
                             particles[i].pos[2] + pos_offset[2]};
            galpy_manager.calcAccPot(particles[i].acc, particles[i].pot, time, pos, &particles[i].pos[0]);

            particles[i].printColumn(std::cout);
            std::cout<<std::endl;
        }
    }    
    return 0;
}
