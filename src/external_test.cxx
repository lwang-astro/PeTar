#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cassert>
#ifdef GALPY
#include "galpy_interface.h"
#elif AGAMA
#include "agama_interface.h"
#endif
#include "io.hpp"

struct Particle{
    double mass;
    double pos[3];
    double vel[3];
    double acc[3];
    double pot;
    double den;

    void readAscii(FILE* fp) {
        int rcount=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf",
                          &mass, &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]);
        if(rcount<7) {
            std::cerr<<"Error: Data reading fails! requiring data number is 7, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void writeBinary(FILE* fp) const {
        fwrite(this, sizeof(*this), 1, fp);
    }

    void writeBinary(std::ostream & _fout) const {
        _fout.write(reinterpret_cast<const char*>(this), sizeof(*this));
    }

    void readBinary(std::istream & _fin) {
        _fin.read(reinterpret_cast<char*>(this), sizeof(*this));
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
             <<std::setw(_width)<<"pot"
             <<std::setw(_width)<<"den";
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
             <<std::setw(_width)<<pot
             <<std::setw(_width)<<den;
    }
};

int main(int argc, char** argv){
    
    int arg_label;

#ifdef GALPY    
    IOParamsGalpy galpy_parameters;
#elif AGAMA
    IOParamsAgama agama_parameters;
#endif
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
    //bool unit_astro_flag = false;
    bool measure_flag = false;
    bool out_binary = true;

    while ((arg_label = getopt_long(argc, argv, "-muAh", long_options, &option_index)) != -1)
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
        //case 'u':
        //    unit_astro_flag = true;
        //    std::cout<<"Use the astronomical unit set (Myr, pc, Msun)\n";
        //    opt_used ++;
        //    break;
        case 'A':
            out_binary = false;
            std::cout<<"Output in ASCII format\n";
            opt_used ++;
            break;
        case 'h':
            std::cout<<"The tool to calculate acceleration, potential and mass density for a given particle list \n"
                     <<"Usage: petar.external [options] [data file]\n"
                     <<"       data file format: if -m, file contains the mesh parameters\n"
                     <<"                             one line: time, dt_evolve, n_step_evolve, dt_output, x_min, x_max, n_x, y_min, y_max, n_y, z_min, z_max, n_z\n"
                     <<"                         else file contains a particle list\n"
                     <<"                             first line: number of particles, time, cm pos offset (3), cm vel offset (30)\n"
                     <<"                             following lines: mass, pos(3), vel(3)\n"
                     <<"Options:\n"
                     <<"    -m    : instead of reading particle list, generate a mesh of points in the x-y plane and x-z plane to create the acceleration, potential and density map.\n"
                     <<"            output filenames are xy[time] and xz[time]\n"
                     <<"            header line:  time nx ny\n"
                     <<"            each line: mass x y z vx vy vz ax ay az pot den\n"
                     <<"            Time-dependent potential is also supported.\n"
                     <<"    -A    : output in ASCII format when -m mode is used (default: BINARY)\n"
                     <<"    -h    : help\n";
            help_flag=true;
            break;
        case '?':
            opt_used +=2;
            break;
        default:
            break;
        }

#ifdef GALPY        
    galpy_parameters.print_flag = true;
    galpy_parameters.read(argc,argv);
#elif AGAMA
    agama_parameters.print_flag = true;
    agama_parameters.read(argc,argv);
#endif
    
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

#ifdef GALPY    
    GalpyManager galpy_manager;
#elif AGAMA
    AgamaManager agama_manager;
#endif
    //if (unit_astro_flag) input_parameters.setStdUnit();

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

#ifdef GALPY        
        galpy_manager.initial(galpy_parameters, time, std::string(), false, true);
        int nset = galpy_manager.getNSet();
#elif AGAMA
        agama_manager.initial(agama_parameters, time, true);
#endif
        time_out = time;

        for (int i=0; i<=n_step; i++) {
            bool out_flag = (time>=time_out);
#ifdef GALPY            
            galpy_manager.updatePotential(time, out_flag);
#endif

            if (out_flag) {
                if (out_binary) {
                    fxy.open(("xy"+std::to_string(i)).c_str(), std::ifstream::binary | std::ifstream::out);
                    fxz.open(("xz"+std::to_string(i)).c_str(), std::ifstream::binary | std::ifstream::out);
                    fxy.write(reinterpret_cast<const char*>(&time), sizeof(time));
                    fxy.write(reinterpret_cast<const char*>(&nx), sizeof(nx));
                    fxy.write(reinterpret_cast<const char*>(&ny), sizeof(ny));
                    fxz.write(reinterpret_cast<const char*>(&time), sizeof(time));
                    fxz.write(reinterpret_cast<const char*>(&nx), sizeof(nx));
                    fxz.write(reinterpret_cast<const char*>(&nz), sizeof(nz));
                }
                else {
                    fxy.open(("xy"+std::to_string(i)).c_str(), std::ifstream::out);
                    fxz.open(("xz"+std::to_string(i)).c_str(), std::ifstream::out);
                    fxy<<time<<" "<<nx<<" "<<ny<<std::endl;
                    fxz<<time<<" "<<nx<<" "<<nz<<std::endl;
                }
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

                        pjk.den = 0.0;
#ifdef GALPY                        
                        galpy_manager.calcAccPot(pjk.acc, pjk.pot, time, 0, pjk.pos, pjk.pos); 
                        for (int k=0; k<nset; k++) {
                            pjk.den += galpy_manager.calcSetDensity(k, time, pjk.pos, pjk.pos);
                        }
#elif AGAMA
                        agama_manager.calcAccPot(pjk.acc, pjk.pot, time, 0, pjk.pos, pjk.pos); 
#endif
                        if (out_binary) {
                            pjk.writeBinary(fxy);
                        }
                        else {
                            pjk.printColumn(fxy);
                            fxy<<std::endl;
                        }
                    }
                
                    for (int k=0; k<nz; k++) {
                        auto& pjk = particle_xz[j][k];
                        pjk.mass = 0;
                        pjk.pos[0] = x;
                        pjk.pos[1] = 0;
                        pjk.pos[2] = zmin + (zmax-zmin)/(nz-1)*k;
                        pjk.vel[0] = pjk.vel[1] = pjk.vel[2] = 0;

                        pjk.den = 0.0;
#ifdef GALPY                        
                        galpy_manager.calcAccPot(pjk.acc, pjk.pot, time, 0, pjk.pos, pjk.pos);  
                        for (int k=0; k<nset; k++) {
                            pjk.den += galpy_manager.calcSetDensity(k, time, pjk.pos, pjk.pos);
                        }
#elif AGAMA
                        agama_manager.calcAccPot(pjk.acc, pjk.pot, time, 0, pjk.pos, pjk.pos); 
#endif
                        if (out_binary) {
                            pjk.writeBinary(fxz);
                        }
                        else {
                            pjk.printColumn(fxz);
                            fxz<<std::endl;
                        }   
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

#ifdef GALPY        
        galpy_manager.initial(galpy_parameters, time, std::string(), false, true);
        int nset = galpy_manager.getNSet();
#elif AGAMA
        agama_manager.initial(agama_parameters, time, true);
#endif

        Particle particles[n];

        Particle::printColumnTitle(std::cout);
        std::cout<<std::endl;

        for (int i=0; i<n; i++) {
            auto& pi = particles[i];
            pi.readAscii(fp);
            double pos[3] = {pi.pos[0] + pos_offset[0],
                             pi.pos[1] + pos_offset[1],
                             pi.pos[2] + pos_offset[2]};
            pi.den = 0.0;
#ifdef GALPY
            galpy_manager.calcAccPot(pi.acc, pi.pot, time, 0, pos, &pi.pos[0]);
            for (int k=0; k<nset; k++) {
                pi.den += galpy_manager.calcSetDensity(k, time, pi.pos, pi.pos);}
#elif AGAMA
            agama_manager.calcAccPot(pi.acc, pi.pot, time, 0, pos, &pi.pos[0]);
#endif
            pi.printColumn(std::cout);
            std::cout<<std::endl;
        }
    }    
    return 0;
}
