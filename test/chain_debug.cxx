#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<unistd.h>
#include<stdio.h>
#include<string>
#ifdef USE_C03
#include<map>
#else //USE_C03
#include<unordered_map>
#endif //USE_C03
#include<particle_simulator.hpp>
#include"AR.h"
#include"hard_force.hpp"
#include"kepler.hpp"
#include"hard.hpp"

#ifndef NAN_CHECK
#define NAN_CHECK(val) assert((val) == (val));
#endif

void chain_print(const ARC::chain<PtclHard> &c, const double ds, const double w, const double pre) {
  // printing digital precision
  std::cout<<std::setprecision(pre);
 
  std::cout<<c.getTime()
           <<std::setw(w)<<(c.getEkin()+c.getPot()+c.getPt())/c.getPt()
           <<std::setw(w)<<c.getEkin()
           <<std::setw(w)<<c.getPot()
           <<std::setw(w)<<c.getPt()
           <<std::setw(w)<<c.getw()
           <<std::setw(w)<<c.getW();
  const int n = c.getN();
  for (int j=0;j<n;j++) {
    std::cout<<std::setw(w)<<c.getP(j).getMass();
    for (int k=0;k<3;k++) {
      std::cout<<std::setw(w)<<c.getP(j).getPos()[k];
    }
    for (int k=0;k<3;k++) {
      std::cout<<std::setw(w)<<c.getP(j).getVel()[k];
    }
  }
  std::cout<<std::setw(w)<<ds;
  std::cout<<std::endl;
}  


int main(int argc, char **argv){
  PS::F64 ds=-1.0, toff=-1.0, rout=-1.0, rin=-1.0;
  PS::S32 n=0;
#ifndef ARC_SYM
  PS::S32 iter=-1,intp=-1;
#endif
  
  int copt;
  int cint=0;
  while ((copt = getopt(argc, argv, "s:n:t:r:R:i:I:h")) != -1)
    switch (copt) {
    case 's':
      ds = atof(optarg);
      cint++;
      break;
    case 'n':
      n = atoi(optarg);
      cint++;
      break;
    case 't':
      toff = atof(optarg);
      cint++;
      break;
    case 'r':
      rin = atof(optarg);
      cint++;
      break;
    case 'R':
      rout = atof(optarg);
      cint++;
      break;
#ifndef ARC_SYM
    case 'i':
      iter = atoi(optarg);
      cint++;
      break;
    case 'I':
      intp = atoi(optarg);
      cint++;
      break;
#endif
    case 'h':
      std::cout<<"options:\n"
               <<"    -s [double]:  step size ("<<ds<<")\n"
               <<"    -n [int]:     number of steps ("<<n<<")\n"
               <<"    -t [double]:  ending time ("<<toff<<")\n"
               <<"    -r [double]:  r_in ("<<rin<<")\n"
               <<"    -R [double]:  r_out ("<<rout<<")\n"
#ifndef ARC_SYM
               <<"    -i [int]:     itermax("<<iter<<")\n"
               <<"    -I [int]:     dense intpmax("<<intp<<")\n"
#endif
               <<"    -h :          help\n"
               <<std::endl;
      return 0;
    default:
      std::cerr<<"Unknown argument. check '-h' for help.\n";
      abort();
    }
  
//  if (argc==1) {
//    std::cerr<<"Please provide particle data filename\n";
//    abort();
//  }
      
  // data file name
  std::string filename="ARC_dump.dat";
  if (argc-cint*2>1) filename=argv[argc-1];
//  FILE* pf = fopen(filename,"r");
//  if (pf==NULL) {
//    std::cerr<<"Error: Filename "<<filename<<" not found\n";
//    abort();
//  }
  // open data file
//  std::fstream fs;
//  fs.open(filename,std::fstream::in);
//  if(!fs.is_open()) {
//    std::cerr<<"Error: Filename "<<filename<<" not found\n";
//    abort();
//  }
  
  // reading particle data
//  PtclHard *p=new PtclHard[n];
//  fread(p,sizeof(PtclHard),n,pf);
//  for (int i=0;i<n;i++) {
//    fs>>p[i].mass>>p[i].pos[0]>>p[i].pos[1]>>p[i].pos[2]>>p[i].vel[0]>>p[i].vel[0]>>p[i].vel[0];
//    if (fs.eof()) {
//      std::cerr<<"Error: data file reach end when reading particles (current loaded particle number is "<<i<<"; required N = "<<n<<std::endl;
//      abort();
//    }
//  }
  

//  const PS::F64 energy_error=1e-12;
//  const PS::F64 dterr=1e-10;
//  const PS::F64 dtmin=1e-24;
//  const PS::S32 exp_method=1;
//  const PS::S32 exp_itermax=20;
//  const PS::S32 exp_fix_iter=0;
  
  std::FILE* fp = std::fopen(filename.c_str(),"r");
  PS::F64 dtmp;
  size_t rc;
  if(toff>0) rc=fread(&dtmp, sizeof(PS::F64), 1, fp);
  else rc=fread(&toff, sizeof(PS::F64), 1, fp);
  if (rc<1) {
      std::cerr<<"Error: read time_end fails!\n";
      abort();
  }
  if(ds>0)  rc=fread(&dtmp, sizeof(PS::F64), 1, fp);
  else rc=fread(&ds, sizeof(PS::F64), 1, fp);
  if (rc<1) {
      std::cerr<<"Error: read ds fails!\n";
      abort();
  }

  ARC_pert_pars int_par; /// ARC integration parameters, rout_, rin_ (L.Wang)
  int_par.read(fp);

  PS::S32 np;
  rc=fread(&np, sizeof(PS::S32), 1, fp);
  if (rc<1) {
      std::cerr<<"Error: read np fails!\n";
      abort();
  }
  
  
  PtclH4* ptr[np+1];
  PtclForce* ftr[np+1];
  PtclH4 pert[np+1];
  PtclForce force[np+1];
  for (int i=0; i<np; i++) {
      pert[i].read(fp);
      force[i].read(fp);
      ptr[i]=&pert[i];
      ftr[i]=&force[i];
  }
 
//  chain_control.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
//  chain_control.read(controlname.c_str());
//  else {
//    chain_control.setabg(0,1,0);
//    chain_control.setEXP(energy_error,dtmin,dterr,exp_itermax,exp_method,3,(bool)exp_fix_iter);
//  }

//  chain_control.print(std::cerr);
  
  ARC::chain<PtclHard> c;
  
  //  c.addP(n,p);
  //  c.Int_pars=ARC_int_pars;
  //  c.init(0);
  c.read(fp);

  ARC::chainpars chain_control;

  chain_control.read(fp);

  if(np>0) chain_control.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_pert<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);
  else chain_control.setA(Newtonian_cut_AW<PtclHard,ARC_pert_pars>,Newtonian_extA_soft<PtclHard,PtclH4*,PtclForce*,ARC_pert_pars>,Newtonian_timescale<ARC_pert_pars>);

  std::fclose(fp);

  if (rout>0) int_par.rout = rout;
  if (rin>0)  int_par.rin = rin;
#ifndef ARC_SYM
  bool err_ignore=false;
  if (iter>0) {
    int seq=chain_control.getSeq();
    chain_control.setIterSeq(iter,seq);
    chain_control.setIterConst(true);
    err_ignore=true;
  }
  if (intp>0) chain_control.setDenIntpmax(intp);
#endif

  chain_control.print(std::cerr);
  std::cerr<<"-------------Int parameters:--------------\n";
  std::cerr<<"rout = "<<int_par.rout<<"; rin = "<<int_par.rin<<std::endl;
  std::cerr<<"time_end= "<<toff<<"; ds_use = "<<ds<<std::endl;

//  double nds = c.calc_next_step_custom(chain_control,&int_par);
//  std::cerr<<"New ds approx ="<<nds<<std::endl
//           <<"Used ds = "<<ds<<std::endl;

  c.init(0.0, chain_control, &int_par);

  c.print(std::cerr);

  int count=0;
#ifdef ARC_SYM
  if(n>0) {
      PS::F64 timetable[8];
#ifdef ARC_OPT_SYM2
      const PS::F64 m1=c.getP(0).getMass();
      const PS::F64 m2=c.getP(1).getMass();
      const PS::F64 m2_mt = m2/(m1+m2);
      const PS::F64 m1_m2_1 = -m1/m2-1.0;
#endif
      for(int i=0; i<n; i++) {
#ifdef ARC_OPT_SYM2
          c.Symplectic_integration_two(ds, chain_control, timetable, m2_mt, m1_m2_1, &int_par, ptr, ftr, np);
#else 
          c.Symplectic_integration(ds, chain_control, timetable, &int_par, ptr, ftr, np);
#endif
          chain_print(c,ds,24,16);
      }
  }
  else {
      bool fix_step_flag=false;
      if(c.slowdown.getkappa()>1.0&&c.getN()==2) fix_step_flag=true;
      count = c.Symplectic_integration_tsyn(ds, chain_control, toff, &int_par, ptr, ftr, np, fix_step_flag);
      chain_print(c,ds,24,16);
      std::cerr<<"Step count="<<count<<std::endl;
  }
#else
  while(toff-c.getTime()>chain_control.dterr||n>0) {
    PS::F64 dsf=c.extrapolation_integration<PtclH4*, PtclForce*, ARC_pert_pars>(ds,chain_control,toff,&int_par,ptr,ftr,np,err_ignore);
    //    std::cerr<<" Time="<<c.getTime()<<" Tdiff="<<toff-c.getTime()<<" ds="<<ds<<" dsf="<<dsf<<std::endl;
    chain_print(c,ds,24,16);

    count++;
    if (dsf<0) ds *= -dsf;
    else if(dsf==0) {
      c.info->ErrMessage(std::cerr);
      abort();
    }

    if (n>0&&count>=n) break;
  }
#endif

  return 0;
}
