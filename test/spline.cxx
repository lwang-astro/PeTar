#include <iomanip>
#include <fstream>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main(int argc, char** argv)
{
  // data file name
  char* filename = argv[argc-1];
  // open data file
  std::fstream fs;
  fs.open(filename,std::fstream::in);
  if(!fs.is_open()) {
    std::cerr<<"Error: Filename "<<filename<<" not found\n";
    abort();
  }

  int N;
  fs>>N;
  double* x=new double[N];
  double* y=new double[N];
  
  for(int i=0;i<N;i++) {
      fs>>x[i]>>y[i];
      if (fs.eof()) {
          std::cerr<<"Error: data file reach end when reading pairs (current loaded pair number is "<<i<<"; required pair number "<<N-1<<std::endl;
          abort();
      }
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  const gsl_interp_type *t = gsl_interp_cspline_periodic; 
  gsl_spline *spline = gsl_spline_alloc (t, N);

  gsl_spline_init (spline, x, y, N);
  
  for (int i=0; i<=100; i++)
  {
      double xi = (1 - i / 100.0) * x[0] + (i / 100.0) * x[N-1];
      double yi = gsl_spline_eval (spline, xi, acc);
      printf ("%g %g\n", xi, yi);
  }
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  delete[] x;
  delete[] y;

  return 0;
}
