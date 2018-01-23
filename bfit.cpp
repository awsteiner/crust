/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2018, Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include <iostream>
#include <vector>

#include <o2scl/nucleus.h>
#include <o2scl/nuclear_mass.h>
#include <o2scl/nuclear_dist.h>
#include <o2scl/mass_fit.h>
#include <o2scl/skyrme_eos.h>
#include <o2scl/ldrop_mass.h>
#include <o2scl/table.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/apr_eos.h>
#include <o2scl/gsl_mmin_simp2.h>
#include <o2scl/cli_readline.h>
#include <o2scl/convert_units.h>
#include <o2scl/cern_minimize.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_fm;
using namespace o2scl_hdf;

class bfit {

public:
  
  static const size_t nfit=9;
  table bdat;
  string col;
  gsl_mmin_simp2<> gms;
  static const size_t ftype=1;

  int load() {
    hdf_file hf;
    hf.open("beard10.o2");
    hdf_input(hf,bdat);
    hf.close();
    return 0;
  }

  double min_func(size_t nv, const ovector_base &x) {
    double ret=0.0;
    for(size_t i=0;i<bdat.get_nlines();i++) {
      double exp=bdat.get(col,i);
      double A1=bdat.get("A1",i);
      double A2=bdat.get("A2",i);
      double r1=bdat.get("Z1",i)/A1;
      double r2=bdat.get("Z2",i)/A2;
      double apx;
      if (ftype==0) {
	apx==x[0]+x[1]*A1+x[2]*A1*A1+
	  x[3]*A2+x[4]*A2*A2+x[5]*r1*r1+x[7]*r2+x[8]*r2*r2+
	  x[6]*r1;
      } else {
	apx=x[0]*sqrt(A1)*sqrt(A2)+x[1]*pow(A1,x[2])+
	  x[3]*pow(A2,x[4])+x[5]*r1*r1+x[7]*r2+x[8]*r2*r2+
	  x[6]*r1;      
      }
      ret+=fabs((exp-apx)/exp);
    }
    ret/=((double)(bdat.get_nlines()));
    return ret;
  }

  int fit() {

    ovector x(nfit);
    
    multi_funct_mfptr<bfit> mfm(this,&bfit::min_func);
    double y;

    gms.verbose=1;
    gms.ntrial*=100;
    double step[1]={10.0};
    gms.set_step(1,step);

    col="B1";

    if (ftype==0) {
      x[0]=-1.834781e+02;
      x[1]=5.431898e+00;
      x[2]=-3.168715e-02;
      x[3]=2.174295e+00;
      x[4]=-2.429285e-03;
      x[5]=-1.634960e+02;
      x[6]=2.392794e+02;
      x[7]=1.930001e+02;
      x[8]=-1.540790e+02;
    } else {
      x[0]=4.664998e+00;
      x[1]=4.466440e+00;
      x[2]=6.728744e-01;
      x[3]=-1.801005e+02;
      x[4]=-3.867800e-02;
      x[5]=-1.366032e+02;
      x[6]=2.072307e+02;
      x[7]=1.582397e+02;
      x[8]=-1.029325e+02;
  
    }

    cout << min_func(nfit,x) << endl;
    
    col="B2";
    
    if (ftype==0) {
      x[0]=-2.516845e+00;
      x[1]=1.188301e-03;
      x[2]=-9.455966e-07;
      x[3]=4.394474e-03;
      x[4]=-3.404338e-05;
      x[5]=-2.418328e+00;
      x[6]=3.242123e+00;
      x[7]=3.332822e+00;
      x[8]=-2.572722e+00;
    } else {
      x[0]=1.188514e-03;
      x[1]=5.023013e+03;
      x[2]=-6.784830e+00;
      x[3]=-2.606291e+00;
      x[4]=-2.540973e-02;
      x[5]=-2.303624e+00;
      x[6]=3.125564e+00;
      x[7]=3.310691e+00;
      x[8]=-2.532734e+00;
    }

    cout << min_func(nfit,x) << endl;

    col="B3";
    
    if (ftype==0) {
      x[0]=-8.771906e-02;
      x[1]=6.020757e-04;
      x[2]=-5.082465e-06;
      x[3]=6.274142e-04;
      x[4]=-4.174494e-06;
      x[5]=-5.381484e-02;
      x[6]=8.802582e-02;
      x[7]=7.954748e-02;
      x[8]=-4.015003e-02;
    } else {
      x[0]=5.520281e-04;
      x[1]=-4.137702e+04;
      x[2]=-7.623811e+00;
      x[3]=-9.309219e-02;
      x[4]=-6.171452e-02;
      x[5]=-5.021372e-02;
      x[6]=8.359441e-02;
      x[7]=8.159595e-02;
      x[8]=-4.174526e-02;
    }

    cout << min_func(nfit,x) << endl;

    char ch2;
    cin >> ch2;

    gms.mmin(nfit,x,y,mfm);

    for(size_t i=0;i<5;i++) {
      char ch;
      cin >> ch;
      gms.mmin(nfit,x,y,mfm);
    }

    for(size_t i=0;i<nfit;i++) {
      cout << "x[" << i << "]=" << x[i] << ";" << endl;
    }
    return 0;
  }

  int run() {
    load();
    col="B1";
    fit();
    return 0;
  }

};

int main(void) {

  cout.setf(ios::scientific);

  bfit b;
  
  b.run();
  

  return 0;
}
