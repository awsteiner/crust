/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2017, Andrew W. Steiner
  
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
#include "crust_fit.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

crust_fit::crust_fit() {
  loaded=false;
  fit_moller=true;
  fit_dir=".";
}

void crust_fit::eval(nucmass &n, double &fmin) {

  std::cout << "Missing 12." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED
  
  fmin=0.0;    
  size_t nn=0;
  for(nuclear_dist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {

    int Z=ndi->Z;
    int N=ndi->N;

    // Alfredo's masses
    if (0) {
      if (Z==21 && N==53-Z) ndi->mex=-38.15/hc_mev_fm;
      if (Z==21 && N==54-Z) ndi->mex=-33.59/hc_mev_fm;
      if (Z==21 && N==55-Z) ndi->mex=-30.32/hc_mev_fm;
      if (Z==22 && N==57-Z) ndi->mex=-33.82/hc_mev_fm;
      if (Z==22 && N==58-Z) ndi->mex=-29.72/hc_mev_fm;
      if (Z==23 && N==60-Z) ndi->mex=-33.05/hc_mev_fm;
      if (Z==23 && N==61-Z) ndi->mex=-30.56/hc_mev_fm;
      if (Z==24 && N==63-Z) ndi->mex=-35.24/hc_mev_fm;
      if (Z==25 && N==65-Z) ndi->mex=-40.73/hc_mev_fm;
      if (Z==25 && N==66-Z) ndi->mex=-36.87/hc_mev_fm;
      if (Z==26 && N==67-Z) ndi->mex=-45.88/hc_mev_fm;
      if (Z==26 && N==68-Z) ndi->mex=-44.01/hc_mev_fm;
      if (Z==27 && N==70-Z) ndi->mex=-46.72/hc_mev_fm;
      if (Z==27 && N==71-Z) ndi->mex=-44.45/hc_mev_fm;
      if (Z==28 && N==74-Z) ndi->mex=-49.18/hc_mev_fm;
      if (Z==29 && N==77-Z) ndi->mex=-46.97/hc_mev_fm;
    }

    if (N>=minN && Z>=minZ && (even_even==false || (N%2==0 && Z%2==0))) {
      fmin+=pow(ndi->mex*hc_mev_fm-n.mass_excess(Z,N),2.0);
      nn++;
    }
  }
  fmin=sqrt(fmin/nn);

  // If the result is not finite, just return an arbitrary large value
  if (!finite(fmin)) {
    fmin=1.0e4;
  }

#endif
  
  return;
}

void crust_fit::load(string dir) {
  o2scl_hdf::ame_load(ame,"",dir);
  o2scl_hdf::mnmsk_load(moller,dir);
  loaded=true;
  return;
}

int crust_fit::read_fit(std::vector<std::string> &sv, bool itive_com) {

  if (loaded==false) {
    o2scl_hdf::ame_load(ame,"","indata");
    o2scl_hdf::mnmsk_load(moller,"indata");
  }

  if (sv.size()<2) {
    cout << "Need filename for read_fit." << endl;
    return 1;
  }
  
  std::cout << "Missing 13." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED
  if (!fit_moller) {
    set_exp_mass(ame);
  } else {
    set_exp_mass(moller);
  }
#endif

  // Load initial guess from mass formula
  ubvector x(lda->nfit);
  lda->guess_fun(lda->nfit,x);

  // If present, load guess from file
  string fn=fit_dir+"/"+sv[1];
  cout << "Reading fit from file named '" << fn << "'." << endl;
  ifstream fin(fn.c_str());
  double temp;
  for(size_t i=0;i<lda->nfit && (fin >> temp);i++) {
    x[i]=temp;
  }
  fin.close();

  // Set parameters from vector 'x'.
  lda->fit_fun(lda->nfit,x);

  // Check fit
  double qual;
  eval(*lda,qual);
  cout << "Fit quality " << qual << "\n" << endl;

  return 0;
}

int crust_fit::perform_fit(std::vector<std::string> &sv, bool itive_com) {
    
  if (loaded==false) {
    o2scl_hdf::ame_load(ame,"","indata");
    o2scl_hdf::mnmsk_load(moller,"indata");
  }

  cout << "Performing mass fit." << endl;

  // Fit quality
  double qual;

  // Mass fitter object
  def_mmin.ntrial*=10;
  even_even=false;

  std::cout << "Missing 14." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED
  if (!fit_moller) {
    set_exp_mass(ame);
  } else {
    set_exp_mass(moller);
  }
#endif

  // Minimizer stepsize
  ubvector step(1);
  step[0]=1.0;
  def_mmin.set_step(1,step);
  def_mmin.verbose=0;
    
  // Load initial guess from mass formula
  ubvector x(lda->nfit);
  lda->guess_fun(lda->nfit,x);
    
  // If present, load guess from file
  if (sv.size()>=2) {
    ifstream fin(sv[1].c_str());
    double temp;
    for(size_t i=0;i<lda->nfit && (fin >> temp);i++) {
      if (i==0) {
	cout << "Reading fit from file named '" << sv[1] << "'." << endl;
      }
      x[i]=temp;
    }
    fin.close();
  }

  // Perform initial fit
  lda->fit_fun(lda->nfit,x);
  eval(*lda,qual);
  std::cout << "Missing 15." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED
  cout << "G " << x << "\n" << qual << endl;
  
  // Try fit a couple more times
  for(size_t i=0;i<3;i++) {
    fit(*lda,qual);
    lda->guess_fun(lda->nfit,x);
    cout << i << " " << x << "\n" << qual << endl;
  }
  cout << endl;

  // Output best parameters
  for(size_t j=0;j<x.size();j++) {
    cout << "    x[" << j << "]=" << x[j] << ";" << endl;
  }
  cout << endl;

  // If specified, store guess in file
  if (sv.size()>=2) {
    cout << "Writing fit to file named '" << sv[1] << "'." << endl;
    ofstream fout(sv[1].c_str());
    fout.setf(ios::scientific);
    fout.precision(10);
    for(size_t i=0;i<lda->nfit;i++) {
      fout << x[i] << endl;
    }
    fout << qual << endl;
    fout.close();
  }
#endif
    
  return 0;
}
