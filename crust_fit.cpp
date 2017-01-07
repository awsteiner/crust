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
using namespace crust;
using namespace o2scl_const;
using namespace o2scl_hdf;

crust_fit::crust_fit() {
  fit_moller=true;
  fit_dir=".";

  def_mmin.ntrial*=10;
  
  o2scl_hdf::ame_load(ame);
  o2scl_hdf::mnmsk_load(moller);
}

int crust_fit::read_fit(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cout << "Need filename for read_fit." << endl;
    return 1;
  }
  
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

  if (!fit_moller) {
    nucdist_set(dist,ame);
  } else {
    nucdist_set(dist,moller);
  }

  // Check fit
  double qual;
  eval(*lda,qual);
  cout << "Fit quality " << qual << "\n" << endl;

  return 0;
}

int crust_fit::perform_fit(std::vector<std::string> &sv, bool itive_com) {
    
  cout << "Performing mass fit." << endl;

  // Fit quality
  double qual;

  // Mass fitter object
  even_even=false;
  
  if (!fit_moller) {
    nucdist_set(dist,ame);
  } else {
    nucdist_set(dist,moller);
  }
  
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

  cout << "qual " << qual << endl;
  
  // Try fit a couple more times
  for(size_t i=0;i<5;i++) {
    fit(*lda,qual);
    lda->guess_fun(lda->nfit,x);
    cout << "qual " << qual << endl;
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
    
  return 0;
}
