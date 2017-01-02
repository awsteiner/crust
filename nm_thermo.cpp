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
#include "nm_thermo.h"

#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int nm_thermo::solve_fun(size_t nv, const ubvector &x, 
			 ubvector &y, matter &nm) {

  nm.n->n=x[0];
  nm.p->n=nm.nb-nm.n->n;
  
  if (nm.n->n<=0.0 || nm.p->n<=0.0 || 
      nm.n->n>nm.nb || nm.p->n>nm.nb) {
    return 1;
  }

  het->calc_temp_e(*nm.n,*nm.p,nm.T,nm.drip_th);
  
  nm.e.mu=nm.n->mu-nm.p->mu;
  relf.calc_mu(nm.e,nm.T);
  y[0]=(nm.p->n-nm.e.n)/nm.p->n;
  
  if (nm.e.n<1.0e-19) {
    return 2;
  }

  return 0;
}

void nm_thermo::calc(matter &nm) {

#ifdef O2SCL_NEVER_DEFINED  
  mm_funct_mfptr_param<nm_thermo,matter> fmp(this,&nm_thermo::solve_fun,nm);

  // Generate good guess for neutron density
  ovector x(1), y(1);
  x[0]=nm.nb*1.0e-1;
  int it=0;
  int ret=solve_fun(1,x,y,nm);
  while (it<20 && ret!=0) {
    if (ret==2) x[0]=fabs(x[0]+9.0*nm.nb)/10.0;
    else x[0]/=2.0;
    it++;
    ret=solve_fun(1,x,y,nm);
  }

  // If guess failed
  if (solve_fun(1,x,y,nm)!=0) {
    O2SCL_ERR("Function free_energy() failed.",o2scl::exc_efailed);
  }

  // Perform solution
  cmr2.msolve(1,x,fmp);
  solve_fun(1,x,y,nm);

  // Determine remaining quantities
  nm.rho=nm.e.n*nm.e.m+nm.p->m*nm.p->n+nm.n->m*nm.n->n;
  nm.rho=cng->convert("1/fm^4","g/cm^3",nm.rho);
  nm.ed=nm.drip_th.ed+nm.e.ed;
  nm.fr=nm.drip_th.ed+nm.e.ed-nm.T*(nm.drip_th.en+nm.e.en);
  nm.pr=nm.drip_th.pr+nm.e.pr;
  nm.gb=nm.fr+nm.pr;
  nm.mun=nm.n->mu;
  nm.mup=nm.p->mu;
#endif

  return;
}

void nm_thermo::calc_dist_x(matter &m, matter &nm) {
    
  vector<nucleus> &dist=m.dist;

  // Work at the same temperature
  nm.T=m.T;

  // We don't need to iterate to solve for the electron, neutron,
  // and proton densities, because even though the mass depends
  // on the number of electrons, the value of Np and Zp do not,
  // so we can temporary use 0.0 for the electron density in the
  // call to nucleus_be() below.
    
  // Add up neutron and proton densities 
  nm.n->n=m.n->n;
  nm.p->n=m.p->n;

  double nptot=0.0;
  for(size_t j=0;j<dist.size();j++) {
    nptot+=dist[j].n*dist[j].Z;
  }

  for(size_t j=0;j<dist.size();j++) {
    
    double Rws, chi;
    lda->nucleus_be(dist[j].Z,dist[j].N,m.p->n,m.n->n,nm.T,nptot,Rws,chi);

    // Compute Nprime
    double Rn3=lda->Rn*lda->Rn*lda->Rn;
    double Nprime=dist[j].N-4.0/3.0*pi*Rn3*m.n->n;

    nm.n->n+=dist[j].n*Nprime;
    nm.p->n+=dist[j].n*dist[j].Z;

  }
    
  // Ensure there's a good electron guess
  nm.e.mu=nm.e.m*1.0001;

  // Compute electron contribution
  nm.e.n=nm.p->n;
  nm.e.mu=nm.e.m*(1.0+1.0e-6);
  relf.calc_density(nm.e,nm.T);
    
  // Add the electron contribution to the free energy
  nm.fr=nm.e.ed-nm.T*nm.e.en;
    
  // Compute nucleonic contribution
  thermo drip_th;

  // Ensure a consistent guess for the chemical potentials
  nm.n->nu=nm.n->m*1.0001;
  nm.p->nu=nm.p->m*1.0001;

  if (nm.T<=0.0) {
    het->calc_e(*nm.n,*nm.p,drip_th);
    nm.fr+=drip_th.ed;
  } else {
    het->calc_temp_e(*nm.n,*nm.p,nm.T,drip_th);
    nm.fr+=drip_th.ed-nm.T*drip_th.en;
  }
  nm.pr=drip_th.pr+nm.e.pr;
  nm.gb=nm.fr+nm.pr;
  nm.mun=nm.n->mu;
  nm.mup=nm.p->mu;

  return;
}

int nm_thermo::check_free_energy_x(matter &m, double T, matter &nm) {
  
  m.T=T;
  calc_dist_x(m,nm);

  if (m.dist.size()!=1) {
    O2SCL_ERR("Wrong distribution size for check_free_energy_x.",
	      o2scl::exc_efailed);
  }
  cout.precision(10);
  cout << endl;
  cout << "Normal    : " << nm.fr << endl;
  double fr1=nm.fr;
      
  // Compute phi
  double Rws, chi;
  lda->nucleus_be(m.dist[0].Z,m.dist[0].N,m.p->n,m.n->n,T,m.e.n,Rws,chi);
  double phi=4.0/3.0*pi*pow(lda->Rn,3.0)*m.dist[0].n;
      
  nm.fr=0.0;
      
  // Add up neutron and proton densities
  nm.n->n=m.dist[0].N*m.dist[0].n;
  nm.p->n=m.dist[0].Z*m.dist[0].n;
  nm.n->n+=(1.0-phi)*m.n->n;
      
  // Ensure there's a good electron guess
  nm.e.mu=nm.e.m*1.0001;

  // Compute electron contribution
  nm.e.n=nm.p->n;
  relf.calc_density(nm.e,T);
  nm.fr+=nm.e.ed-T*nm.e.en;
      
  // Ensure a consistent guess for the chemical potentials
  nm.n->nu=nm.n->m*1.0001;
  nm.p->nu=nm.p->m*1.0001;

  // Compute nucleonic contribution
  thermo drip_th2;
  het->calc_temp_e(*nm.n,*nm.p,T,drip_th2);
  nm.fr+=drip_th2.ed-T*drip_th2.en;
  cout << "Alternate : " << nm.fr << endl;
    
  cout << "Difference: "
       << fabs(nm.fr-fr1)/fabs(nm.fr) << endl;

  test_mgr t;
  t.set_output_level(2);
  t.test_rel(nm.fr,fr1,1.0e-12,"free_energy_nm_x()");
  t.report();

  return 0;
}

