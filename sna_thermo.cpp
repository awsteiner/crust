/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2020, Andrew W. Steiner
  
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
#include "sna_thermo.h"

#include <o2scl/test_mgr.h>

using namespace std;
using namespace crust;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int sna_thermo::free_energy_sna(matter &m, double T) {
    
  double rest=0.0;
  vector<nucleus> &dist=m.dist;

  // This doesn't work with free protons at the moment
  if (m.p->n!=0.0) {
    O2SCL_ERR2("Function sna_thermo::free_energy_sna() doesn't work with ",
	       "free protons.",o2scl::exc_eunimpl);
  }

  if (m.dist.size()!=1) {
    O2SCL_ERR2("Function sna_thermo::free_energy_sna() doesn't work with ",
	       "a distribution size not equal to 1.",o2scl::exc_efailed);
  }

  // Drip contribution 
  if (m.n->n>0.0) {
    thermo drip_th;

    // Ensure a consistent guess for the chemical potentials
    m.n->nu=m.n->m*1.0001;
    m.p->nu=m.p->m*1.0001;

    if (T<=0.0) {
      het->calc_e(*m.n,*m.p,drip_th);
      // Drip contribution 
      m.fr=drip_th.ed;
      part1=drip_th.ed;
      rest+=m.n->n*m.n->m;
    } else {
      het->calc_temp_e(*m.n,*m.p,T,drip_th);
      m.fr=drip_th.ed-T*drip_th.en;
      part1=drip_th.ed-T*drip_th.en;
      rest+=m.n->n*m.n->m;
    }
  } else {
    m.fr=0.0;
    part1=0.0;
  }

  // The total proton density
  double nptot=m.dist[0].n*m.dist[0].Z;
      
  // Ensure there's a good electron guess
  m.e.mu=m.e.m*1.0001;

  // Compute electron contribution
  m.e.n=nptot;
  relf.calc_density(m.e,T);
  rest+=m.e.n*m.e.m;
  
  // If there's a magnetic field, correct accordingly
  if (mag_field>1.0e-20) {
    elec_B.n=nptot;
    elec_B.mu=m.e.mu;
    mfz.calc_density_zerot_mag(elec_B,-mag_field*o2scl_const::ec_gauss_fm2_f<double>());
    m.e.ed=elec_B.ed;
    m.e.en=0.0;
    m.e.pr=elec_B.pr;
    m.e.ed=elec_B.ed;
    m.e.nu=elec_B.nu;
    m.e.mu=elec_B.mu;
  }
      
  // Add electron contribution to free energy
  m.fr+=m.e.ed-T*m.e.en;
  part2=m.e.ed-T*m.e.en;

  // Compute nuclear contribution
  double Rws, chi, be;
  be=lda->nucleus_be(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,T,m.e.n,
		     Rws,chi);
      
  // Compute the thermodynamics of the nucleus, without the rest
  // mass which will be added in later
  if ((m.dist[0].Z+m.dist[0].N)%2==0) m.dist[0].g=1.0;
  else m.dist[0].g=3.0;
  double mass_neutron=cng.convert("kg","1/fm",o2scl_const::mass_neutron_f<double>());
  double mass_proton=cng.convert("kg","1/fm",o2scl_const::mass_proton_f<double>());
    
  m.dist[0].m=m.dist[0].Z*mass_proton+m.dist[0].N*mass_neutron+be;
  m.dist[0].non_interacting=true;
  m.dist[0].inc_rest_mass=false;
  cla.calc_density(m.dist[0],T);

  // Compute Nprime
  double Rn3=lda->Rn*lda->Rn*lda->Rn;
  double Nprime=m.dist[0].N-4.0/3.0*pi*Rn3*m.n->n;
  
  // Add nuclear contribution to free energy, including the rest
  // mass energy density. Note that only N' and Z' are used to
  // compute the nuclear rest mass, since the remaining part of
  // the rest mass is in the drip.
  m.fr+=m.dist[0].n*(be+m.dist[0].Z*mass_proton+Nprime*mass_neutron);
  if (inc_nuc_trans) {
    m.fr+=m.dist[0].ed-T*m.dist[0].en;
  }
  rest+=m.dist[0].n*(m.dist[0].Z*mass_proton+Nprime*mass_neutron);
  part3=m.dist[0].n*be;
  part4=m.dist[0].n*(m.dist[0].Z*mass_proton+Nprime*mass_neutron);

  m.rho=cng.convert("1/fm^4","g/cm^3",rest);

  return 0;
}

int sna_thermo::free_energy_sna_fix_nb_nn(double nb, matter &m, double T) {
    
  // Compute initial guess for the number density of nuclei by first
  // computing the neutron radius
  lda->drip_binding_energy_full_d(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,0.0,T);
  m.dist[0].n=(nb-m.n->n)/(m.dist[0].Z+m.dist[0].N);
  if (lda->exc_volume) {
    m.dist[0].n=(nb-m.n->n)/(m.dist[0].Z+m.dist[0].N-
			     4.0/3.0*pi*pow(lda->Rn,3.0)*m.n->n);
  }
    
  baryon_density_sna(m,T);
  double nbx=m.nb;
  size_t i;
  for(i=0;i<40 && fabs((nb-nbx)/nb)>1.0e-10;i++) {
    m.dist[0].n+=(nb-nbx)/(m.dist[0].Z+m.dist[0].N);
    baryon_density_sna(m,T);
    nbx=m.nb;
  }
  if (false && fabs((nb-nbx)/nb)>1.0e-10) {
    cout << "Function free_energy_sna_fix_nb() failed." << endl;
    cout << nb << " " << nbx << " " << nb-nbx << " "
	 << (nb-nbx)/nb << endl;
    exit(-1);
  }

  free_energy_sna(m,T);

  return 0;
}

int sna_thermo::free_energy_sna_fix_nb_nnuc(double nb, matter &m, double T) {
    
  // Compute initial guess for the number density of nuclei by first
  // computing the neutron radius
  lda->drip_binding_energy_full_d(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,0.0,T);
  m.n->n=nb-m.dist[0].n*(m.dist[0].Z+m.dist[0].N);
  //if (lda->exc_volume) {
  //m.dist[0].n=(nb-m.n->n)/(m.dist[0].Z+m.dist[0].N-
  //4.0/3.0*pi*pow(lda->Rn,3.0)*m.n->n);
  //}
    
  baryon_density_sna(m,T);
  double nbx=m.nb;
  size_t i;
  for(i=0;i<40 && fabs((nb-nbx)/nb)>1.0e-10;i++) {
    m.n->n+=nb-m.nb;
    baryon_density_sna(m,T);
    nbx=m.nb;
  }

  free_energy_sna(m,T);

  return 0;
}

int sna_thermo::check_free_energy_sna(dist_thermo &dt, test_mgr &t) {
  
  cout.precision(10);

  double T=cng.convert("K","1/fm",2.0e8);

  matter m;
  m.dist.resize(1);
  m.dist[0].Z=40;
  m.dist[0].N=60;
  m.dist[0].n=0.0001;
  m.n->n=0.01;
  m.p->n=0.0;
  vector<nucleus> &dist=m.dist;
  
  double rest=0.0, be=0.0, Rws=0.0, chi=0.0;
  
  free_energy_sna(m,T);
  double fr1=m.fr;
  dt.free_energy_dist(m,T);
  double fr2=m.fr;

  dt.mass_density(m,T);

  double rest1=m.rho;

  m.fr=m.e.ed-T*m.e.en;
  rest=m.e.n*m.e.m;

  lda->exc_volume=false;

  be=lda->nucleus_be(dist[0].Z,dist[0].N,m.p->n,m.n->n,T,m.e.n,Rws,chi);
      
  double phi=4.0/3.0*pi*pow(lda->Rn,3.0)*dist[0].n;
  
  if ((dist[0].Z+dist[0].N)%2==0) dist[0].g=1.0;
  else dist[0].g=3.0;
  double mass_neutron=cng.convert("kg","1/fm",o2scl_const::mass_neutron_f<double>());
  double mass_proton=cng.convert("kg","1/fm",o2scl_const::mass_proton_f<double>());
  dist[0].m=dist[0].Z*mass_proton+dist[0].N*mass_neutron+be;
  dist[0].non_interacting=true;
  dist[0].inc_rest_mass=false;
  cla.calc_density(dist[0],T);
  m.fr+=dist[0].n*(be+dist[0].Z*mass_proton+dist[0].N*mass_neutron);
  if (inc_nuc_trans) {
    dist[0].ed-T*dist[0].en;
  }
  rest+=dist[0].n*(dist[0].Z*mass_proton+dist[0].N*mass_neutron);

  {
    // Compute neutron drip free energy
    thermo drip_th;

    // Ensure a consistent guess for the chemical potentials
    m.n->nu=m.n->m*1.0001;
    m.p->nu=m.p->m*1.0001;
    
    het->calc_temp_e(*m.n,*m.p,T,drip_th);
    m.fr+=(1.0-phi)*(drip_th.ed-T*drip_th.en);
    rest+=(1.0-phi)*(m.n->n*m.n->m);
  }

  m.rho=cng.convert("1/fm^4","g/cm^3",rest);

  t.test_rel(m.rho,rest1,1.0e-15,"free_energy_sna() 1");
  t.test_rel(m.fr,fr1,4.0e-14,"free_energy_sna() 2");

  if (!t.report()) {
    exit(-1);
  }

  return 0;
}

int sna_thermo::baryon_density_sna(matter &m, double T) {
    
  vector<nucleus> &dist=m.dist;
    
  // Compute total proton density
  double nptot=0.0;
  for(size_t j=0;j<dist.size();j++) {
    nptot+=dist[j].n*dist[j].Z;
  }

  double nbtot=m.n->n+m.p->n;

  // Add nuclear contribution to the baryon density
  double Rws, chi, be;
  be=lda->nucleus_be(dist[0].Z,dist[0].N,m.p->n,m.n->n,T,m.e.n,
		     Rws,chi);
  
  // Compute Nprime
  double Nprime=dist[0].N-4.0/3.0*pi*lda->Rn*lda->Rn*lda->Rn*m.n->n;
  
  nbtot+=dist[0].n*(dist[0].Z+Nprime);
  
  m.nb=nbtot;

  return 0;
}

double sna_thermo::get_phi_sna(matter &m, double T) {
    
  vector<nucleus> &dist=m.dist;
    
  // Compute total electron density
  double ne=0.0;
  for(size_t j=0;j<dist.size();j++) {
    ne+=dist[j].n*dist[j].Z;
  }
  
  double phi=0.0;
  double Rws, chi_i;
  lda->nucleus_be(dist[0].Z,dist[0].N,m.p->n,m.n->n,T,ne,Rws,chi_i);
  
  phi+=dist[0].n*4.0/3.0*pi*pow(lda->Rn,3.0);
  
  return phi;
}


