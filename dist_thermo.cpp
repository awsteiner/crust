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
#include "dist_thermo.h"

#include <o2scl/test_mgr.h>
#include <o2scl/mroot_cern.h>

using namespace std;
using namespace o2scl;
using namespace crust;
using namespace o2scl_const;
using namespace o2scl_hdf;

int dist_thermo::free_energy_dist_alpha(matter &m, double &T, double alpha) {
  free_energy_dist(m,T,false);
  double phi=get_phi(m,T);
  double alpha2=m.n->n*(1.0-phi);

  size_t i;
  for(i=0;i<20 && fabs(alpha2-alpha)/fabs(alpha)>1.0e-8;i++) {
    m.n->n=alpha/(1.0-phi);
    free_energy_dist(m,T,false);
    phi=get_phi(m,T);
    alpha2=m.n->n*(1.0-phi);
  }
  if (i==20) {
    O2SCL_ERR2("Iteration failed in ",
	       "dist_thermo::free_energy_dist_alpha().",o2scl::exc_efailed);
  }

  return 0;
}

int dist_thermo::free_energy_dist(matter &m, double T, bool eval_chem_pots) {

  bool inc_nuc_trans=false;

  if (!finite(m.n->n) || !finite(m.e.n)) {
    cout << "T: " << T << endl;
    cout << m << endl;
    O2SCL_ERR("Not finite in free_energy_dist() [point 1].",
	      o2scl::exc_efailed);
  }

  vector<nucleus> &dist=m.dist;

  // Compute the rest mass energy density as well
  m.rho=0.0;

  m.fr=0.0;

  if (m.n->n>0.0 || m.p->n>0.0) {
      
    // Ensure a consistent guess for the chemical potentials
    m.n->nu=m.n->m*1.0001;
    m.p->nu=m.p->m*1.0001;

    // Compute drip free energy (with rest mass contribution)
    if (T<=0.0) {
      het->calc_e(*m.n,*m.p,m.drip_th);
      m.fr+=m.drip_th.ed;
      part1=m.drip_th.ed;
    } else {
      het->calc_temp_e(*m.n,*m.p,T,m.drip_th);
      m.fr+=m.drip_th.ed-T*m.drip_th.en;
      part1=m.drip_th.ed-T*m.drip_th.en;
    }

    // Contribution to rest mass energy density
    m.rho=m.n->n*m.n->m;

    // Contribution to baryon density
    m.nb=m.n->n+m.p->n;

  } else {
    m.drip_th.ed=0.0;
    m.drip_th.en=0.0;
    m.nb=0.0;
    m.rho=0.0;
    part1=0.0;
  }
      
  // Compute total proton density
  double nptot=m.p->n;
  for(size_t i=0;i<dist.size();i++) {
    nptot+=dist[i].n*dist[i].Z;
  }

  // Compute electron contribution
  m.e.mu=m.e.m+1.0e-10;
  m.e.n=nptot;
  relf.calc_density(m.e,T);

  // If there's a magnetic field, correct the electrons accordingly
  if (mag_field>1.0e-20) {
    elec_B.n=nptot;
    elec_B.mu=m.e.mu;
    mfz.calc_density_zerot_mag(elec_B,-mag_field*o2scl_const::ec_gauss_fm2);
    m.e.ed=elec_B.ed;
    m.e.en=0.0;
    m.e.pr=elec_B.pr;
    m.e.ed=elec_B.ed;
    m.e.nu=elec_B.nu;
    m.e.mu=elec_B.mu;
  }

  m.fr+=m.e.ed-T*m.e.en;
  part2=m.e.ed-T*m.e.en;

  // Electron contribution to rest mass energy density
  m.rho+=m.e.m*m.e.n;

  part3=0.0;
  part4=0.0;
    
  // Compute nuclear contributions
  for(size_t i=0;i<dist.size();i++) {
      
    if (dist[i].n<0.0) {
      cout << "T: " << T << endl;
      cout << m << endl;
      O2SCL_ERR("Density negative in free_energy_dist().",o2scl::exc_efailed);
    }
    double Rws, chi;
    double be;
      
    if (!finite(m.n->n) || !finite(m.e.n)) {
      cout << "T: " << T << endl;
      cout << m << endl;
      O2SCL_ERR("Not finite in free_energy_dist() [point 2].",
		o2scl::exc_efailed);
    }
      
    //if (use_pasta) {
    //be=lda->nucleus_be_pasta(dist[i].Z,dist[i].N,m.p->n,
    //m.n->n,T,m.e.n,Rws,chi);
    //} else {
    be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
		       Rws,chi);
    //}
    
    // Compute the thermodynamics of the nucleus, without the rest
    // mass which will be added in later
    if ((dist[i].Z+dist[i].N)%2==0) dist[i].g=1.0;
    else dist[i].g=3.0;
    double mass_neutron=cng.convert
      ("kg","1/fm",o2scl_mks::mass_neutron);
    double mass_proton=cng.convert
      ("kg","1/fm",o2scl_mks::mass_proton);
    dist[i].m=dist[i].Z*mass_proton+dist[i].N*mass_neutron+be;
    dist[i].non_interacting=true;
    dist[i].inc_rest_mass=false;

    if (eval_chem_pots==true) {
      if (dist[i].Z==1 && dist[i].N==0) {
	m.zeta[i]=0.0;
      } else if (dist[i].Z==2 && dist[i].N==2) {
	m.zeta[i]=dist[i].N/(dist[i].Z*0.08);
      } else {
	m.zeta[i]=dist[i].N/(dist[i].Z*lda->nn+dist[i].Z*chi*lda->df_dchi);
      }
      //cout << "Here: " << i << " " << lda->nn << " " << chi << " "
      //<< lda->df_dchi << " " << m.zeta[i] << " " << lda->Rn << endl;
    }
      
    // Compute Nprime
    double Rn3=lda->Rn*lda->Rn*lda->Rn;
    double Rws3=Rws*Rws*Rws;
    double Nprime=dist[i].N-4.0/3.0*pi*Rn3*m.n->n;

    if (check==check_free_energy_dist) {
      cout.precision(10);
      cout << "i,exc,N': ";
      cout.width(2);
      cout << i << " " << lda->exc/hc_mev_fm << " " << Nprime << endl;
    }

    // Contribution to baryon density
    m.nb+=dist[i].n*(dist[i].Z+Nprime);

    // Correct the density by the appropriate factor first
    double fact=Rws3/(Rws3-Rn3);
    dist[i].n*=fact;

    cla.calc_density(dist[i],T);

    // Return density to normal value
    dist[i].n/=fact;
      
    // Correct translational energy
    double Tc=17.0/hc_mev_fm;
    double Wtrans=8.0/hc_mev_fm;
    double A0trans=12.0;
    double z=T/Tc;
    double gamma=4.0;
    double ftrans;

    if (inc_nuc_trans && z<1.0) {
      // FIXME is dist[i].n correct in the second term here?
      ftrans=dist[i].ed-T*dist[i].en-
	dist[i].n*T*T/Wtrans*(1.0+A0trans/(dist[i].N+dist[i].Z));
      ftrans*=(1.0-exp(-1.0*gamma*pow(1.0-z*z,2.0)))/(1.0-exp(-gamma));
    } else {
      ftrans=0.0;
    }
      
    // Add nuclear contribution to free energy, including the rest
    // mass energy density. Note that only N' and Z' are used to
    // compute the nuclear rest mass, since the remaining part of
    // the rest mass is in the drip.
    m.fr+=dist[i].n*(be+dist[i].Z*mass_proton+Nprime*mass_neutron)+ftrans;
    m.rho+=dist[i].n*(dist[i].Z*mass_proton+Nprime*mass_neutron);
    part3+=dist[i].n*be;
    part4+=dist[i].n*(dist[i].Z*mass_proton+Nprime*mass_neutron);
    
  }

  // Convert the result to g/cm^3
  m.rho=cng.convert("1/fm^4","g/cm^3",m.rho);

  if (eval_chem_pots) {

    // Perform second pass for chemical potentials

    double geo=4.0/3.0*pi, phi=0.0;

    // The infinite part of the neutron chemical potential

    m.mun=m.n->mu;
    
    ubvector dnnout_dni(dist.size());

    // Initialize dnnout_dni to zero
    for(size_t j=0;j<dist.size();j++) {
      dnnout_dni[j]=0.0;
    }
    
    // First pass: compute dnnout_dni, the nuclear part of the
    // neutron chemical potential, and phi

    for(size_t j=0;j<dist.size();j++) {

      // Compute the nuclear derivatives
      double Rws, chi, be;
      be=lda->nucleus_be(dist[j].Z,dist[j].N,m.p->n,m.n->n,T,m.e.n,
			 Rws,chi);
      double Rn3=lda->Rn*lda->Rn*lda->Rn;

      // Contribution to phi
      phi+=dist[j].n*geo*Rn3;
      
      // The contribution from each nucleus to the neutron 
      // chemical potential
      m.mun+=dist[j].n*((lda->dexc_dnn+lda->dshell_dnn)/hc_mev_fm-
			geo*Rn3*m.n->m);

      // Compute each contribution to dnnout_dni[i]
      for(size_t i=0;i<dist.size();i++) {
	double delta_ij=0.0;
	if (i==j) delta_ij=1.0;
	double dchij_dni=m.zeta[j]*dist[i].Z;
	double dphij_dni=dist[j].N/lda->nn*
	  (delta_ij-dist[j].n/lda->nn*lda->df_dchi*dchij_dni);
	dnnout_dni[i]+=dphij_dni;
      }

    }
    
    // Finish dnnout_dni and initialize nuclear chemical potential
    // to zero

    for(size_t i=0;i<dist.size();i++) {
      dnnout_dni[i]*=m.n->n/(1.0-phi);
      m.mu[i]=0.0;
    }

    // Compute the part of the nuclear chemical potential inside 
    // the sum
    
    for(size_t j=0;j<dist.size();j++) {

      // Compute the nuclear derivatives
      double Rws, chi, be;
      be=lda->nucleus_be(dist[j].Z,dist[j].N,m.p->n,m.n->n,T,m.e.n,
			 Rws,chi);
      double Rn3=lda->Rn*lda->Rn*lda->Rn;
      double dnuc_dchij=(lda->dbulk_dchi+lda->dsurf_dchi+lda->dcoul_dchi+
			 lda->dexc_dchi+lda->dshell_dchi)/hc_mev_fm;
      double dnuc_dnnout=(lda->dexc_dnn+lda->dshell_dnn)/hc_mev_fm;

      for(size_t i=0;i<dist.size();i++) {
	double dchij_dni=m.zeta[j]*dist[i].Z;
	double dNprime_dni=-4.0*pi*lda->Rn*lda->Rn*
	  (lda->Rn/3.0*dnnout_dni[i]+m.n->n*lda->dRn_dchi*dchij_dni);
	m.mu[i]+=dist[j].n*(dnuc_dchij*dchij_dni+
			    dnuc_dnnout*dnnout_dni[i]+m.n->m*dNprime_dni);
	if (!finite(m.mu[i])) {
	  cout << "i,j,distj.n,dnuc_dchij,dchij_dni: " 
	       << i << " " << j << " " << dist[j].n << " " 
	       << dnuc_dchij << " " << dchij_dni << endl;
	  cout << "zeta,Z: " << m.zeta[j] << " " << dist[i].Z << endl;
	  cout << "bulk,surf,coul,exc,shell: " 
	       << lda->dbulk_dchi << " " << lda->dsurf_dchi << " "
	       << lda->dcoul_dchi << " " << lda->dexc_dchi << " "
	       << lda->dshell_dchi << endl;
	  O2SCL_ERR2("Chemical potentials not finite at point 2 in ",
		     "free_energy_dist().",o2scl::exc_efailed);
	}
      }
    }
    
    // Assemble the full chemical potential for nuclei

    for(size_t i=0;i<dist.size();i++) {
      
      double Rws, chi, be;
      be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
			 Rws,chi);
      
      // Compute Nprime
      double Rn3=lda->Rn*lda->Rn*lda->Rn;
      double Nprime=dist[i].N-geo*Rn3*m.n->n;
      
      double mass_neutron=cng.convert
	("kg","1/fm",o2scl_mks::mass_neutron);
      double mass_proton=cng.convert
	("kg","1/fm",o2scl_mks::mass_proton);
      m.mu[i]+=m.n->mu*dnnout_dni[i]+dist[i].Z*m.e.mu+
	be+dist[i].Z*mass_proton+Nprime*mass_neutron;

      if (!finite(m.mu[i])) {
	cout << "i,m.n.mu,dnnout_dni,Z: " << i << " " 
	     << m.n->mu << " " << dnnout_dni[i] << " " << dist[i].Z << endl;
	cout << "mue,be,Nprime,mui: " << m.e.mu << " " << be << " " 
	     << Nprime << " " << m.mu[i] << endl;
	O2SCL_ERR2("Chemical potentials not finite at point 3 in ",
		   "free_energy_dist().",o2scl::exc_efailed);
      }
    }

  }
    
  if (check==check_free_energy_dist) {
      
    cout << endl;

    // Store old version and convert back to 1/fm^4
    double rest1=cng.convert("g/cm^3","1/fm^4",m.rho);
    double fr1=m.fr;
      
    m.rho=0.0;
    m.fr=m.e.ed-T*m.e.en;
      
    lda->exc_volume=false;

    double phi=0.0;

    for(size_t i=0;i<dist.size();i++) {
	
      double Rws, chi;
      double be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
				Rws,chi);
	
      phi+=4.0/3.0*pi*pow(lda->Rn,3.0)*dist[i].n;
      cout << "i, N, N': ";
      cout.width(2);
      cout << i << " ";
      cout.width(2);
      cout << dist[i].N << " "
	   << dist[i].N-4.0/3.0*pi*pow(lda->Rn,3.0)*m.n->n << endl;
      
      if ((dist[i].Z+dist[i].N)%2==0) dist[i].g=1.0;
      else dist[i].g=3.0;
      double mass_neutron=cng.convert
	("kg","1/fm",o2scl_mks::mass_neutron);
      double mass_proton=cng.convert
	("kg","1/fm",o2scl_mks::mass_proton);
      dist[i].m=dist[i].Z*mass_proton+dist[i].N*mass_neutron+be;
      dist[i].non_interacting=true;
      dist[i].inc_rest_mass=false;
      cla.calc_density(dist[i],T);
      m.fr+=dist[i].n*(be+dist[i].Z*mass_proton+dist[i].N*mass_neutron)+
	dist[i].ed-T*dist[i].en;
      m.rho+=dist[i].n*(dist[i].Z*mass_proton+dist[i].N*mass_neutron);
    }

    {
      // Compute neutron drip free energy
      thermo drip_th;

      // Ensure a consistent guess for the chemical potentials
      m.n->nu=m.n->m*1.0001;
      m.p->nu=m.p->m*1.0001;

      het->calc_temp_e(*m.n,*m.p,T,drip_th);
      m.fr+=(1.0-phi)*(drip_th.ed-T*drip_th.en);
      m.rho+=(1.0-phi)*(m.n->n*m.n->m);
    }
    
    cout << endl;

    cout << "Normal: rest, fr: " << rest1 << " " << fr1 << endl;
    cout << "Alt   : rest, fr: " << m.rho << " " << m.fr << endl;
    cout << "Difference      : " 
	 << fabs(m.rho-rest1)/fabs(m.rho) << " "
	 << fabs(m.fr-fr1)/fabs(m.fr) << endl;

    test_mgr t;
    t.set_output_level(2);
    t.test_rel(m.rho,rest1,4.0e-6,"free_energy_dist() 1");
    t.test_rel(m.fr,fr1,3.0e-8,"free_energy_dist() 2");
    t.report();
    
  }

  return 0;
}

int dist_thermo::mass_density(matter &m, double T) {
    
  vector<nucleus> &dist=m.dist;

  m.rho=0.0;

  // Compute neutron drip mass density
  if (m.n->n>0.0) {
    m.rho+=m.n->n*m.n->m;
  }
      
  // Compute total proton density to obtain electron density
  double nptot=m.p->n;
  for(size_t i=0;i<dist.size();i++) {
    nptot+=dist[i].Z*dist[i].n;
  }

  // Ensure there's a good electron guess
  m.e.mu=m.e.m*1.0001;

  // Compute electron contribution
  m.e.n=nptot;
  relf.calc_density(m.e,T);

  // If there's a magnetic field, correct accordingly
  if (mag_field>1.0e-20) {
    elec_B.n=nptot;
    elec_B.mu=m.e.mu;
    mfz.calc_density_zerot_mag(elec_B,-mag_field*o2scl_const::ec_gauss_fm2);
    m.e.ed=elec_B.ed;
    m.e.en=0.0;
    m.e.pr=elec_B.pr;
    m.e.ed=elec_B.ed;
    m.e.nu=elec_B.nu;
    m.e.mu=elec_B.mu;
  }

  // Add electron contribution
  m.rho+=m.e.n*m.e.m;
      
  // Add nuclear contribution to the rest mass density
  for(size_t i=0;i<dist.size();i++) {

    double Rws, chi, be;
    //if (use_pasta) {
    //be=lda.nucleus_be_pasta(dist[i].Z,dist[i].N,m.p->n,
    //m.n->n,T,m.e.n,Rws,chi);
    //} else {
    be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
		       Rws,chi);
    //}

    // Compute Nprime
    double Nprime=dist[i].N-4.0/3.0*pi*lda->Rn*lda->Rn*lda->Rn*m.n->n;
      
    double mass_neutron=cng.convert("kg","1/fm",o2scl_mks::mass_neutron);
    double mass_proton=cng.convert("kg","1/fm",o2scl_mks::mass_proton);
    m.rho+=dist[i].n*(dist[i].Z*mass_proton+Nprime*mass_neutron);
  }
    
  // Convert the result to g/cm^3
  m.rho=cng.convert("1/fm^4","g/cm^3",m.rho);

  if (check==check_mass_density) {

    cout.precision(10);
    cout << "rho: " << m.rho << endl;
    double temp=m.rho;
      
    // Compute electron contribution
    m.rho=m.e.n*m.e.m;

    double phi=0.0;

    // Add nuclear contribution to the rest mass density
    for(size_t i=0;i<dist.size();i++) {
      double Rws, chi;
      double be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
				Rws,chi);
      phi+=4.0/3.0*pi*pow(lda->Rn,3.0)*dist[i].n;
      double mass_neutron=cng.convert("kg","1/fm",o2scl_mks::mass_neutron);
      double mass_proton=cng.convert("kg","1/fm",o2scl_mks::mass_proton);
      m.rho+=dist[i].n*(dist[i].Z*mass_proton+dist[i].N*mass_neutron);
    }
      
    // Compute neutron drip mass density
    if (m.n->n>0.0) {
      m.rho+=(1.0-phi)*m.n->n*m.n->m;
    }
      
    // Convert the result to g/cm^3
    m.rho=cng.convert("1/fm^4","g/cm^3",m.rho);

    cout << "rho: " << m.rho << endl;
    cout << "Difference: "
	 << fabs(temp-m.rho)/m.rho << endl;

    test_mgr t;
    t.set_output_level(2);
    t.test_rel(temp,m.rho,1.0e-10,"mass_density");
    t.report();

    return 0;
  }

  return 0;
}

int dist_thermo::check_pressure(matter &m, double T) {

  test_mgr t;
  t.set_output_level(1);

  lda->use_ame=false;

  double mu1, mu2, mu3, mu4, pr1, gb1;
  deriv_gsl<free_dist_deriv_alpha> gd;
  double xn, xnuc, xnuc2, xnuc3;

  m.dist.clear();
  nucleus nuc;
  nuc.Z=32;
  nuc.N=126;
  nuc.n=1.0e-5;
  m.dist.push_back(nuc);
  m.n->n=7.0e-3;

  cout.precision(6);
  cout << "------------------------------------------------------" << endl;
  cout << "Test distribution with one nucleus:\n" << endl;
  cout << "Z,N,T: " << m.dist[0].Z << " " << m.dist[0].N << " " 
       << T << endl;
  cout << endl;

  gibbs_energy_dist(m,T);
  
  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  mu1=m.mun;
  mu2=m.mu[0];
  pr1=m.pr;
  gb1=m.gb;
  
  gd.h=1.0e-4;
  free_dist_deriv_alpha fda(*this,m,T);
  double xx=m.n->n*(1.0-get_phi(m,T)), der, dere;
  gd.deriv_err(xx,fda,der,dere);
  cout << "alpha,dfda,err: " << xx << " " << der << " " << dere << endl;
  m.mun=der*(1.0-get_phi(m,T));
  cout << "dfdn: " << m.mun << endl;

  // Return to original values
  m.n->n=7.0e-3;
  m.dist[0].n=1.0e-5;
  
  deriv_gsl<free_dist_deriv_nuc> gd2;
  gd2.h=1.0e-8;
  free_dist_deriv_nuc fdi(*this,m,T,m.n->n*(1.0-get_phi(m,T)),0);
  xx=m.dist[0].n;
  gd2.deriv_err(xx,fdi,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[0]=der;
  cout << endl;
  
  // Return to original values
  m.n->n=7.0e-3;
  m.dist[0].n=1.0e-5;
  free_energy_dist(m,T);

  t.test_rel(mu1,m.mun,1.0e-10,"mu1");
  t.test_rel(mu2,m.mu[0],1.0e-10,"mu2");

  // Combine to form the gibbs energy density and pressure
  m.gb=m.n->n*m.mun;
  for(size_t j=0;j<m.dist.size();j++) {
    m.gb+=m.dist[j].n*m.mu[j];
  }
  m.pr=m.gb-m.fr;

  t.test_rel(m.pr,pr1,1.0e-9,"pr1");
  t.test_rel(m.gb,gb1,1.0e-10,"gb1");

  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  cout << "------------------------------------------------------" << endl;
  cout << "Test distribution with one nucleus:\n" << endl;

  xn=5.0e-13;
  xnuc=3.0e-9;

  m.dist[0].Z=28;
  m.dist[0].N=36;
  m.n->n=xn;
  m.dist[0].n=xnuc;

  cout << "Z,N,T: " << m.dist[0].Z << " " << m.dist[0].N << " " 
       << T << endl;
  cout << endl;

  gibbs_energy_dist(m,T);
  
  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  mu1=m.mun;
  mu2=m.mu[0];
  pr1=m.pr;
  gb1=m.gb;
  
  gd.h=1.0e-14;
  xx=m.n->n*(1.0-get_phi(m,T));
  gd.deriv_err(xx,fda,der,dere);
  cout << "alpha,dfda,err: " << xx << " " << der << " " << dere << endl;
  m.mun=der*(1.0-get_phi(m,T));
  cout << "dfdn: " << m.mun << endl;

  // Return to original values
  m.n->n=xn;
  m.dist[0].n=xnuc;
  
  gd2.h=1.0e-10;
  xx=m.dist[0].n;
  free_dist_deriv_nuc fdi2(*this,m,T,m.n->n*(1.0-get_phi(m,T)),0);
  gd2.deriv_err(xx,fdi2,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[0]=der;
  cout << endl;
  
  // Return to original values
  m.n->n=xn;
  m.dist[0].n=xnuc;
  free_energy_dist(m,T);

  t.test_rel(mu1,m.mun,1.0e-8,"mu1");
  t.test_rel(mu2,m.mu[0],1.0e-10,"mu2");

  // Combine to form the gibbs energy density and pressure
  m.gb=m.n->n*m.mun;
  for(size_t j=0;j<m.dist.size();j++) {
    m.gb+=m.dist[j].n*m.mu[j];
  }
  m.pr=m.gb-m.fr;

  t.test_rel(m.pr,pr1,1.0e-7,"pr1");
  t.test_rel(m.gb,gb1,1.0e-10,"gb1");

  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  cout << "------------------------------------------------------" << endl;
  cout << "Test distribution with two nuclei:\n" << endl;

  xn=7.0e-3;
  xnuc=1.0e-5;
  xnuc2=2.1e-5;

  m.n->n=xn;

  m.dist[0].Z=32;
  m.dist[0].N=126;
  m.dist[0].n=xnuc;

  nucleus nuc2;
  nuc2.Z=50;
  nuc2.N=40;
  nuc2.n=xnuc2;
  m.dist.push_back(nuc2);

  gibbs_energy_dist(m,T);
  
  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  mu1=m.mun;
  mu2=m.mu[0];
  mu3=m.mu[1];
  pr1=m.pr;
  gb1=m.gb;
  
  xx=m.n->n*(1.0-get_phi(m,T));
  gd.h=1.0e-4;
  gd.deriv_err(xx,fda,der,dere);
  cout << "alpha,dfda,err: " << xx << " " << der << " " << dere << endl;
  m.mun=der*(1.0-get_phi(m,T));
  cout << "dfdn: " << m.mun << endl;

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;

  gd2.h=1.0e-8;

  free_dist_deriv_nuc fdi3(*this,m,T,m.n->n*(1.0-get_phi(m,T)),0);
  xx=m.dist[0].n;
  gd2.deriv_err(xx,fdi3,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[0]=der;

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;

  free_dist_deriv_nuc fdi4(*this,m,T,m.n->n*(1.0-get_phi(m,T)),1);
  xx=m.dist[1].n;
  gd2.deriv_err(xx,fdi4,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[1]=der;
  cout << endl;

  t.test_rel(mu1,m.mun,1.0e-10,"mu1");
  t.test_rel(mu2,m.mu[0],1.0e-10,"mu2");
  t.test_rel(mu3,m.mu[1],1.0e-10,"mu3");

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;
  free_energy_dist(m,T);

  // Combine to form the gibbs energy density and pressure
  m.gb=m.n->n*m.mun;
  for(size_t j=0;j<m.dist.size();j++) {
    m.gb+=m.dist[j].n*m.mu[j];
  }
  m.pr=m.gb-m.fr;

  t.test_rel(m.pr,pr1,1.0e-8,"pr1");
  t.test_rel(m.gb,gb1,1.0e-10,"gb1");

  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  cout << "------------------------------------------------------" << endl;
  cout << "Test distribution with three nuclei:\n" << endl;

  xn=7.0e-3;
  xnuc=1.0e-5;
  xnuc2=2.1e-5;
  xnuc3=1.0e-6;

  m.n->n=xn;

  m.dist[0].Z=32;
  m.dist[0].N=126;
  m.dist[0].n=xnuc;
  m.dist[1].Z=50;
  m.dist[1].N=40;
  m.dist[1].n=xnuc2;

  nucleus nuc3;
  nuc3.Z=40;
  nuc3.N=80;
  nuc3.n=xnuc3;
  m.dist.push_back(nuc3);
  
  gibbs_energy_dist(m,T);
  
  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  mu1=m.mun;
  mu2=m.mu[0];
  mu3=m.mu[1];
  mu4=m.mu[2];
  pr1=m.pr;
  gb1=m.gb;
  
  xx=m.n->n*(1.0-get_phi(m,T));
  gd.h=1.0e-4;
  gd.deriv_err(xx,fda,der,dere);
  cout << "alpha,dfda,err: " << xx << " " << der << " " << dere << endl;
  m.mun=der*(1.0-get_phi(m,T));
  cout << "dfdn: " << m.mun << endl;

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;
  m.dist[2].n=xnuc3;

  gd2.h=1.0e-8;

  free_dist_deriv_nuc fdi50(*this,m,T,m.n->n*(1.0-get_phi(m,T)),0);
  xx=m.dist[0].n;
  gd2.deriv_err(xx,fdi50,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[0]=der;

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;
  m.dist[2].n=xnuc3;

  free_dist_deriv_nuc fdi51(*this,m,T,m.n->n*(1.0-get_phi(m,T)),1);
  xx=m.dist[1].n;
  gd2.deriv_err(xx,fdi51,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[1]=der;

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;
  m.dist[2].n=xnuc3;

  free_dist_deriv_nuc fdi52(*this,m,T,m.n->n*(1.0-get_phi(m,T)),2);
  xx=m.dist[2].n;
  gd2.deriv_err(xx,fdi52,der,dere);
  cout << "n_nuc,dfdi,err: " << xx << " " << der << " " << dere << endl;
  m.mu[2]=der;
  cout << endl;

  t.test_rel(mu1,m.mun,1.0e-10,"mu1");
  t.test_rel(mu2,m.mu[0],1.0e-10,"mu2");
  t.test_rel(mu3,m.mu[1],1.0e-10,"mu3");
  t.test_rel(mu4,m.mu[2],1.0e-10,"mu4");

  m.n->n=xn;
  m.dist[0].n=xnuc;
  m.dist[1].n=xnuc2;
  m.dist[2].n=xnuc3;
  free_energy_dist(m,T);

  // Combine to form the gibbs energy density and pressure
  m.gb=m.n->n*m.mun;
  for(size_t j=0;j<m.dist.size();j++) {
    m.gb+=m.dist[j].n*m.mu[j];
  }
  m.pr=m.gb-m.fr;

  t.test_rel(m.pr,pr1,1.0e-9,"pr1");
  t.test_rel(m.gb,gb1,1.0e-10,"gb1");

  cout << "Neutrons (n,nu): " << m.n->n << " " << m.mun << endl;
  cout << "Nuclei: " << endl;
  for(size_t i=0;i<m.dist.size();i++) {
    cout << "\t" << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
  }
  cout << "gb,fr,pr: " << m.gb << " " << m.fr << " " << m.pr << endl;
  cout << endl;

  t.report();
  
  return 0;
}

int dist_thermo::gibbs_energy_dist(matter &m, double T) {
    
  if (!finite(m.n->n) || !finite(m.e.n)) {
    cout << "T: " << T << endl;
    cout << m << endl;
    O2SCL_ERR("Not finite in gibbs_energy_dist() [point 1].",
	      o2scl::exc_efailed);
  }

  // Ensure the matter arrays have the correct size
  if (m.dist.size()!=m.mu.size() ||
      m.dist.size()!=m.zeta.size()) {
    m.mu.resize(m.dist.size());
    m.zeta.resize(m.dist.size());
  }

  for(size_t i=0;i<m.dist.size();i++) {
    if (m.dist[i].n<0.0) {
      cout << "T: " << T << endl;
      cout << m << endl;
      O2SCL_ERR("Density negative in free_energy_dist().",o2scl::exc_efailed);
    }
  }

  // Compute the chemical potentials in free_energy_dist()
  free_energy_dist(m,T,true);
  
  // Combine to form the gibbs energy density
  m.gb=m.n->n*m.mun;
  for(size_t j=0;j<m.dist.size();j++) {
    m.gb+=m.dist[j].n*m.mu[j];
  }

  // Compute pressure
  m.pr=m.gb-m.fr;

  if (!finite(m.pr)) {
    cout << "Neutrons: " << m.n->n << " " << m.mun << endl;
    cout << m.n->pr << " " << m.n->ed << endl;
    cout << endl;
    for(size_t i=0;i<m.dist.size();i++) {
      cout << i << " " << m.dist[i].n << " " << m.mu[i] << endl;
    }
    cout << endl;
    cout << "gb,fr: " << m.gb << " " << m.fr << endl;
    O2SCL_ERR2("Pressure not finite in ",
	       "dist_thermo::gibbs_energy_dist().",o2scl::exc_efailed);
  }

  return 0;
}

int dist_thermo::free_energy_cell(matter &m, double T, int ix, double &fr, 
				  double &vol, double &Nprime) {
  
  double fr_neut=0.0, fr_elec;

  vector<nucleus> &dist=m.dist;

  if (m.n->n>0.0 || m.p->n>0.0) {
      
    // Ensure a consistent guess for the chemical potentials
    m.n->nu=m.n->m*1.0001;
    m.p->nu=m.p->m*1.0001;
      
    // Compute drip free energy (with rest mass contribution)
    thermo drip_th;
    het->calc_temp_e(*m.n,*m.p,T,drip_th);
    fr_neut=drip_th.ed-T*drip_th.en;
  }

  // Compute total proton density
  double nptot=m.p->n;
  for(size_t i=0;i<dist.size();i++) {
    nptot+=dist[i].n*dist[i].Z;
  }

  // Compute electron contribution
  m.e.mu=m.e.m+1.0e-10;
  m.e.n=nptot;
  relf.calc_density(m.e,T);
  fr_elec=m.e.ed-T*m.e.en;

  double Rws, chi;
  double be;
  
  be=lda->nucleus_be(dist[ix].Z,dist[ix].N,m.p->n,m.n->n,T,m.e.n,
		     Rws,chi);
  
  // Compute the thermodynamics of the nucleus, without the rest
  // mass which will be added in later
  if ((dist[ix].Z+dist[ix].N)%2==0) dist[ix].g=1.0;
  else dist[ix].g=3.0;
  double mass_neutron=cng.convert
    ("kg","1/fm",o2scl_mks::mass_neutron);
  double mass_proton=cng.convert
    ("kg","1/fm",o2scl_mks::mass_proton);
  dist[ix].m=dist[ix].Z*mass_proton+dist[ix].N*mass_neutron+be;
  dist[ix].non_interacting=true;
  dist[ix].inc_rest_mass=false;
    
  // Compute Nprime
  double Rn3=lda->Rn*lda->Rn*lda->Rn;
  double Rws3=Rws*Rws*Rws;
  vol=4.0/3.0*pi*Rws3;
  Nprime=dist[ix].N-4.0/3.0*pi*Rn3*m.n->n;
  
  fr=be+dist[ix].Z*mass_proton+Nprime*mass_neutron+vol*(fr_neut+fr_elec);
    
  return 0;
}

double dist_thermo::get_phi(matter &m, double T) {
    
  vector<nucleus> &dist=m.dist;
    
  // Compute total electron density
  double ne=0.0;
  for(size_t j=0;j<dist.size();j++) {
    ne+=dist[j].n*dist[j].Z;
  }

  double phi=0.0;
  for(size_t j=0;j<dist.size();j++) {
      
    double Rws, chi_i;
    //if (use_pasta) {
    //lda.nucleus_be_pasta(dist[j].Z,dist[j].N,m.p->n,m.n->n,T,ne,Rws,chi_i);
    //} else {
    lda->nucleus_be(dist[j].Z,dist[j].N,m.p->n,m.n->n,T,ne,Rws,chi_i);
    //}
    
    phi+=dist[j].n*4.0/3.0*pi*pow(lda->Rn,3.0);

  }

  return phi;
}
  
int dist_thermo::baryon_density(matter &m, double T) {
    
  vector<nucleus> &dist=m.dist;
    
  // Compute total proton density
  double nptot=0.0;
  for(size_t j=0;j<dist.size();j++) {
    nptot+=dist[j].n*dist[j].Z;
  }

  double nbtot=m.n->n+m.p->n;

  // Add nuclear contribution to the baryon density
  for(size_t i=0;i<dist.size();i++) {

    double Rws, chi, be;
    //if (use_pasta) {
    //be=lda.nucleus_be_pasta(dist[i].Z,dist[i].N,m.p->n,
    //m.n->n,T,m.e.n,Rws,chi);
    //} else {
    be=lda->nucleus_be(dist[i].Z,dist[i].N,m.p->n,m.n->n,T,m.e.n,
		       Rws,chi);
    //}
      
    // Compute Nprime
    double Nprime=dist[i].N-4.0/3.0*pi*lda->Rn*lda->Rn*lda->Rn*m.n->n;
      
    nbtot+=dist[i].n*(dist[i].Z+Nprime);
  }
    
  m.nb=nbtot;

  return 0;
}

int dist_thermo::gibbs_press_neut_func(double pr, double nnf, double alpha,
				       double nn_new, matter &m_new, double T,
				       double &y1, double &y2) {
  
  if (alpha<=0.0) return 1;
  if (nn_new<0.0) return 2;

  // Scale densities by alpha and set neutron density
  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n*=alpha;
  }
  m_new.n->n=nn_new;

  if (!finite(m_new.n->n) || !finite(m_new.e.n)) {
    cout << "alpha,nn_new,T: " 
	 << alpha << " " << nn_new << " " << T << endl;
    cout << m_new << endl;
    O2SCL_ERR("Not finite in gibbs_pressure_funct2().",o2scl::exc_efailed);
  }

  // Compute pressure and phi
  gibbs_energy_dist(m_new,T);
  double phi=get_phi(m_new,T);

  // This happens if one of the nuclei is unphysical and returned
  // a binding energy of 1.0e100
  if (m_new.pr==0.0) {
    return 3;
  }

  // Return the nuclear distribution densities back to the original value
  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n/=alpha;
    if (!finite(m_new.dist[i].n)) {
      cout << "alpha,nn_new,T: " 
	   << alpha << " " << nn_new << " " << T << endl;
      cout << m_new << endl;
      O2SCL_ERR2("Error in post density.",
		 "in gibbs_pressure_funct2().",o2scl::exc_efailed);
    }
  }
  
  y1=(m_new.pr-pr)/pr;
  y2=(m_new.n->n*(1.0-phi)-alpha*nnf)/(alpha*nnf);

  if (!finite(y1) || !finite(y2)) {
    cout << "nn_new,phi,alpha,nnf: " 
	 << nn_new << " " << phi << " " << alpha << " " << nnf << endl;
    O2SCL_ERR2("Variables y1 or y2 not finite in ",
	       "gibbs_press_neut_func().",o2scl::exc_efailed);
  }

  return 0;
}

int dist_thermo::gibbs_fixp_neutron(double pr_old, double nn_full_old,
				    matter &m_new, double T, double nn_shift) {

  if (nn_full_old==0.0 && nn_shift==0.0) {
    return gibbs_fixp(pr_old,m_new,T);
  }
  
  if (!finite(pr_old) || !finite(nn_full_old)) {
    cout << "nn_full, pr_old: " << nn_full_old << " " << pr_old << endl;
    O2SCL_ERR("Error in gibbs_fixp_neutron().",o2scl::exc_efailed);
  }
    
  funct_solve_pressure2 fsp2(*this,m_new,T,pr_old,nn_full_old);
  mm_funct11 fmp2=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&funct_solve_pressure2::operator()),&fsp2,
     std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  funct_solve_pressure3 fsp3(*this,m_new,T,pr_old,nn_full_old);
  mm_funct11 fmp3=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&funct_solve_pressure3::operator()),&fsp3,
     std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
    
  mroot_hybrids<> gmh;
  gmh.ntrial*=10;
  ubvector x(2), y(2);
  x[0]=1.0;
  x[1]=nn_full_old+nn_shift;

  // Try first, and double check that the pressure isn't zero
  // (which means that the new nucleus doesn't exist)
  int ret=fsp2(2,x,y);
  if (ret!=0) return -2;

  gmh.err_nonconv=false;
  int retx=gmh.msolve(2,x,fmp2);
  gmh.err_nonconv=true;

  // 1/7/17: I'm not sure if this is a good idea, but this
  // just bypasses the solver if the neutron density is very small
  if (retx!=0) {
    if (nn_full_old<1.0e-10) {
      return -3;
    } else {
      O2SCL_ERR("Solver failed in gibbs_fixp_neutron().",o2scl::exc_efailed);
    }
  }

  /*
  if (retx!=0 && nn_full_old<1.0e-8) {
    cout << nn_full_old+nn_shift << " " << x[1] << endl;
    x[0]=1.0;
    if (x[1]<0.0) x[1]=log(1.0e-13);
    x[1]=log(nn_full_old+nn_shift);
    gmh.err_nonconv=false;
    int retx2=gmh.msolve(2,x,fmp3);
    x[1]=exp(x[1]);
    cout << "retx2: " << retx2 << " " << x[1] << endl;
    cout << x[1] << endl;
    exit(-1);
  }
  */

  if (x[0]<0.0) {
    cout << "alpha, nn: " << x[0] << " " << x[1] << endl;
    O2SCL_ERR2("Variable alpha is negative in",
	       " gibbs_fixp_neutron().",o2scl::exc_efailed);
  }

  // Ensure m_new has the correct new composition
  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n*=x[0];
  }
  m_new.n->n=x[1];


  return 0;
}

double dist_thermo::gibbs_pressure_func
(double pr0, double n_fact, matter &m_new, double T) {
    
  if (n_fact<0.0 || !finite(n_fact)) {
    cout << n_fact << endl;
    O2SCL_ERR("n_fact negative in gibbs_pressure_func().",
	      o2scl::exc_efailed);
  }

  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n*=n_fact;
  }
  m_new.n->n*=n_fact;

  gibbs_energy_dist(m_new,T);
  double pr=m_new.pr;
  
  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n/=n_fact;
  }
  m_new.n->n/=n_fact;

  return (pr-pr0)/pr0;
}

int dist_thermo::gibbs_fixp(double pr_target, matter &m_new, double T) {
    
  funct_solve_pressure fsp(*this,m_new,T,pr_target);

  root_cern<funct_solve_pressure> cmr;
  root_brent_gsl<funct_solve_pressure> grb;

  double ll=0.5;
  double ul=2.0;

  // The variable 'corr' is the ratio of old to new densities
  double corr=ll;

  // 1/7/17: Just try the CERN solver first, then fall back to the
  // GSL one if the CERN solver fails
  cmr.err_nonconv=false;
  corr=1.0;
  int cret=cmr.solve(corr,fsp);
  
  if (cret!=0) {
    for(size_t i=0;i<3;i++) {
      if (fsp(ll)*fsp(ul)>0.0) {
	ll=sqrt(ll);
	ul=sqrt(ll);
      }
    }
    if (fsp(ll)*fsp(ul)>0.0) {
      for(corr=ll;corr<=ul;corr*=1.1) {
	cout << corr << " " << fsp(corr) << endl;
      }
      cout << ul << " " << fsp(ul) << endl;
      O2SCL_ERR("Couldn't bracket root in dist_thermo::gibbs_fixp().",
		o2scl::exc_efailed);
    }
    corr=ll;
    grb.solve_bkt(corr,ul,fsp);
    if (corr<ll || corr>ul) {
      O2SCL_ERR("Out of range in gibbs_fixp().",o2scl::exc_efailed);
    }
  }

  // Evaluate the gibbs energy at the chosen value of corr
  for(size_t i=0;i<m_new.dist.size();i++) {
    m_new.dist[i].n*=corr;
  }
  m_new.n->n*=corr;
  gibbs_energy_dist(m_new,T);

  return 0;
}

double dist_thermo::dist_sna_solve(double nb, matter &m, double T) {
  baryon_density(m,T);
  double nbx=m.nb;
  return (nbx-nb)/nb;
}

int dist_thermo::free_energy_dist_sna(double nb, matter &m, double T) {

  // Using the solver fails at about rho = 2e8 g/cm^3 currently.
  if (false) {
      
    root_cern<dist_thermo::funct_dist_sna> cmr;
    dist_thermo::funct_dist_sna fda(*this,m,T,nb);
    cmr.solve(m.dist[0].n,fda);
      
  } else {

    baryon_density(m,T);
    double nbx=m.nb;
    for(size_t i=0;i<10;i++) {
      m.dist[0].n+=(nb-nbx)/(m.dist[0].Z+m.dist[0].N);
      baryon_density(m,T);
      nbx=m.nb;
    }

  }

  return 0;
}

int dist_thermo::gibbs_energy_per_baryon_cell
(matter &m, double T, double P, int ix, double &gpb) {
  double fr, vol, Nprime;
  free_energy_cell(m,T,ix,fr,vol,Nprime);
  gpb=(fr+vol*P)/(m.dist[ix].Z+Nprime+vol*m.n->n);
  return 0;
}

int dist_thermo::gibbs_energy_cell
(matter &m, double T, double P, int ix, double &g_cell) {
  double fr, vol, Nprime;
  free_energy_cell(m,T,ix,fr,vol,Nprime);
  g_cell=fr+vol*P;
  return 0;
}

