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
#include <o2scl/test_mgr.h>
#include <o2scl/nstar_cold.h>

#include "crust.h"

using namespace std;
using namespace o2scl;
using namespace crust;
using namespace o2scl_const;
using namespace o2scl_hdf;

bool crust::compare_density(const nucleus &n1, const nucleus &n2) {
  return (n1.n>n2.n);
}

bool crust::compare_Z(const nucleus &n1, const nucleus &n2) {
  if (n1.Z==n2.Z) return n1.N<n2.N;
  return (n1.Z<n2.Z);
}

crust_driver::crust_driver() : 
  cng(o2scl_settings.get_convert_units()),
  dt(cng), nmt(lda,sk,cng), snat(lda,sk) {
  
  model_set=false;

  rho_summary=1.0e20;

  use_pasta=false;
  simple_pyc=false;
  more_reactions=true;
    
  n_lda.init(cng.convert("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  p_lda.init(cng.convert("kg","1/fm",o2scl_mks::mass_proton),2.0);
	     
  n_lda.inc_rest_mass=true;
  p_lda.inc_rest_mass=true;
  n_lda.non_interacting=false;
  p_lda.non_interacting=false;

  het=&sk;

  lda.full_surface=true;
  lda.exc_volume=true;
  lda.new_skin_mode=false;
  lda.set_eos_had_temp_base(*het);
  lda.set_n_and_p(n_lda,p_lda);
  lda.hd_exp=5.0;
  lda.hd_coeff=0.5;
    
  Tptr=cng.convert("K","1/fm",1.0e8);
  
  // Ensure the SNA minimizer is sufficiently accurate
  cm.tol_abs=1.0e-12;

  ubvector x(lda.nfit);

  // New version with Moller et al. data for SLy4
  x[0]=9.208247e-01;
  x[1]=1.164440e+00;
  x[2]=1.964531e+00;
  x[3]=9.041767e-01;
  x[4]=-1.577306e-02;
  x[5]=1.740151e-01;
  x[6]=5.277076e+00;
  x[7]=-1.217548e+00;
  x[8]=2.566948e-02;
  x[9]=3.870849e-03;
  x[10]=3.575317e-02;

  lda.fit_fun(lda.nfit,x);

  check=check_none;
  calc_heating=true;
  ec_heating=0.25;
  verbose=1;
  dist_type="nickel";
  allow_palpha=false;
  acc_inc_factor=1.01;

  cf.lda=&lda;
  dt.lda=&lda;
  dt.het=het;

  feq_nb=2.0e-10;
  feq_Z=28;
  feq_N=34;
  feq_nn=1.0e-10;
  feq_end=0.09;
  feq_dnb=1.1;
}

crust_driver::~crust_driver() {
}

int crust_driver::init_dist(string mode, matter &m) {

  m.n->n=0.0;
  m.p->n=0.0;

  vector<nucleus> &current=m.dist;

  size_t ninit;
  if (mode=="nickel") {

    if (verbose>0) cout << "Nickel-56" << endl;
    ninit=1;
    int Z[1]={26};
    int N[1]={30};
    double frac[1]={1.0};
      
    for(size_t i=0;i<ninit;i++) {
      nucleus nuc;
      nuc.Z=Z[i];
      nuc.N=N[i];
      nuc.n=frac[i]*1.1e-11;
      current.push_back(nuc);
    }

  } else if (mode=="heavy") {

    if (verbose>0) cout << "Palladium-56" << endl;
    ninit=1;
    int Z[1]={46};
    int N[1]={60};
    double frac[1]={1.0};
      
    for(size_t i=0;i<ninit;i++) {
      nucleus nuc;
      nuc.Z=Z[i];
      nuc.N=N[i];
      nuc.n=frac[i]*5.5e-12;
      current.push_back(nuc);
    }

  } else if (mode=="schatz") {

    if (verbose>0) cout << "Ashes from Schatz et al." << endl;
    string fn="data/schatz_edit.dat";
    ifstream fin(fn.c_str());
    int sZ, sA;
    double abun;

    while (fin >> sZ) {
      fin >> sA >> abun;

      nucleus nuc;
      nuc.Z=sZ;
      nuc.N=sA-sZ;
      nuc.n=abun*6.0e-10;
      current.push_back(nuc);
    }

    fin.close();

  } else if (mode=="ashes") {

    // From Table 1 in Horowitz (2007)
    if (verbose>0) cout << "X-ray burst ashes." << endl;
    ninit=17;
    int Z[17], N[17];
    double den[17];
    Z[0]=8;
    Z[1]=10;
    Z[2]=12;
    Z[3]=14;
    Z[4]=15;
    Z[5]=20;
    Z[6]=22;
    Z[7]=24;
    Z[8]=26;
    Z[9]=27;
    Z[10]=28;
    Z[11]=30;
    Z[12]=32;
    Z[13]=33;
    Z[14]=34;
    Z[15]=36;
    Z[16]=47;
    den[0]=0.0301;
    den[1]=0.0116;
    den[2]=0.0023;
    den[3]=0.0023;
    den[4]=0.0023;
    den[5]=0.0046;
    den[6]=0.0810;
    den[7]=0.0718;
    den[8]=0.1019;
    den[9]=0.0023;
    den[10]=0.0764;
    den[11]=0.0856;
    den[12]=0.0116;
    den[13]=0.1250;
    den[14]=0.3866;
    den[15]=0.0023;
    den[16]=0.0023;
    for(size_t i=0;i<17;i++) {
      int Nbest=0;
      double Emin=1.0e12;
      for(int N0=Z[i]/2-1;N0<Z[i]*2+1;N0++) {
	if (cf.ame.binding_energy(Z[i],N0)/(N0+Z[i])<Emin) {
	  // Find the minimum binding energy per nucleon
	  Emin=cf.ame.binding_energy(Z[i],N0)/(N0+Z[i]);
	  Nbest=N0;
	}
      }
      N[i]=Nbest;
    }

    for(size_t i=0;i<ninit;i++) {
      nucleus nuc;
      nuc.Z=Z[i];
      nuc.N=N[i];
      nuc.n=den[i]*9.0e-12;
      current.push_back(nuc);
    }
  } else {
    O2SCL_ERR("Unknown initial distribution type.",o2scl::exc_efailed);
  }

  // Begin with no dripped neutrons or protons
  m.n->n=0.0;
  m.p->n=0.0;

  // Output baryon density
  dt.baryon_density(m,Tptr);
  if (verbose>0) cout << "Initial baryon density: " 
		      << m.nb << " fm^{-3}" << endl;

  // Output mass density
  dt.mass_density(m,Tptr);
  if (verbose>0) {
    cout << "Initial mass density: " << m.rho << " g/cm^3" << endl;
  }
    
  // Ensure there's enough space in the matter arrays
  m.mu.resize(current.size());
  m.zeta.resize(current.size());

  return 0;
}
  
int crust_driver::delta_ZN(int &Z, int &N, int &deltaZ_lo, int &deltaZ_hi, 
			   int &deltaN_lo, int &deltaN_hi) {

  deltaZ_lo=8;
  deltaZ_hi=8;
  deltaN_lo=16;
  deltaN_hi=16;

  // Ensure protons take us to next shell closure
  if (Z>=28) {
    deltaZ_lo=50-28;
    deltaZ_hi=50-28;
  }
  if (Z>=50) {
    deltaZ_lo=82-50;
    deltaZ_hi=82-50;
  }

  // Ensure neutrons take us to next shell closure, or sometimes
  // even higher neutron numbers
  if (N>=50) {
    deltaN_lo=82-50;
    deltaN_hi=82-50;
  }
  if (N>=82) {
    deltaN_lo=126-82;
    deltaN_hi=308-82;
  }
  if (N>=126) {
    deltaN_lo=184-126;
    deltaN_hi=406-126;
  }

  return 0;
}

int crust_driver::check_free_energy_cell_fun() {

  cout << "Checking free_energy_cell()." << endl;
  
  cout.precision(8);

  // Distribution
  cout << "Initializing for distribution." << endl;
  matter m;
  init_dist("ashes",m);

  double T=Tptr;
    
  vector<nucleus> &dist=m.dist;
  double tot=0.0;
  for(size_t i=0;i<dist.size();i++) {
    double tmp,tmp2,tmp3;
    dt.free_energy_cell(m,T,((int)i),tmp,tmp2,tmp3);
    tot+=tmp*dist[i].n;
    cout.width(2);
    cout << i << " "<< tmp*dist[i].n << endl;
  }
  cout << endl;
  cout << "These numbers should be equal: " << endl;
  cout << tot << endl;
  dt.free_energy_dist(m,T,false);
  cout << m.fr << endl;
  cout << endl;
  cout << fabs(m.fr-tot)/fabs(m.fr) << endl;

  test_mgr t;
  t.set_output_level(2);
  t.test_rel(m.fr,tot,1.0e-10,"free_energy_call().");
  t.report();

  return 0;
}

int crust_driver::compute_sna(double nb, double T, matter &m, bool debug) {
    
  if (m.dist.size()!=1) {
    O2SCL_ERR("Wrong distribution size in compute_sna().",
	      o2scl::exc_esanity);
  }

  if (feq_fix_mode) {
    m.dist[0].Z=fix_Z;
    m.dist[0].N=fix_N;
  }

  // Initial guess for Z and N
  int Z_guess=m.dist[0].Z;
  int N_guess=m.dist[0].N;

  // Function to minimize
  sna_thermo::free_energy_sna_neut func(snat,m,T,nb);

  // Minimizer
  min_cern<sna_thermo::free_energy_sna_neut> cm2;
  cm2.tol_abs=1.0e-12;

  // Set nn and full_min to the minimum with
  // the previous nucleus
  double full_min=1.0e100;
  cm2.min_bkt(m.n->n,0.0,nb,full_min,func);

  // Storage for the current best minimum
  int Zmin=Z_guess, Nmin=N_guess;
  double nn_min=m.n->n;
    
  // The initial guess for the neutron density at each nucleus
  double nn_guess=m.n->n;
    
  // The current minimum
  double min=0.0;
    
  // Compute loop sizes
  int deltaZ_lo, deltaZ_hi, deltaN_lo, deltaN_hi;
  delta_ZN(Z_guess,N_guess,deltaZ_lo,deltaZ_hi,deltaN_lo,deltaN_hi);

  if (!feq_fix_mode) {

    for(int tZ=Z_guess-deltaZ_lo;tZ<=Z_guess+deltaZ_hi;tZ++) {
      for(int tN=N_guess-deltaN_lo;tN<=N_guess+deltaN_hi;tN++) {
	
	if (tZ<10) tZ=10;
	if (tN<4) tN=4;
	
	// Minimize to determine neutron density
	m.n->n=nn_guess;
	m.dist[0].Z=tZ;
	m.dist[0].N=tN;

	cm2.min_bkt(m.n->n,0.0,nb,min,func);

	/*
	  if (debug) {
	  dt.baryon_density(m,T);
	  cout << "Z,N,nn,n_nuc,fr: " 
	  << tZ << " " << tN << " " << m.n->n << " " 
	  << m.dist[0].n << " " << min << endl;
	  }
	*/
	if (min<full_min) {
	  Zmin=tZ;
	  Nmin=tN;
	  nn_min=m.n->n;
	  full_min=min;
	}

	if (verbose>1) {
	  cout << "crust_driver::compute_sna(): Z,N,nn,n_nuc,fr,Zmin,Nmin:\n"
	       << "  " << tZ << " " << tN << " " << m.n->n << " " 
	       << m.dist[0].n << " " << min << " " << Zmin << " "
	       << Nmin << " " << nn_min << " " << full_min << endl;
	}
	
      }
    }

  } else {

    Zmin=fix_Z;
    Nmin=fix_N;
    nn_min=m.n->n;

    cm2.min_bkt(m.n->n,0.0,nb,min,func);
    full_min=min;
  }

  // Set user variables
  m.dist[0].Z=Zmin;
  m.dist[0].N=Nmin;
  m.n->n=nn_min;

  // Double-check the minimization
  cm2.min_bkt(m.n->n,0.0,nb,m.fr,func);
    
  // Compute the mass density
  dt.mass_density(m,T);

  // Compute number density of nuclei
  m.dist[0].n=(nb-m.n->n)/
    (m.dist[0].Z+m.dist[0].N-4.0/3.0*pi*pow(lda.Rn,3.0)*m.n->n);
  
  // Compute pressure
  dt.gibbs_energy_dist(m,T);

  // Reset free energy, etc. to original values 
  snat.free_energy_sna_fix_nb_nn(nb,m,T);

  if (true) {

    // Compute c_e^2
    sna_thermo::free_press_sna_fixed_ye gesy_fr(snat,m,T,dt,true);
    sna_thermo::free_press_sna_fixed_ye gesy_pr(snat,m,T,dt,false);
    
    o2scl::deriv_gsl<sna_thermo::free_press_sna_fixed_ye> gd;
    gd.h=nb/1.0e4;
    m.dPdnb_Ye=gd.deriv(nb,gesy_pr);
    m.dfdnb_Ye=gd.deriv(nb,gesy_fr);
    
  }

  // This checks to ensure the minimum in m.n->n*(1.0-phi) is
  // the same as the minimum in m.n->n. One can test this with,
  // e.g. acc -feq dir_x 0.07 100 40 0.065 0.09. This
  // appears to work somewhat, but is limited by the accuracy
  // of free_energy_sna_fix_nb_nn() and other functions.
  if (false) {
    cout.precision(8);

    double cent=m.n->n;
    for(m.n->n=cent/(1.0+1.0e-3);m.n->n<cent*(1.0+1.0e-3);
	m.n->n*=(1.0+1.0e-4)) {
      snat.free_energy_sna_fix_nb_nn(nb,m,T);
      cout << m.n->n << " " << dt.get_phi(m,T) << " "
	   << m.n->n*(1.0-dt.get_phi(m,T)) << " " << m.fr << endl;
    }
    // Return to original values
    m.n->n=cent;
    snat.free_energy_sna_fix_nb_nn(nb,m,T);
    cout << endl;

    cout << "nN,nn,fr: " 
	 << m.dist[0].n << " " << m.n->n << " " << full_min << endl;
    double alpha=m.n->n*(1.0-dt.get_phi(m,T));
    cout << "alpha,phi,nb: " << alpha << " " << dt.get_phi(m,T) << " "
	 << nb << endl;
    double alphaL=alpha/(1.0+1.0e-3);
    //double alphaL=alpha;
    double alphaR=alpha*(1.0+2.0e-3);
    
    for(alpha=alphaL;alpha<alphaR;alpha*=(1.0+1.0e-4)) {
      sna_thermo::free_energy_sna_fix_nb_nnhat fmna(snat,nb,m,T);
      double ret=fmna(alpha);
      cout << m.n->n << " " << alpha << " " << ret << endl;
    }
    min_cern<sna_thermo::free_energy_sna_fix_nb_nnhat> cm3;
    double min2=0.0;
    //cm3.verbose=2;
    sna_thermo::free_energy_sna_fix_nb_nnhat fmna2(snat,nb,m,T);
    cm3.min_bkt(alpha,alphaL,alphaR,min2,fmna2);
    cout << alpha << " " << m.dist[0].n << " " << m.n->n << " " 
	 << min2 << endl;
    exit(-1);
  }

  if (debug) {
    cout << endl;
    cout << "Final: " 
	 << m.dist[0].Z << " " << m.dist[0].N << " " << m.n->n << " " 
	 << m.dist[0].n << " " << m.fr << endl;
  }
  
  return 0;
}

int crust_driver::compute_sna_dist(double nb, matter &m, double T, 
				   bool debug) {
    
  if (m.dist.size()!=1) {
    O2SCL_ERR("Wrong distribution size in compute_sna_dist().",
	      o2scl::exc_esanity);
  }

  // Initial guess for Z and N, and the number density of nuclei
  int Z_guess=m.dist[0].Z;
  int N_guess=m.dist[0].N;
  m.dist[0].n=nb/(Z_guess+N_guess)/2.0;

  // Set nn and full_min to the minimum with
  // the previous nucleus
  double full_min=1.0e100;
  dist_thermo::free_energy_dist_neut func(dt,m,T,nb);
  min_cern<dist_thermo::free_energy_dist_neut> cm2;
  cm2.min_bkt(m.n->n,0.0,nb,full_min,func);

  // Storage for the current best minimum
  int Zmin=Z_guess, Nmin=N_guess;
  double nn_min=m.n->n;
    
  // The initial guess for the neutron density at each nucleus
  double nn_guess=m.n->n;
    
  // The current minimum
  double min=0.0;
    
  if (debug) cout.precision(8);
    
  // Compute loop sizes
  int deltaZ_lo, deltaZ_hi, deltaN_lo, deltaN_hi;
  delta_ZN(Z_guess,N_guess,deltaZ_lo,deltaZ_hi,deltaN_lo,deltaN_hi);

  for(int tZ=Z_guess-deltaZ_lo;tZ<=Z_guess+deltaZ_hi;tZ++) {
    for(int tN=N_guess-deltaN_lo;tN<=N_guess+deltaN_hi;tN++) {
	
      if (tZ<10) tZ=10;
      if (tN<4) tN=4;

      // Minimize to determine neutron density
      m.n->n=nn_guess;
      m.dist[0].Z=tZ;
      m.dist[0].N=tN;
      cm2.min_bkt(m.n->n,0.0,nb,min,func);

      if (debug) {
	dt.baryon_density(m,T);
	cout << tZ << " " << tN << " " << m.n->n << " " << min << " " 
	     << m.dist[0].n << " " << m.nb << endl;
      }

      if (min<full_min) {
	Zmin=tZ;
	Nmin=tN;
	nn_min=m.n->n;
	full_min=min;
      }
    }
  }
    
  // Set user variables
  m.dist[0].Z=Zmin;
  m.dist[0].N=Nmin;
  m.n->n=nn_min;

  // Double-check the minimization
  cm2.min_bkt(m.n->n,0.0,nb,m.fr,func);
    
  // Compute the mass density
  dt.mass_density(m,T);

  // Reset free energy, etc. to original values 
  dt.free_energy_dist(m,T,false);

  if (debug) {
    cout << endl;
    cout << m.dist[0].Z << " " << m.dist[0].N << " " << m.n->n << " " 
	 << m.fr << " " << m.rho << endl;
    exit(-1);
  }

  return 0;
}

bool crust_driver::pycno_allowed(double Z1, double Z2, double rho) {
  if (Z1*Z2<16.1) return true;
  else if (rho>1.5e12 && Z1*Z2<100.1) return true;
  else if (rho>1.0e13 && Z1*Z2<144.1) return true;
  //else if (rho>2.0e13 && Z1*Z2<200.0) return true;
  else if (rho>2.0e13) return true;
  return false;
}

int crust_driver::prune_distribution(vector<nucleus> &dist) {
    
  // Compute total number density
  typedef vector<nucleus>::iterator iter_t;
  double ntot=0.0;
  for(iter_t it=dist.begin();it!=dist.end();it++) {
    ntot+=it->n;
  }
    
  // Remove nuclei with negligable density
  bool clear_done=false;
  while (clear_done==false) {
    clear_done=true;
    for(iter_t it=dist.begin();clear_done==true && 
	  it!=dist.end();it++) {
      if (it->n==0.0 || it->n<ntot*1.0e-10) {
	//cout << "Pruning: " << it->Z << " " << it->N << " " << it->n << endl;
	dist.erase(it);
	clear_done=false;
      }
    }
  }

  return 0;
}

int crust_driver::update_flow(int Z1, int N1, int Z2, int N2,
			      std::string type, double dheat, double rho) {
  if (type==((string)"ec")) heat_full+=dheat/ec_heating;
  else heat_full+=dheat;
  ofstream fout(flow_fn.c_str(),ios::app);
  fout.setf(ios::scientific);
  fout << Z1 << " " << N1 << " " << Z2 << " " << N2 << " " 
       << type << " " << rho << " " << dheat << endl;
  fout.close();
  return 0;
}
  
bool crust_driver::dist_switch_gb(matter &m, matter &m_new, 
				  int &cnt, bool &restart, double &heat, 
				  double T, string type, int debug) {
  
  if (calc_heating) {
    // Compute the heating as the difference in the gibbs
    // energy per baryon.
    heat+=(m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm;
  }
    
  cnt++;
    
  // Copy over distribution
  if (true) {
    m.dist=m_new.dist;
  } else {
    m.dist.clear();
    nucleus nuc;
    for(size_t i=0;i<m_new.dist.size();i++) {
      m.dist.push_back(nuc);
      m.dist[i].Z=m_new.dist[i].Z;
      m.dist[i].N=m_new.dist[i].N;
      m.dist[i].n=m_new.dist[i].n;
    }
    if (m.dist.size()!=m_new.dist.size()) {
      cout << "Problem in dist_switch_gb()." << endl;
      exit(-1);
    }
  }

  // Copy over neutron density
  m.n->n=m_new.n->n;
  
  // Make sure proton density is set
  // (I'm not sure if this is necessary)
  m.p->n=0.0;

  m.mu.resize(m.dist.size());
  m.zeta.resize(m.dist.size());
    
  restart=true;
  
  prune_distribution(m.dist);
    
  return true;
}

int crust_driver::make_table(std::vector<std::string> &sv, bool itive_com) {

  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  cout << "Making table." << endl;
  vector<double> den_grid;
  thermo tab_th;
  den_grid.push_back(0.0);
  den_grid.push_back(1.0e-6);
  den_grid.push_back(1.0e-4);
  den_grid.push_back(2.0e-4);
  den_grid.push_back(5.0e-4);
  den_grid.push_back(1.0e-3);
  den_grid.push_back(1.5e-3);
  den_grid.push_back(2.3e-3);
  den_grid.push_back(3.5e-3);
  den_grid.push_back(5.3e-3);
  den_grid.push_back(8.0e-3);
  for(double nbt=0.01;nbt<=0.24001;nbt+=0.001) {
    den_grid.push_back(nbt);
  }
  den_grid.push_back(0.36);
  den_grid.push_back(0.50);
  den_grid.push_back(0.75);
  den_grid.push_back(0.83);
  den_grid.push_back(0.91);
  den_grid.push_back(0.94);
  den_grid.push_back(0.97);
  size_t n=den_grid.size();
  t3d.tab.set_xy("nn",n,den_grid,"np",n,den_grid);
  t3d.tab.new_slice("ed");
  t3d.tab.new_slice("pr");
  t3d.tab.new_slice("mun");
  t3d.tab.new_slice("mup");

  fermion ne_mt(cng.convert("kg","1/fm",o2scl_mks::mass_neutron),2.0);
		
  fermion pmt(cng.convert("kg","1/fm",o2scl_mks::mass_proton),2.0);
	      
  ne_mt.inc_rest_mass=true;
  pmt.inc_rest_mass=true;
  ne_mt.non_interacting=false;
  pmt.non_interacting=false;

  for(size_t i=0;i<n;i++) {
    for(size_t j=0;j<n;j++) {
      ne_mt.n=den_grid[i];
      pmt.n=den_grid[j];
      het->calc_e(ne_mt,pmt,tab_th);
      t3d.tab.set_val(ne_mt.n,pmt.n,"ed",tab_th.ed);
      t3d.tab.set_val(ne_mt.n,pmt.n,"pr",tab_th.pr);
      t3d.tab.set_val(ne_mt.n,pmt.n,"mun",ne_mt.mu);
      t3d.tab.set_val(ne_mt.n,pmt.n,"mup",pmt.mu);
    }
  }
  t3d.tab.set_interp_type(itp_linear);
  cout << "Done. Setting EOS pointer to new table EOS." << endl;

  ne_mt.n=2.0e-3;
  pmt.n=0.165;
  het->calc_e(ne_mt,pmt,tab_th);
  cout << "Saturation (before table): " 
       << (tab_th.ed-0.08*ne_mt.m-0.08*pmt.m)/0.16*hc_mev_fm << " ";
  cout << ne_mt.mu << " " << pmt.mu << " " << tab_th.pr << endl;

  het=&t3d;
  lda.set_eos_had_temp_base(*het);
  
  ne_mt.n=2.0e-3;
  pmt.n=0.165;
  het->calc_e(ne_mt,pmt,tab_th);
  cout << "Saturation (after table) : " 
       << (tab_th.ed-0.08*ne_mt.m-0.08*pmt.m)/0.16*hc_mev_fm << " ";
  cout << ne_mt.mu << " " << pmt.mu << " " << tab_th.pr << endl;

  return 0;
}

int crust_driver::model(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cout << "No model specified." << endl;
    return o2scl::exc_efailed;
  }

  if (sv[1]=="APR") {

    cout << "Selected APR model." << endl;
    het=&apr;
    dt.het=&apr;
    snat.het=&apr;
    nmt.het=&apr;

  } else if (sv[1]=="RMF") {
    
    if (sv.size()<3) {
      cout << "RMF model specified but no data file." << endl;
      return o2scl::exc_efailed;
    }

    cout << "Loading model from " << sv[2] << " ." << endl;
    hdf_file hf;
    hf.open(sv[2]);
    
    rmf.mnuc=939.0/hc_mev_fm;
    rmf.mw=782.501/hc_mev_fm;
    rmf.mr=763.0/hc_mev_fm;

    hf.getd("cs",rmf.cs);
    hf.getd("cw",rmf.cw);
    hf.getd("cr",rmf.cr);
    hf.getd("b",rmf.b);
    hf.getd("c",rmf.c);
    hf.getd("zeta",rmf.zeta);
    hf.getd("xi",rmf.xi);
    hf.getd("a1",rmf.a1);
    hf.getd("a2",rmf.a2);
    hf.getd("a3",rmf.a3);
    hf.getd("a4",rmf.a4);
    hf.getd("a5",rmf.a5);
    hf.getd("a6",rmf.a6);
    hf.getd("b1",rmf.b1);
    hf.getd("b2",rmf.b2);
    hf.getd("b3",rmf.b3);
    hf.getd("ms",rmf.ms);

    het=&rmf;
    dt.het=&rmf;
    snat.het=&rmf;
    nmt.het=&rmf;

    hf.close();

  } else if (sv[1]=="skyrme") {

    cout << "Selecting Skyrme with specified coefficients." << endl;
    skyrme_load(sk,"NRAPR");
    sk.t0=o2scl::function_to_double(sv[2])/hc_mev_fm;
    sk.t1=o2scl::function_to_double(sv[3])/hc_mev_fm;
    sk.t2=o2scl::function_to_double(sv[4])/hc_mev_fm;
    sk.t3=o2scl::function_to_double(sv[5])/hc_mev_fm;
    sk.x0=o2scl::function_to_double(sv[6]);
    sk.x1=o2scl::function_to_double(sv[7]);
    sk.x2=o2scl::function_to_double(sv[8]);
    sk.x3=o2scl::function_to_double(sv[9]);
    sk.alpha=o2scl::function_to_double(sv[10]);
    
  } else {

    cout << "Selected Skyrme model." << endl;
    skyrme_load(sk,sv[1]);
    model_set=true;
    het=&sk;
    dt.het=&sk;

  }

  lda.set_eos_had_temp_base(*het);

  return 0;
}

int crust_driver::full_eq(std::vector<std::string> &sv, bool itive_com) {

  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  if (tsh.mag_field>0.0) dt.mag_field=tsh.mag_field;

  if (sv.size()<2) {
    cout << "No filename prefix given." << endl;
    return o2scl::exc_efailed;
  }
  string prefix=sv[1];

  if (feq_fix_mode) {
    feq_Z=fix_Z;
    feq_N=fix_N;
  }

  if (sv.size()>2) feq_nb=o2scl::stod(sv[2]);
  if (sv.size()>3) feq_Z=o2scl::stoi(sv[3]);
  if (sv.size()>4) feq_N=o2scl::stoi(sv[4]);
  if (sv.size()>5) feq_nn=o2scl::stod(sv[5]);
  if (sv.size()>6) feq_end=o2scl::stod(sv[6]);
    
  cout.precision(3);

  double T=0.0;
    
  matter m;
    
  // Set initial values
  nucleus nuc;
  nuc.Z=feq_Z;
  nuc.N=feq_N;
  m.dist.push_back(nuc);
  m.mu.resize(1);
  m.n->n=feq_nn;
  m.p->n=0.0;

  table_units<> feq;
  feq.line_of_names("nb rho N Z nn bulk surf coul pair shell exc");
  feq.line_of_names("BEoA ede ne fr fr_check nnuc chi mnuc fr_x nnL npL");
  feq.line_of_names("kfn nun nuN pr Rws Rn Rp Ncell pr_x nb_x dPdnb_Ye");
  feq.line_of_names("dfdnb_Ye");

  feq.set_unit("nb","1/fm^3");
  feq.set_unit("rho","g/cm^3");
  feq.set_unit("nn","1/fm^3");
  feq.set_unit("bulk","MeV");
  feq.set_unit("surf","MeV");
  feq.set_unit("coul","MeV");
  feq.set_unit("pair","MeV");
  feq.set_unit("shell","MeV");
  feq.set_unit("exc","MeV");
  feq.set_unit("BEoA","MeV");
  feq.set_unit("ede","1/fm^4");
  feq.set_unit("ne","1/fm^3");
  feq.set_unit("fr","1/fm^4");
  feq.set_unit("fr_check","1/fm^4");
  feq.set_unit("nnuc","1/fm^3");
  feq.set_unit("mnuc","1/fm");
  feq.set_unit("fr_x","1/fm^4");
  feq.set_unit("nnL","1/fm^3");
  feq.set_unit("npL","1/fm^3");
  feq.set_unit("kfn","1/fm");
  feq.set_unit("nun","1/fm^3");
  feq.set_unit("nuN","1/fm^3");
  feq.set_unit("pr","1/fm^4");
  feq.set_unit("Rws","fm");
  feq.set_unit("Rn","fm");
  feq.set_unit("Rp","fm");
  feq.set_unit("pr_x","1/fm^4");
  feq.set_unit("nb_x","1/fm^3");
  feq.set_unit("dPdnb_Ye","1/fm");
  feq.set_unit("dfdnb_Ye","1/fm");

  int it=0;

  // For output file
  hdf_file hf;
  string fname=prefix+"_feq.o2";

  for(double nbx=feq_nb;nbx<=feq_end;nbx*=feq_dnb,it++) {

    compute_sna(nbx,T,m);
      
    // This is probably unnecessary since it's computed in 
    // compute_sna().
    dt.mass_density(m,T);
      
    if (it%20==0) {
      cout.setf(ios::left);
      cout.width(10); cout << "nb";
      cout.width(10); cout << "rho";
      cout.unsetf(ios::left);
      cout.width(4); cout << "Z ";
      cout.width(4); cout << "N ";
      cout.setf(ios::left);
      cout.width(10); cout << "nn";
      cout.width(11); cout << " BE/A";
      cout.width(10); cout << "chi";
      cout.width(10); cout << "fr";
      cout.width(10); cout << "pr";
      cout.width(10); cout << "fr_x";
      cout << "  ";
      cout.width(10); cout << "nn";
      cout.width(10); cout << "np";
      cout.unsetf(ios::left);
      cout << endl;
    }

    cout << nbx << " " << m.rho << " "; 
    cout.width(3);
    cout << m.dist[0].Z << " ";
    cout.width(3);
    cout << m.dist[0].N << " " << m.n->n << " ";
      
    // Binding energy, without excluded volume correction in fm^{-1}
    double Rws, chi, be, Rneut, Rprot;
    if (use_pasta) {
      be=lda.nucleus_be_pasta(m.dist[0].Z,m.dist[0].N,0.0,
			      m.n->n,T,m.e.n,Rws,chi)-lda.exc/hc_mev_fm;
    } else {
      be=lda.nucleus_be(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,T,m.e.n,Rws,chi)-
	lda.exc/hc_mev_fm;
    }
    Rneut=lda.Rn;
    Rprot=lda.Rp;
    cout << be*hc_mev_fm/(m.dist[0].N+m.dist[0].Z) << " ";
      
    // Total mass in fm^{-1}.
    double mass_neutron=cng.convert("kg","1/fm",o2scl_mks::mass_neutron);
      
    double mass_proton=cng.convert("kg","1/fm",o2scl_mks::mass_proton);
      
    double mnuc=be+m.dist[0].Z*mass_proton+m.dist[0].N*mass_neutron;
      
    cout << chi << " ";

    double fr_neutron;
    {
      // Compute neutron drip free energy
      thermo drip_th;

      // Ensure a consistent guess for the chemical potentials
      m.n->nu=m.n->m*1.0001;
      m.p->nu=m.p->m*1.0001;

      if (T<=0.0) {
	het->calc_e(*m.n,*m.p,drip_th);
	fr_neutron=drip_th.ed-m.n->n*m.n->m;
      } else {
	het->calc_temp_e(*m.n,*m.p,T,drip_th);
	fr_neutron=drip_th.ed-T*drip_th.en-m.n->n*m.n->m;
      }
    }

    cout << m.fr << " ";

    double fr_check=mnuc*m.dist[0].n+(1.0-chi)*
      (m.n->n*mass_neutron+fr_neutron)+m.e.ed;
    if (false) {
      if (fabs(fr_check-m.fr)/fabs(m.fr)>1.0e-4) {
	O2SCL_ERR2("Free energy check failed in ",
		   "crust_driver::full_eq().",o2scl::exc_efailed);
      }
    }
    cout << m.pr << " ";
      
    matter nm;
    if (m.rho>1.0e9) {
      m.T=T;
      nmt.calc_nm_from_dist(m,nm);
    } else {
      nm.fr=0.0;
      nm.pr=0.0;
      nm.n->n=0.0;
      nm.p->n=0.0;
    }
      
    cout << nm.fr << " " << (m.fr>nm.fr) << " " << lda.nn << " "
	 << lda.np << endl;

    static const size_t ncols=34;
    if (feq.get_ncolumns()!=ncols) {
      O2SCL_ERR("Column mismatch.",o2scl::exc_efailed);
    }
    
    double A=(m.dist[0].N+m.dist[0].Z);
    double line[ncols]={nbx,m.rho,((double)m.dist[0].N),
			((double)m.dist[0].Z),
			m.n->n,lda.bulk/A,lda.surf/A,lda.coul/A,
			lda.pair/A,lda.shell/A,lda.exc/A,be*hc_mev_fm/A,
			m.e.ed,m.e.n,m.fr,fr_check,m.dist[0].n,chi,
			mnuc,nm.fr,lda.nn,lda.np,m.n->kf,m.mun,m.mu[0],
			m.pr,Rws,Rneut,Rprot,m.dist[0].N+
			4.0/3.0*pi*m.n->n*(pow(Rws,3.0)-pow(Rneut,3.0)),
			nm.pr,nm.n->n+nm.p->n,m.dPdnb_Ye,m.dfdnb_Ye};
    feq.line_of_data(ncols,line);
    
    hf.open_or_create(fname);
    hdf_output(hf,feq,"feq");
    hf.close();
  }

  if (feq.get_nlines()<=2) {

    cout << "Not enough lines in table for pressure." << endl;

  } else {

    table_units<> *t2=(table_units<> *)(&feq);

    //def_interp_mgr<ubvector_const_view,linear_interp> lint1;
    //def_interp_mgr<ubvector_const_subvector,linear_interp> lint2;
    //t2->set_interp(lint1,lint2);

    t2->function_column("fr/nb","EoA");
    t2->set_unit("EoA","1/fm");
    t2->deriv("nb","EoA","dEoAdnb");
    t2->set_unit("dEoAdnb","fm^2");
    t2->function_column("nb*nb*dEoAdnb","pr_int");
    t2->set_unit("pr_int","1/fm^4");

    o2scl::nstar_cold cns;
    cns.nb_start=0.02;
    cns.include_muons=false;

    fermion nx(cng.convert("kg","1/fm",o2scl_mks::mass_neutron),2.0);
	       
    fermion px(cng.convert("kg","1/fm",o2scl_mks::mass_proton),2.0);
	       
    nx.non_interacting=false;
    px.non_interacting=false;
    //cns.set_n_and_p(nx,px);
    cns.set_eos(*het);
    cns.calc_eos();
    shared_ptr<table_units<> > te=cns.get_eos_results();

    t2->new_column("ed_nm");
    t2->new_column("pr_nm");
    t2->set_unit("ed_nm","1/fm^4");
    t2->set_unit("pr_nm","1/fm^4");

    for(size_t i=0;i<t2->get_nlines();i++) {
      double nbx=t2->get("nb",i);
      if (nbx<0.02) {
	t2->set("ed_nm",i,0.0);
	t2->set("pr_nm",i,0.0);
      } else {
	t2->set("ed_nm",i,te->interp("nb",nbx,"ed"));
	t2->set("pr_nm",i,te->interp("nb",nbx,"pr"));
      }
    }
  }

  hf.open_or_create(fname);
  hdf_output(hf,feq,"feq");
  hf.close();
  
  return 0;
}

int crust_driver::full_eq2(std::vector<std::string> &sv, bool itive_com) {

  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  if (sv.size()<2) {
    cout << "No filename prefix given." << endl;
    return o2scl::exc_efailed;
  }
  string prefix=sv[1];

  if (sv.size()>2) feq_nb=o2scl::stod(sv[2]);
  if (sv.size()>3) feq_Z=o2scl::stoi(sv[3]);
  if (sv.size()>4) feq_N=o2scl::stoi(sv[4]);
  if (sv.size()>5) feq_nn=o2scl::stod(sv[5]);
  if (sv.size()>6) feq_end=o2scl::stod(sv[6]);
    
  cout.precision(3);

  double T=0.0;
    
  matter m;
    
  // Set initial values
  nucleus nuc;
  nuc.Z=feq_Z;
  nuc.N=feq_N;
  m.dist.push_back(nuc);
  m.mu.resize(1);
  m.dist[0].n=feq_nb;
  m.n->n=feq_nn;
  m.p->n=0.0;

  table_units<> feq;
  feq.line_of_names("nb rho N Z nn bulk surf coul pair shell exc");
  feq.line_of_names("BEoA ede ne fr fr_check nnuc chi mnuc fr_nm nnL npL ");
  feq.line_of_names("kfn nun nuN pr2");

  feq.set_unit("nb","1/fm^3");
  feq.set_unit("rho","g/cm^3");
  feq.set_unit("nn","1/fm^3");
  feq.set_unit("bulk","MeV");
  feq.set_unit("surf","MeV");
  feq.set_unit("coul","MeV");
  feq.set_unit("pair","MeV");
  feq.set_unit("shell","MeV");
  feq.set_unit("exc","MeV");
  feq.set_unit("BEoA","MeV");
  feq.set_unit("ede","1/fm^4");
  feq.set_unit("ne","1/fm^3");
  feq.set_unit("fr","1/fm^4");
  feq.set_unit("fr_check","1/fm^4");
  feq.set_unit("nnuc","1/fm^3");
  feq.set_unit("mnuc","1/fm");
  feq.set_unit("fr_nm","1/fm^4");
  feq.set_unit("nnL","1/fm^3");
  feq.set_unit("npL","1/fm^3");
  feq.set_unit("kfn","1/fm");
  feq.set_unit("nun","1/fm^3");
  feq.set_unit("nuN","1/fm^3");
  feq.set_unit("pr2","1/fm^4");
    
  int it=0;
  for(double nbx=feq_nb;nbx<=feq_end;nbx*=feq_dnb,it++) {
      
    compute_sna_dist(nbx,m,T);

    m.pr=0.0;
      
    // This is probably unnecessary since it's computed in 
    // compute_sna().
    dt.mass_density(m,T);
      
    if (it%20==0) {
      cout << "nb rho Z N nn be/A chi fr fr_nm phase pr nn np" << endl;
    }
      
    cout << nbx << " " << m.rho << " "; 
    cout.width(3);
    cout << m.dist[0].Z << " ";
    cout.width(3);
    cout << m.dist[0].N << " " << m.n->n << " ";
      
    double Rws, chi, be;
    if (use_pasta) {
      be=lda.nucleus_be_pasta(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,T,m.e.n,
			      Rws,chi)-lda.exc/hc_mev_fm;
    } else {
      be=lda.nucleus_be(m.dist[0].Z,m.dist[0].N,0.0,m.n->n,T,m.e.n,
			Rws,chi)-lda.exc/hc_mev_fm;
    }
    cout << be*hc_mev_fm/(m.dist[0].N+m.dist[0].Z) << " ";
      
    double mass_neutron=cng.convert("kg","1/fm",o2scl_mks::mass_neutron);
      
    double mass_proton=cng.convert("kg","1/fm",o2scl_mks::mass_proton);
      
    double mnuc=be+m.dist[0].Z*mass_proton+m.dist[0].N*mass_neutron;
      
    cout << chi << " ";

    double fr_neutron;
    {
      // Compute neutron drip free energy
      thermo drip_th;

      // Ensure a consistent guess for the chemical potentials
      m.n->nu=m.n->m*1.0001;
      m.p->nu=m.p->m*1.0001;

      het->calc_temp_e(*m.n,*m.p,T,drip_th);
      fr_neutron=drip_th.ed-T*drip_th.en-m.n->n*m.n->m;
    }
	
    cout << m.fr << " ";
    double fr_check=mnuc*m.dist[0].n+(1.0-chi)*
      (m.n->n*mass_neutron+fr_neutron)+m.e.ed;

    matter nm;
    m.T=T;
    nmt.calc_nm_from_dist(m,nm);
      
    cout << nm.fr << " " << (m.fr>nm.fr) << " ";

    // Compute pressure
    dt.gibbs_energy_dist(m,T);
    cout << m.pr << " ";
      
    // Output neutron and proton densities
    cout << lda.nn << " "
	 << lda.np << endl;
      
    static const size_t ncols=26;
    if (feq.get_ncolumns()!=ncols) {
      O2SCL_ERR("Column mismatch.",o2scl::exc_efailed);
    }

    double A=(m.dist[0].N+m.dist[0].Z);
    double line[ncols]={nbx,m.rho,((double)m.dist[0].N),
			((double)m.dist[0].Z),m.n->n,lda.bulk/A,
			lda.surf/A,lda.coul/A,lda.pair/A,lda.shell/A,
			lda.exc/A,be*hc_mev_fm/A,m.e.ed,m.e.n,m.fr,
			fr_check,m.dist[0].n,chi,mnuc,nm.fr,
			lda.nn,lda.np,m.n->kf,m.mun,m.mu[0],m.pr};
    feq.line_of_data(ncols,line);
  }

  if (feq.get_nlines()<=1) {
    cout << "Not enough lines in table." << endl;
    return 0;
  }

  table_units<> *t2=(table_units<> *)(&feq);
  t2->function_column("fr/nb","EoA");
  t2->set_unit("EoA","1/fm");
  t2->deriv("nb","EoA","dEoAdnb");
  t2->set_unit("dEoAdnb","fm^2");
  t2->function_column("nb*nb*dEoAdnb","pr");
  t2->set_unit("pr","1/fm^4");

  nstar_cold cns;
  cns.nb_start=0.02;
  cns.include_muons=false;
  
  fermion nx(cng.convert("kg","1/fm",o2scl_mks::mass_neutron),2.0);
	     
  fermion px(cng.convert("kg","1/fm",o2scl_mks::mass_proton),2.0);
	     
  nx.non_interacting=false;
  px.non_interacting=false;
  //cns.set_n_and_p(nx,px);
  cns.set_eos(*het);
  cns.calc_eos();
  shared_ptr<table_units<> > te=cns.get_eos_results();
  
  t2->new_column("ed_nm");
  t2->new_column("pr_nm");
  t2->set_unit("ed_nm","1/fm^4");
  t2->set_unit("pr_nm","1/fm^4");
  
  for(size_t i=0;i<t2->get_nlines();i++) {
    double nbx=t2->get("nb",i);
    if (nbx<0.02) {
      t2->set("ed_nm",i,0.0);
      t2->set("pr_nm",i,0.0);
    } else {
      t2->set("ed_nm",i,te->interp("nb",nbx,"ed"));
      t2->set("pr_nm",i,te->interp("nb",nbx,"pr"));
    }
  }
  
  hdf_file hf;
  string fname=prefix+"_feq.o2";
  hf.open(fname);
  hdf_output(hf,feq,"feq");
  hf.close();
  
  return 0;
}

int crust_driver::full_eq_dist(std::vector<std::string> &sv, bool itive_com) {
    
  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  if (sv.size()<2) {
    cout << "No filename prefix given." << endl;
    return o2scl::exc_efailed;
  }
  string prefix=sv[1];

  size_t i_start=0;
  if (sv.size()==3) {
    i_start=o2scl::stoi(sv[2]);
  }
    
  double T=Tptr;
  //double T=0.0;
    
  matter m, m_new;
  m.p->n=0.0;
  m.n->n=0.0;
  m_new.p->n=0.0;
  m_new.n->n=0.0;

  // Load single nucleus table
  table_units<> feq;
  hdf_file hf;
  string fname=prefix+"_feq.o2";
  hf.open(fname);
  hdf_input(hf,feq,"feq");
  hf.close();

  // Output table
  table_units<> feqd;
  feqd.line_of_names("nb rho nn ede ne fr a A Z");

  feqd.set_unit("nb","1/fm^3");
  feqd.set_unit("rho","g/cm^3");
  feqd.set_unit("nn","1/fm^3");
  feqd.set_unit("ede","1/fm^4");
  feqd.set_unit("ne","1/fm^3");
  feqd.set_unit("fr","1/fm^4");
  feqd.set_unit("a","fm");

  // The nuclear distributions
  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;
    
  int cnt;

  for(size_t i=i_start;i<feq.get_nlines();i++) {

    double heat;

    double nb=feq.get("nb",i);

    // Start with fresh distribution
    current.clear();
      
    nucleus nuc;
    nuc.Z=((int)(feq.get("Z",i)+1.0e-10));
    nuc.N=((int)(feq.get("N",i)+1.0e-10));
    nuc.n=feq.get("nnuc",i);
    current.push_back(nuc);
    m.n->n=feq.get("nn",i);
    m.mu.resize(1);

    cnt=0;
    {
      rn.delta_n=0.1;
      rn.gen_reaction(this,m,m_new,T,cnt);
      rn.emit_neutron(this,m,m_new,T,cnt,heat);
      rn.neut_capture(this,m,m_new,T,cnt,heat);
    }
    {
      rn.delta_n=0.01;
      rn.gen_reaction(this,m,m_new,T,cnt);
      rn.emit_neutron(this,m,m_new,T,cnt,heat);
      rn.neut_capture(this,m,m_new,T,cnt,heat);
    }
      
    dt.mass_density(m,T);
    nb=dt.baryon_density(m,T);

    cout.width(3);
    cout << i << " ";
      
    cout << nb << " " << feq.get("nb",i) << " " << m.rho << " "
	 << feq.get("Z",i) << " " << feq.get("N",i) << " ";
    cout << m.fr << " " << feq.get("fr",i) << " ";

    // Total baryon density (without neutrons)
    double nb_nuclei=0.0;
    for(size_t j=0;j<current.size();j++) {
      nb_nuclei+=(current[j].Z+current[j].N)*current[j].n;
    }

    cout.width(2);
    cout << current.size() << " ";
      
    for(size_t j=0;j<current.size() && j<2;j++) {
      cout << "(" << current[j].Z << "," << current[j].N << ") ";
      cout << current[j].n*(current[j].Z+current[j].N)/nb_nuclei << " ";
    }
      
    cout << cnt << endl;

    double line[9]={nb,m.rho,m.n->n,m.e.ed,m.e.en,m.fr,
		    m.average_a(),m.average_A(),m.average_Z()};
    feqd.line_of_data(9,line);

  }

  cout << "Going to save." << endl;
  fname=prefix+"_feqd.o2";
  hf.open(fname);
  hdf_output(hf,feqd,"feqd");
  hf.close();
    
  return 0;
}

int crust_driver::test_ndrip(std::vector<std::string> &sv, bool itive_com) {

  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  matter m;

  cout << "T: " << Tptr << endl;
  cout << endl;

  //double n_nuc=2.0e-6, shift=2.0e-8;
  //nucleus nuc;
  //int tZ=36;
  //int tN=84;
  double n_nuc=5.0e-6, shift=5.0e-8;
  nucleus nuc;
  int tZ=16;
  int tN=40;

  nuc.Z=tZ;
  nuc.N=tN;
  m.dist.push_back(nuc);
  m.dist[0].Z=tZ;
  m.dist[0].N=tN;
  m.dist[0].n=n_nuc;

  m.n->n=0.0;

  cout.precision(8);

  dt.gibbs_energy_dist(m,Tptr);
  cout << "No neutrons: " << endl;
  cout << "nb,rho,fr,pr: " << m.nb << " " << m.rho << " " 
       << m.fr << " " << m.pr << endl;
  cout << "gb: " << m.gb << endl;
  cout << m.fr << " " << dt.part1 << " " << dt.part2 << " " 
       << dt.part3 << " " << dt.part4 << endl;
  double t1=m.fr;
  double t2=dt.part1;
  double t3=dt.part2;
  double t4=dt.part3;
  double t5=dt.part4;
  cout << "nn: " << m.n->n << endl;
  cout << "nN: " << m.dist[0].n << endl;
  cout << "nb2: " 
       << m.dist[0].n*(m.dist[0].Z+m.dist[0].N)+
    m.n->n << endl;
  cout << "gpb: " << m.gb/m.nb << endl;
  cout << endl;
  
  // -----------------------------------------------

  matter m_new;

  nuc.Z=tZ;
  nuc.N=tN;
  m_new.dist.push_back(nuc);
  m_new.dist[0].Z=tZ;
  m_new.dist[0].N=tN;
  m_new.dist[0].n=n_nuc-shift;

  nuc.Z=tZ;
  nuc.N=tN-1;
  m_new.dist.push_back(nuc);
  m_new.dist[1].Z=tZ;
  m_new.dist[1].N=tN-1;
  m_new.dist[1].n=shift;

  m_new.n->n=shift;
  
  cout << "Add neutrons, and compare gb/nb at fixed pr: " << endl;
  dt.gibbs_fixp_neutron(m.pr,shift,m_new,Tptr,shift);
  cout << "nb,rho,fr,pr: " << m_new.nb << " " << m_new.rho << " " 
       << m_new.fr << " " << m_new.pr << endl;
  cout << "gb: " << m_new.gb << endl;
  cout << m_new.fr << " " << dt.part1 << " " << dt.part2 << " " 
       << dt.part3 << " " << dt.part4 << endl;
  cout << m_new.fr-t1 << " " << dt.part1-t2 << " " 
       << dt.part2-t3 << " " << dt.part3-t4 << " "
       << dt.part4-t5 << endl;
  cout << "nN: " << m_new.dist[0].n << " " << m_new.dist[1].n << endl;
  cout << "nn: " << m_new.n->n << endl;
  cout << "nb2: " 
       << m_new.dist[0].n*(m_new.dist[0].Z+m_new.dist[0].N)+
    m_new.dist[1].n*(m_new.dist[1].Z+m_new.dist[1].N)+
    m_new.n->n << endl;
  cout << "gpb: " << m_new.gb/m_new.nb << endl;
  cout << endl;

  cout << m.gb/m.nb-m_new.gb/m_new.nb << endl;

  return 0;
}

int crust_driver::acc(std::vector<std::string> &sv, bool itive_com) {

  if (model_set==false) {
    O2SCL_ERR("Model not set.",o2scl::exc_efailed);
  }

  // Matter objects
  matter m, m_new, nm;

  // Store reaction summaries
  int Zmax=100;
  int Nmax=300;
  ubmatrix ec_gpb(Zmax,Nmax), en_gpb(Zmax,Nmax);
  ubmatrix ec_m(Zmax,Nmax), en_m(Zmax,Nmax);
  for(int i=0;i<Zmax;i++) {
    for(int j=0;j<Nmax;j++) {
      ec_gpb(i,j)=0.0;
      en_gpb(i,j)=0.0;
      ec_m(i,j)=0.0;
      en_m(i,j)=0.0;
    }
  }

  bool wrote_summary=false;

  // The current nuclear distribution
  vector<nucleus> &current=m.dist;
  // The trial nuclear distribution
  vector<nucleus> &trial=m_new.dist;
    
  double T=Tptr;
    
  if (sv.size()<2) {
    cout << "No filename prefix given." << endl;
    return o2scl::exc_efailed;
  }
  string prefix=sv[1];

  // Setup flow filename
  if (flow_fn.length()<1) {
    flow_fn=prefix+"_flow.txt";
  }

  cout.precision(3);

  // Output
  table_units<> dist_tab;

  // Heating
  double heat=0.0;
  heat_full=0.0;

  if (sv.size()>3) {
    cout << "File loading doesn't work right now "
	 << "because of first_col." << endl;
    exit(-1);

    string file=sv[2];
    double density=o2scl::stod(sv[3]);
    cout << "Loading file: " << file 
	 << " and retrieving composition at density rho="
	 << density << endl;
    hdf_file hf_in;
    hf_in.open(file);
    string name;
    hdf_input(hf_in,dist_tab,name);
    hf_in.close();
      
    size_t row=dist_tab.lookup("rho",density);
    cout << "Found density " << dist_tab["rho"][row] << " at row "
	 << row << endl;

    // Set neutron and proton densities. Initialize electron density
    m.n->n=dist_tab["nn"][row];
    m.p->n=0.0;
    m.e.n=0.0;
      
    for(size_t i=13;i<dist_tab.get_ncolumns();i++) {
      if (dist_tab.get(i,row)>0.0) {
	int tZ, tN, tA;
	nmi.parse_elstring(dist_tab.get_column_name(i),tZ,tN,tA);

	nucleus nuc;
	nuc.Z=tZ;
	nuc.N=tN;
	nuc.n=dist_tab.get(i,row);
	m.dist.push_back(nuc);
	  
	cout << "Added nucleus (" << tZ << "," << tN 
	     << ") to distribution with density " << nuc.n << endl;
	  
	m.e.n+=nuc.n*tZ;
      }
    }
      
    cout << "n,p,e: " << m.n->n << " " << m.p->n << " " << m.e.n << endl;
      
    // Ensure there's enough space in the matter arrays
    m.mu.resize(current.size());
    m.zeta.resize(current.size());
      
    heat=dist_tab["heat"][row];
      
    cout << "Removing end rows of table." << endl;
    while (dist_tab.get_nlines()>row) {
      dist_tab.delete_row(row);
    }

    dt.gibbs_energy_dist(m,T);
    cout << "Initial pressure: " << m.pr << endl;

  } else {
      
    // Initialize
    init_dist(dist_type,m);

    string names=((string)"nb rho Ye Yn fr Z A Q fr_x a heat ")+
      "heat_full nn pr gb ne mun mue phi t_acc t_C t_O t_Ne t_Mg Nd";
    dist_tab.line_of_names(names);
    dist_tab.set_unit("nb","1/fm^3");
    dist_tab.set_unit("rho","g/cm^3");
    dist_tab.set_unit("fr","1/fm^4");
    dist_tab.set_unit("pr","1/fm^4");
    dist_tab.set_unit("gb","1/fm^4");
    dist_tab.set_unit("fr_x","1/fm^4");
    dist_tab.set_unit("a","fm");
    dist_tab.set_unit("heat","MeV");
    dist_tab.set_unit("heat_full","MeV");
    dist_tab.set_unit("nn","1/fm^3");
    dist_tab.set_unit("ne","1/fm^3");
    dist_tab.set_unit("mun","1/fm");
    dist_tab.set_unit("mue","1/fm");
    dist_tab.set_unit("t_acc","s");
    dist_tab.set_unit("t_C","s");
    dist_tab.set_unit("t_O","s");
    dist_tab.set_unit("t_Ne","s");
    dist_tab.set_unit("t_Mg","s");
  }

  if (lda.exc_volume) {
    cout << "Excluded volume correction." << endl;
  } else {
    cout << "No excluded volume correction." << endl;
  }
  if (lda.inc_shell) {
    cout << "Shell effects." << endl;
  } else {
    cout << "No shell effects." << endl;
  }
  if (simple_pyc) {
    cout << "Simple pycnonuclear reactions." << endl;
  } else {
    cout << "Full pycnonuclear reactions." << endl;
  }
  if (calc_heating) {
    cout << "Calculate heating." << endl;
  } else {
    cout << "No heating." << endl;
  }

  // Count number of reaction blocks
  int cnt=0;
    
  double rho=0.0;
  double max_rho=1.6e14;
  bool first_column=false;

  double pr_last=0.0;

  dist_tab.add_constant("worst_dec",0.0);

  double mub_last=0.0;

  if (check==check_acc_numbers) {
    cout << "Checking acc numbers." << endl;
    for(size_t i=0;i<current.size();i++) {
      current[i].n*=1.0e6;
    }
    max_rho=1.0e12;
  }
  
  // Increase the density slowly
  for(size_t i=0;rho<max_rho;i++) {

    // Sort (This is important for pynonuclear reactions,
    // since low Z fusions happen first)
    std::sort(current.begin(),current.end(),compare_Z);
      
    // Perform non-fusion reactions first
    cnt=0;
    bool done=false;
    int iix=0;
    while (!done) {
      int cnt2=0;
      rn.emit_neutron(this,m,m_new,T,cnt2,heat);
      rn.neut_capture(this,m,m_new,T,cnt2,heat);
      rn.elec_capture(this,m,m_new,T,cnt2,heat);
      rn.beta_decay(this,m,m_new,T,cnt2,heat);
      cnt+=cnt2;
      if (iix==50) {
	O2SCL_ERR("Non-fusion reactions never stopping in acc().",
		  o2scl::exc_efailed);
      }
      iix++;
      if (cnt2==0) done=true;
    }

    // Add fusions
    done=false;
    iix=0;
    while (!done) {
      int cnt2=0;
      rn.pyc_fusion(this,m,m_new,T,cnt2,heat);
      if (more_reactions) {
	rn.elec_capture(this,m,m_new,T,cnt2,heat);
	rn.beta_decay(this,m,m_new,T,cnt2,heat);
	rn.emit_neutron(this,m,m_new,T,cnt2,heat);
	rn.neut_capture(this,m,m_new,T,cnt2,heat);
      }
      cnt+=cnt2;
      if (iix>50) {
	for(size_t j=0;j<m.dist.size();j++) {
	  cout << "(" << m.dist[j].Z << "," << m.dist[j].N << ") ";
	  cout << m.dist[j].n << " ";
	}
	cout << cnt << endl;
      }
      if (iix==500) {
	O2SCL_ERR("Fusions never stopping in acc().",
		  o2scl::exc_efailed);
      }
      iix++;
      if (cnt2==0) done=true;
    }

    if (rho>1.0e13 && allow_palpha) {
	
      done=false;
      iix=0;
      while (!done) {
	int cnt2=0;
	rn.emit_fragment(this,m,m_new,T,cnt2,heat);
	cnt+=cnt2;
	if (iix==50) {
	  O2SCL_ERR("Alpha emission never stopping in acc().",
		    o2scl::exc_efailed);
	}
	iix++;
	if (cnt2==0) done=true;
      }
	
    }

    if (i%20==0) {

      cout.setf(ios::left);
      cout.width(10); cout << "nb";
      cout.width(10); cout << "rho";
      cout.width(10); cout << "<A>";
      cout.width(10); cout << "heat";
      cout.width(11); cout << " delta f";
      cout.width(10); cout << "fr";
      cout.width(11); cout << " pr";
      cout.width(10); cout << "nn";
      //cout.width(10); cout << "Y_n";
      cout.width(10); cout << "<Z>";
      cout.width(10); cout << "Q";
      cout << " Nd Distribution" << endl;
      cout.unsetf(ios::left);
    }

    // Compute current electron density
    double netmp=0.0;
    for(size_t j=0;j<current.size();j++) {
      netmp+=current[j].n*current[j].Z;
    }

    cout << m.nb << " ";
      
    dt.mass_density(m,T);
    rho=m.rho;
    cout << rho << " ";

    double avgA=m.average_A();
    cout << avgA << " ";
      
    double avga=m.average_a();
    //cout << avga << " ";
      
    cout << heat << " ";
    //cout << m.e.pr << " ";
      
    m.T=T;
    nmt.calc_nm_from_dist(m,nm);

    dt.free_energy_dist(m,T,false);
    cout.setf(ios::showpos);
    cout << nm.fr-m.fr << " ";
    cout.unsetf(ios::showpos);
    cout << m.fr << " ";
 
    dt.gibbs_energy_dist(m,T);
    cout.setf(ios::showpos);
    cout << m.pr << " ";
    cout.unsetf(ios::showpos);

    if (i!=0) {
      if (m.pr<pr_last && 
	  fabs(m.pr-pr_last)>dist_tab.get_constant("worst_dec")) {
	dist_tab.set_constant("worst_dec",fabs(m.pr-pr_last));
      }
    }
    pr_last=m.pr;
      
    cout << m.n->n << " ";
    // cout << m.n->n/m.nb << " ";

    double aveZ=m.average_Z();
    cout << aveZ << " ";

    double Q=m.impurity();
    cout << Q << " ";
      
    cout.width(3);
    cout << current.size() << " ";

    {
      // Sort with decreasing density so we can output only
      // the most populous nuclei
      vector<nucleus> sorted=current;
      std::sort(sorted.begin(),sorted.end(),compare_density);
      
      double nbnuc=0.0;
      for(size_t j=0;j<sorted.size();j++) {
	nbnuc+=sorted[j].n*(sorted[j].Z+sorted[j].N);
      }
      for(size_t j=0;j<sorted.size() && j<2;j++) {
	cout << "(" << sorted[j].Z << "," << sorted[j].N << ") ";
	cout << sorted[j].n*(sorted[j].Z+sorted[j].N)/(nbnuc) << " ";
      }
    }
    cout << cnt << endl;

    if (i==0) {
      cout << "Setting initial heating (" << heat << " and " << heat_full 
	   << ") to zero." << endl;
      heat=0.0;
      heat_full=0.0;
    }

    // Add a row to the table and zero it out, store the row index
    // in 'row'
    double line[1]={0.0};
    dist_tab.line_of_data(1,line);
    size_t row=dist_tab.get_nlines()-1;
    for(size_t j=0;j<dist_tab.get_ncolumns();j++) {
      dist_tab.set(j,row,0.0);
    }

    // For each nucleus in the distribution, set the density
    // in the current table row
    for(size_t j=0;j<current.size();j++) {
	
      string col;
      if (current[j].Z>119) {
	col=itos(current[j].Z)+"-"+itos(current[j].N);
      } else {
	col=nmi.tostring(current[j].Z,current[j].N);
      }
	  
      // Test if we should add a column to the table
      if (!dist_tab.is_column(col)) {
	dist_tab.new_column(col);
	if (first_column==false) {
	  size_t ix;
	  ix=dist_tab.lookup_column(col);
	  dist_tab.add_constant("first_col",((double)ix));
	  first_column=true;
	}
	dist_tab.set_unit(col,"1/fm^3");
	// Zero all data for the new column
	for(size_t k=0;k<dist_tab.get_nlines();k++) {
	  dist_tab.set(col,k,0.0);
	}
      }
	
      //dist_tab[col][row]=current[j].n;
      dist_tab.set(col,row,current[j].n);
    }
	
    dist_tab.set("nb",row,m.nb);
    dist_tab.set("rho",row,rho);
    dist_tab.set("Ye",row,m.e.n/m.nb);
    dist_tab.set("Yn",row,m.n->n/m.nb);
    dist_tab.set("fr",row,m.fr);
    dist_tab.set("pr",row,m.pr);
    dist_tab.set("gb",row,m.gb);
    dist_tab.set("fr_x",row,nm.fr);
    dist_tab.set("Z",row,aveZ);
    dist_tab.set("A",row,avgA);
    dist_tab.set("a",row,avga);
    dist_tab.set("heat",row,heat);
    dist_tab.set("heat_full",row,heat_full);
    dist_tab.set("nn",row,m.n->n);
    dist_tab.set("ne",row,m.e.n);
    dist_tab.set("Q",row,Q);
    if (m.n->n>0.0) {
      dist_tab.set("mun",row,m.n->mu);
    } else {
      dist_tab.set("mun",row,0.0);
    }
    dist_tab.set("mue",row,m.e.mu);
    dist_tab.set("phi",row,dt.get_phi(m,T));
    dist_tab.set("Nd",row,((double)current.size()));
    {
      double t_acc, t_c_fusion, t_o_fusion, t_ne_fusion, t_mg_fusion;
      pyc.fusion_times(dt,m,T,t_acc,t_c_fusion,t_o_fusion,
		       t_ne_fusion,t_mg_fusion);
      dist_tab.set("t_acc",row,t_acc);
      dist_tab.set("t_C",row,t_c_fusion);
      dist_tab.set("t_O",row,t_o_fusion);
      dist_tab.set("t_Ne",row,t_ne_fusion);
      dist_tab.set("t_Mg",row,t_mg_fusion);
    }

    // Now that the table has been updated, check the multizone
    // code
    if (false && i!=0) {
      //if (m.gb/m.nb<mub_last && m.n->n>1.0e-4) {
      if (rho>2.7e12) {
	mz.test_stability(T,dt,9,dist_tab,i-8);
	exit(-1);
      }
    }
    mub_last=m.gb/m.nb;

    // Write to output file periodically
    if (rho>max_rho*0.9 || rho>1.0e12 || i%10==0) {
      string fname2=prefix+"_acc.o2";
      hdf_file hf;
      hf.open_or_create(fname2);
      hdf_output(hf,dist_tab,"acc");
      hf.close();
    }

    // Write flow to file
    if (i==0) {
      // Truncate the file before appending
      ofstream fout;
      fout.open(flow_fn.c_str());
      fout.setf(ios::scientific);
      fout << "iter: " << i << " " << rho << endl;
      fout.close();
    } else {
      ofstream fout;
      fout.open(flow_fn.c_str(),ios::app);
      fout.setf(ios::scientific);
      fout << "iter: " << i << " " << rho << endl;
      fout.close();
    }
    
    // Add mass density to the grid and 
    // summarize reactions 
    if (wrote_summary==false && m.rho>rho_summary) {
      cout << "Writing ec_summary: " << m.rho << " " << rho_summary << endl;
      rn.ec_summary(this,m,m_new,T,ec_gpb,ec_m);
      rn.en_summary(this,m,m_new,T,en_gpb,en_m);
      string fnamex=prefix+"_rxns.o2";
      hdf_file hf;
      hf.open(fnamex);
      hf.setd_mat_copy("ec_gpb",ec_gpb);
      hf.setd_mat_copy("ec_m",ec_m);
      hf.setd_mat_copy("en_gpb",en_gpb);
      hf.setd_mat_copy("en_m",en_m);
      hf.close();
      wrote_summary=true;
    }

    // Increase the pressure for the next iteration
    dt.gibbs_energy_dist(m,T);
    if (m.n->n>0.0) {
      double phi_old=dt.get_phi(m,T);
      dt.gibbs_fixp_neutron(m.pr*acc_inc_factor,
			    m.n->n*(1.0-phi_old),m,T,0.0);
    } else {
      dt.gibbs_fixp(m.pr*acc_inc_factor,m,T);
    }
    dt.gibbs_energy_dist(m,T);
  }
    
  // ------------------------------------------------------
  // Add nuclear matter to the table

  cout << "Starting nuclear matter." << endl;
  dist_tab.new_column("fr_nm");
  dist_tab.new_column("ed_nm");
  dist_tab.new_column("pr_nm");
  dist_tab.set_unit("fr_nm","1/fm^4");
  dist_tab.set_unit("ed_nm","1/fm^4");
  dist_tab.set_unit("pr_nm","1/fm^4");

  for(size_t i=0;i<dist_tab.get_nlines();i++) {
    double nbx=dist_tab.get("nb",i);
    if (nbx<6.0e-3) {
      dist_tab.set("fr_nm",i,0.0);
      dist_tab.set("ed_nm",i,0.0);
      dist_tab.set("pr_nm",i,0.0);
    } else {
      nm.nb=nbx;
      nm.T=T;
      nmt.calc(nm);
      dist_tab.set("fr_nm",i,nm.fr);
      dist_tab.set("ed_nm",i,nm.ed);
      dist_tab.set("pr_nm",i,nm.pr);
    }
  }
  cout << "Done with nuclear matter." << endl;

  // ------------------------------------------------------
  // Write the final output file

  string fname2x=prefix+"_acc.o2";
  hdf_file hf;
  hf.open_or_create(fname2x);
  hdf_output(hf,dist_tab,"acc");
  hf.close();

  return 0;
}

int crust_driver::set_tptr(std::vector<std::string> &sv, bool itive_com) {
  if (sv.size()<2) {
    cout << "Need temperature argument." << endl;
    return o2scl::exc_efailed;
  }
  double argt=o2scl::stod(sv[1]);
  if (argt<0.0) {
    cout << "Negative temperature. Using zero." << endl;
    argt=0.0;
  }
  Tptr=cng.convert("K","1/fm",argt);
  if (verbose>0) {
    cout << "Temperature is set to " << argt << " K or "
	 << Tptr << " 1/fm." << endl;
  }
  return 0;
}

int crust_driver::check_fun(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<2) {
    cout.setf(ios::left);
    size_t width=24;
    cout.width(width);
    cout << "check_none: " << check_none << endl;
    cout.width(width);
    cout << "check_mass_density: " << check_mass_density << endl;
    cout.width(width);
    cout << "check_free_energy_sna: " << check_free_energy_sna << endl;
    cout.width(width);
    cout << "check_free_energy_dist: " << check_free_energy_dist << endl;
    cout.width(width);
    cout << "check_free_energy_nm: " << check_nm_from_dist << endl;
    cout.width(width);
    cout << "check_ldrop_derivs: " << check_ldrop_derivs << endl;
    cout.width(width);
    cout << "check_mixture: " << check_mixture << endl;
    cout.width(width);
    cout << "check_sfactor: " << check_sfactor << endl;
    cout.width(width);
    cout << "check_free_energy_cell: " << check_free_energy_cell << endl;
    cout.width(width);
    cout << "check_pressure: " << check_pressure << endl;
    cout.width(width);
    cout << "check_rate2: " << check_rate2 << endl;
    cout.width(width);
    cout << "check_mass_fit: " << check_mass_fit << endl;
    cout.width(width);
    cout << "check_feq_numbers: " << check_feq_numbers << endl;
    cout.width(width);
    cout << "check_acc_numbers: " << check_acc_numbers << endl;
    cout.unsetf(ios::left);
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_mass_density) {
    
    cout << "Checking mass_density()." << endl;
    check=check_mass_density;
    
    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    matter m;
    init_dist("ashes",m);

    vector<string> a;
    a.push_back("");
    a.push_back("2.0e8");
    set_tptr(a,false);

    dt.check=check;
    dt.mass_density(m,Tptr);

    return 0;
  }

  if (o2scl::stoi(sv[1])==check_free_energy_sna) {

    test_mgr t;
    t.set_output_level(2);
    
    cout << "Checking free_energy_sna()." << endl;
    check=check_free_energy_sna;

    {
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
      
      snat.free_energy_sna(m,T);
      double fr1=m.fr;
      dt.free_energy_dist(m,T);
      double fr2=m.fr;
      t.test_rel(snat.part1,dt.part1,1.0e-12,"part1");
      t.test_rel(snat.part2,dt.part2,1.0e-12,"part2");
      t.test_rel(snat.part3,dt.part3,1.0e-12,"part3");
      t.test_rel(snat.part4,dt.part4,1.0e-12,"part4");
      t.test_rel(fr1,fr2,1.0e-12,"fr");
    }
    
    snat.check_free_energy_sna(dt,t);

    return 0;
  }
  
  if (o2scl::stoi(sv[1])==check_free_energy_dist) {

    cout << "\nChecking free_energy_dist().\n" << endl;
    check=check_free_energy_dist;
    
    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    cout << "Initializing for X-ray bust ashes distribution." << endl;
    matter m;
    init_dist("ashes",m);
    for(size_t i=0;i<m.dist.size();i++) {
      m.dist[i].n*=1.0e5;
      // Make sure these nuclei aren't in the mass table
      m.dist[i].N+=2*m.dist[i].Z;
    }
    
    m.n->n=0.01;
    cout << "Set external neutron density to " << m.n->n << ".\n" << endl;
    
    vector<string> a;
    a.push_back("");
    a.push_back("2.0e8");
    set_tptr(a,false);

    dt.check=check;
    dt.free_energy_dist(m,Tptr,false);
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_nm_from_dist) {
    
    cout << "\nChecking nm_thermo::calc_nm_from_dist().\n" << endl;

    check=check_nm_from_dist;

    test_mgr t;
    t.set_output_level(2);

    matter m, nm;
    nucleus nuc;
    nuc.Z=26;
    nuc.N=30;
    nuc.n=1.1e-6;
    m.dist.push_back(nuc);
    
    vector<string> a;
    a.push_back("");
    a.push_back("2.0e10");
    set_tptr(a,false);

    dt.check=check;
    m.T=Tptr;
    m.n->n=0.0;
    m.p->n=0.0;
    nmt.calc_nm_from_dist(m,nm);

    t.test_rel(nm.n->n,30.0*1.1e-6,1.0e-12,"nm_from_dist nn");
    t.test_rel(nm.p->n,26.0*1.1e-6,1.0e-12,"nm_from_dist np");

    t.report();
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_ldrop_derivs) {

    cout << "\nChecking ldrop_crust derivatives().\n" << endl;

    check=check_ldrop_derivs;

    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    lda.test_derivatives();
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_mixture) {

    cout << "\nChecking pycnonuclear fusion in a mixture().\n" << endl;

    check=check_mixture;
    pyc.test_mixture();
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_sfactor) {
    check=check_sfactor;
    pyc.test_Sfactor();
    pyc.test_rate();
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_free_energy_cell) {

    check=check_free_energy_cell;

    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);
    
    check_free_energy_cell_fun();
    
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_pressure) {

    check=check_pressure;

    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    matter m;
    dt.check_pressure(m,0.0);
    
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_rate2) {
    check=check_rate2;
    pyc.test_rate2();
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_mass_fit) {
    
    test_mgr t;
    t.set_output_level(2);
    
    check=check_mass_fit;
    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    // Load initial guess from mass formula
    ubvector x(lda.nfit), y(lda.nfit);
    lda.guess_fun(lda.nfit,x);

    std::vector<std::string> sv3={"mf"};
    cf.perform_fit(sv3,false);
    lda.guess_fun(lda.nfit,y);
    t.test_rel_vec(lda.nfit,x,y,2.0e-5,"mass fit");

    t.report();
  }
  
  if (o2scl::stoi(sv[1])==check_feq_numbers) {

    test_mgr t;
    t.set_output_level(2);
    
    check=check_feq_numbers;
    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);
    std::vector<std::string> sv3={"feq","check","0.02","20","28",
				 "0.01","0.0201"};
    full_eq(sv3,false);

    o2scl_hdf::hdf_file hf;
    hf.open("check_feq.o2");
    string name;
    table_units<> tab;
    hdf_input(hf,tab,name);
    hf.close();

    t.test_rel(tab.get("rho",0),3.349758e13,1.0e-6,"feq rho");
    t.test_rel(tab.get("coul",0),3.270691e-1,1.0e-6,"feq coul");
    t.test_rel(tab.get("fr_x",0),9.577481e-2,1.0e-6,"feq fr_x");
    t.test_rel(tab.get("nun",0),4.703074,1.0e-6,"feq nun");
    t.test_rel(tab.get("pr",0),2.499303e-4,1.0e-6,"feq pr");

    /*
      std::vector<std::string> sv4={"feqd","check","0.02","12","44",
      "0.01","0.0201"};
      full_eq_dist(sv4,false);
      
      hf.open("check_feqd.o2");
      hdf_input(hf,tab,name);
      hf.close();
      cout << "Here4." << endl;
      
      t.test_rel(tab.get("rho",0),3.349758e13,1.0e-6,"feq rho");
      t.test_rel(tab.get("coul",0),3.270779e-1,1.0e-6,"feq coul");
      t.test_rel(tab.get("fr_x",0),9.571164e-2,1.0e-6,"feq fr_x");
      t.test_rel(tab.get("nun",0),4.703074,1.0e-6,"feq nun");
      t.test_rel(tab.get("pr",0),2.499303e-4,1.0e-6,"feq pr");
    */

    t.report();
    
    return 0;
  }

  if (o2scl::stoi(sv[1])==check_acc_numbers) {

    test_mgr t;
    t.set_output_level(2);
    
    check=check_acc_numbers;

    std::vector<std::string> sv1={"model","SLy4"};
    model(sv1,false);
    std::vector<std::string> sv2={"rf","data/SLy4_moller.fit"};
    cf.read_fit(sv2,false);

    std::vector<std::string> sv3={"acc","check"};
    acc(sv3,false);

    o2scl_hdf::hdf_file hf;
    hf.open("check_acc.o2");
    string name;
    table_units<> tab;
    hdf_input(hf,tab,name);
    hf.close();

    t.test_rel(tab.get("rho",0),1.882493e12,1.0e-6,"acc rho");
    t.test_rel(tab.get("nb",0),1.124159e-3,1.0e-6,"acc nb");
    t.test_rel(tab.get("fr_x",0),5.395385e-3,1.0e-6,"acc fr_x");
    t.test_rel(tab.get("pr",0),1.404249e-5,1.0e-6,"acc pr");
    t.test_rel(tab.get("t_C",0),8.618467e-12,1.0e-6,"acc t_C");
    t.test_rel(tab.get("Si56",0),2.007427e-5,1.0e-6,"acc Si56");
    t.test_gen(tab.get_nlines()==1,"acc nlines");

    t.report();

    return 0;
  }
  
  return 0;
}

int crust_driver::run(int argc, char *argv[]) {
  
  // ---------------------------------------
  // Specify command-line option object
    
#ifdef O2SCL_READLINE    
  cli_readline cl;
#else
  cli cl;
#endif
  cl.prompt="crust> ";
  cl.gnu_intro=false;
    
  // ---------------------------------------
  // Set options

  static const int nopt=13;
  comm_option_s options[nopt]={
    {0,"feq","Construct single-nucleus full_equilibrium crust.",
     1,6,"<filename prefix> <nb> <Z> <N> <nn> <nb_end>",
     ((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::full_eq),
     cli::comm_option_both},
    {0,"feq2","Construct single-nucleus full_equilibrium crust.",
     1,6,"<filename prefix> <nb> <Z> <N> <nn> <nb_end>",
     ((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::full_eq2),
     cli::comm_option_both},
    {0,"tov","Construct star from crust.",
     2,2,"<crust file> <output prefix>",((string)"Help here."),
     new comm_option_mfptr<tov_shear>(&tsh,&tov_shear::tov),
     cli::comm_option_both},
    {0,"shear","Construct star from crust.",
     2,2,"<crust file> <prefix>",((string)"Help here."),
     new comm_option_mfptr<tov_shear>(&tsh,&tov_shear::shear),
     cli::comm_option_both},
    {0,"feqd","full_equilibrium crust with distribution.",
     1,2,"<filename prefix>",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::full_eq_dist),
     cli::comm_option_both},
    {0,"model","Select model EOS (default is SLy4).",
     1,10,"<model>",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::model),
     cli::comm_option_both},
    {0,"make-table","",
     0,0,"",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::make_table),
     cli::comm_option_both},
    {0,"mf","Perform a mass fit.",
     0,1,"[file]",((string)"Help here."),
     new comm_option_mfptr<crust_fit>(&cf,&crust_fit::perform_fit),
     cli::comm_option_both},
    {0,"rf","Read mass fit.",
     1,1,"<file>",((string)"Help here."),
     new comm_option_mfptr<crust_fit>(&cf,&crust_fit::read_fit),
     cli::comm_option_both},
    {0,"acc","Construct accreted crust with distribution",
     1,3,"<prefix> [file density]",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::acc),
     cli::comm_option_both},
    {0,"test-ndrip","",0,0,"","",
     new comm_option_mfptr<crust_driver>(this,&crust_driver::test_ndrip),
     cli::comm_option_both},
    {0,"tptr","Set temperature",
     1,1,"<prefix>",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::set_tptr),
     cli::comm_option_both},
    {0,"check","Check.",
     0,1,"[index]",((string)"Help here."),
     new comm_option_mfptr<crust_driver>(this,&crust_driver::check_fun),
     cli::comm_option_both}
  };

  cl.set_comm_option_vec(nopt,options);

  // ---------------------------------------
  // Set parameters

  // Liquid drop parameters
    
  cli::parameter_bool p_exc_volume;
  p_exc_volume.b=&lda.exc_volume;
  p_exc_volume.help="Excluded volume correction (default true).";
  cl.par_list.insert(make_pair("lda.exc_volume",&p_exc_volume));

  cli::parameter_bool p_extra_corr;
  p_extra_corr.b=&lda.extra_corr;
  p_extra_corr.help="Add dense matter corrections to nuclei (default true).";
  cl.par_list.insert(make_pair("lda.extra_corr",&p_extra_corr));

  cli::parameter_bool p_use_ame;
  p_use_ame.b=&lda.use_ame;
  p_use_ame.help="Use AME nuclei where available (default true).";
  cl.par_list.insert(make_pair("lda.use_ame",&p_use_ame));

  cli::parameter_bool p_use_moller;
  p_use_moller.b=&lda.use_moller;
  p_use_moller.help="Use MNMSK nuclei where available (default true).";
  cl.par_list.insert(make_pair("lda.use_moller",&p_use_moller));

  cli::parameter_bool p_inc_shell;
  p_inc_shell.b=&lda.inc_shell;
  p_inc_shell.help="Shell effects (default true).";
  cl.par_list.insert(make_pair("lda.inc_shell",&p_inc_shell));

  cli::parameter_bool p_full_surface;
  p_full_surface.b=&lda.full_surface;
  p_full_surface.help="Full surface energy expression (default true).";
  cl.par_list.insert(make_pair("lda.full_surface",&p_full_surface));

  cli::parameter_bool p_new_skin_mode;
  p_new_skin_mode.b=&lda.new_skin_mode;
  p_new_skin_mode.help=((string)"New skin in bulk energy ")
    +"(probably not working; default false).";
  cl.par_list.insert(make_pair("lda.new_skin_mode",&p_new_skin_mode));

  cli::parameter_bool p_use_pasta;
  p_use_pasta.b=&use_pasta;
  p_use_pasta.help="Pasta (probably not working; default false).";
  cl.par_list.insert(make_pair("lda.use_pasta",&p_use_pasta));

  cli::parameter_double p_hd_exp;
  p_hd_exp.d=&lda.hd_exp;
  p_hd_exp.help="High-density correction exponent from Eq. 6 (default 5.0).";
  cl.par_list.insert(make_pair("lda.hd_exp",&p_hd_exp));

  cli::parameter_double p_hd_coeff;
  p_hd_coeff.d=&lda.hd_coeff;
  p_hd_coeff.help=((string)"High-density correction coefficient ")
    +"from Eq. 6 (default 0.5).";
  cl.par_list.insert(make_pair("lda.hd_coeff",&p_hd_coeff));

  // Equilibrium crust parameters

  cli::parameter_double p_feq_nb;
  p_feq_nb.d=&feq_nb;
  p_feq_nb.help="full eq initial nb (default 2.0e-10).";
  cl.par_list.insert(make_pair("feq_nb",&p_feq_nb));

  cli::parameter_int p_feq_N;
  p_feq_N.i=&feq_N;
  p_feq_N.help="full eq initial N (default 28).";
  cl.par_list.insert(make_pair("feq_N",&p_feq_N));

  cli::parameter_int p_feq_Z;
  p_feq_Z.i=&feq_Z;
  p_feq_Z.help="full eq initial Z (default 34).";
  cl.par_list.insert(make_pair("feq_Z",&p_feq_Z));

  cli::parameter_double p_feq_nn;
  p_feq_nn.d=&feq_nn;
  p_feq_nn.help="full eq initial neutron density (default 1.0e-10).";
  cl.par_list.insert(make_pair("feq_nn",&p_feq_nn));

  cli::parameter_double p_feq_end;
  p_feq_end.d=&feq_end;
  p_feq_end.help="full eq ending baryon density (default 0.09).";
  cl.par_list.insert(make_pair("feq_end",&p_feq_end));

  cli::parameter_double p_feq_dnb;
  p_feq_dnb.d=&feq_dnb;
  p_feq_dnb.help="full eq density step (default 1.1).";
  cl.par_list.insert(make_pair("feq_dnb",&p_feq_dnb));

  // Fusion parameters

  cli::parameter_bool p_allow_highZ;
  p_allow_highZ.b=&pyc.allow_highZ;
  p_allow_highZ.help="Allow fusion for high Z (default true)";
  cl.par_list.insert(make_pair("pyc.allow_highZ",&p_allow_highZ));

  cli::parameter_bool p_use_fit;
  p_use_fit.b=&pyc.use_fit;
  p_use_fit.help="Use Yakovlev's S factor fit (default false).";
  cl.par_list.insert(make_pair("pyc.use_fit",&p_use_fit));

  cli::parameter_bool p_simple_pyc;
  p_simple_pyc.b=&simple_pyc;
  p_simple_pyc.help="Simple pycnonuclear fusion (default false).";
  cl.par_list.insert(make_pair("simple_pyc",&p_simple_pyc));

  cli::parameter_double p_Mdot;
  p_Mdot.d=&rn.mdot;
  p_Mdot.help="Accretion rate (1.0e-10).";
  cl.par_list.insert(make_pair("Mdot",&p_Mdot));

  // Fit parameters

  cli::parameter_bool p_fit_moller;
  p_fit_moller.b=&cf.fit_moller;
  p_fit_moller.help="Fit moller (default true).";
  cl.par_list.insert(make_pair("fit_moller",&p_fit_moller));

  // Accreted crust parameters

  cli::parameter_bool p_calc_heating;
  p_calc_heating.b=&calc_heating;
  p_calc_heating.help="Calculate heating (default true).";
  cl.par_list.insert(make_pair("calc_heating",&p_calc_heating));

  cli::parameter_bool p_more_reactions;
  p_more_reactions.b=&more_reactions;
  p_more_reactions.help="More reactions.";
  cl.par_list.insert(make_pair("more_reactions",&p_more_reactions));

  cli::parameter_bool p_allow_palpha;
  p_allow_palpha.b=&allow_palpha;
  p_allow_palpha.help="Output flow (default false).";
  cl.par_list.insert(make_pair("allow_palpha",&p_allow_palpha));

  cli::parameter_bool p_ec_one_cell;
  p_ec_one_cell.b=&rn.ec_one_cell;
  p_ec_one_cell.help="Simple pyc (default true).";
  cl.par_list.insert(make_pair("ec_one_cell",&p_ec_one_cell));

  cli::parameter_double p_ec_heating;
  p_ec_heating.d=&ec_heating;
  p_ec_heating.help="Electron capture heating fraction (default 0.25).";
  cl.par_list.insert(make_pair("ec_heating",&p_ec_heating));

  cli::parameter_double p_acc_inc_factor;
  p_acc_inc_factor.d=&acc_inc_factor;
  p_acc_inc_factor.help="Pressure increase in crust (default 1.01).";
  cl.par_list.insert(make_pair("acc_inc_factor",&p_acc_inc_factor));

  cli::parameter_string p_dist_type;
  p_dist_type.str=&dist_type;
  p_dist_type.help="Distribution type (default nickel).";
  cl.par_list.insert(make_pair("dist_type",&p_dist_type));
    
  cli::parameter_double p_delta_n;
  p_delta_n.d=&rn.delta_n;
  p_delta_n.help="Fractional chunk size (default 0.01).";
  cl.par_list.insert(make_pair("delta_n",&p_delta_n));

  cli::parameter_double p_rho_summary;
  p_rho_summary.d=&rho_summary;
  p_rho_summary.help=((std::string)"If the mass density is larger ")
    +"than this value, summarize electron captures (default 1.0e20).";
  cl.par_list.insert(make_pair("rho_summary",&p_rho_summary));

  // Shear parameters
    
  cli::parameter_int p_hd_flag;
  p_hd_flag.i=&tsh.hd_flag;
  p_hd_flag.help="High-density EOS flag for shear modes (default 0).";
  cl.par_list.insert(make_pair("hd_flag",&p_hd_flag));
    
  cli::parameter_double p_mag_field;
  p_mag_field.d=&tsh.mag_field;
  p_mag_field.help="Magnetic field (in Gauss, default 0.0).";
  cl.par_list.insert(make_pair("mag_field",&p_mag_field));

  cli::parameter_double p_freq_mass;
  p_freq_mass.d=&tsh.freq_mass;
  p_freq_mass.help="Mass for frequency calculation (default 1.4).";
  cl.par_list.insert(make_pair("freq_mass",&p_freq_mass));

  cli::parameter_bool p_freq_gr;
  p_freq_gr.b=&tsh.freq_gr;
  p_freq_gr.help="Include GR in frequency calculation (default false).";
  cl.par_list.insert(make_pair("freq_gr",&p_freq_gr));

  // Other parameters

  cli::parameter_int p_check;
  p_check.i=&check;
  p_check.help="Check (default 0).";
  cl.par_list.insert(make_pair("check",&p_check));

  cli::parameter_int p_verbose;
  p_verbose.i=&verbose;
  p_verbose.help="verbose (default 1).";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  cli::parameter_bool p_feq_fix_mode;
  p_feq_fix_mode.b=&feq_fix_mode;
  p_feq_fix_mode.help="feq_fix_mode (default false).";
  cl.par_list.insert(make_pair("feq_fix_mode",&p_feq_fix_mode));

  cli::parameter_int p_fix_Z;
  p_fix_Z.i=&fix_Z;
  p_fix_Z.help="fix_Z (default 1).";
  cl.par_list.insert(make_pair("fix_Z",&p_fix_Z));

  cli::parameter_int p_fix_N;
  p_fix_N.i=&fix_N;
  p_fix_N.help="fix_N (default 1).";
  cl.par_list.insert(make_pair("fix_N",&p_fix_N));

  // ---------------------------------------
  // Process command-line arguments and run

  cl.run_auto(argc,argv);
    
  return 0;
}
