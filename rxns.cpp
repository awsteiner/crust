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
#include "rxns.h"
#include "crust.h"

using namespace std;
using namespace o2scl;
using namespace crust;
using namespace o2scl_const;

rxns::rxns() {
  delta_n=0.01;
  ec_one_cell=false;
  mdot=1.0e-10;

  fused_6=false;
  fused_8=false;
  fused_10=false;
  fused_12=false;
}

int rxns::add_missing(vector<nucleus> &trial, int newZ, int newN,
		      size_t &ix) {
  
  // Look for new nucleus in dist
  bool found=false;
  ix=0;
  for(size_t j=0;j<trial.size() && found==false;j++) {
    if (trial[j].Z==newZ && trial[j].N==newN) {
      found=true;
      ix=j;
    }
  }
    
  // If not found, add to trial distribution
  if (found==false) {
    nucleus new_nuc;
    trial.push_back(new_nuc);
    ix=trial.size()-1;
      
    trial[ix].Z=newZ;
    trial[ix].N=newN;
    // Ensure the new density is zero
    trial[ix].n=0.0;
  }

  return 0;
}

int rxns::ec_summary(crust_driver *a, matter &m, matter &m_new, 
		     double T, ubmatrix &gpb_store,
		     ubmatrix &m_store) {

  int Zmax=gpb_store.size1();
  int Nmax=gpb_store.size2();
  int Zmin=6;

  for(int tZ=Zmin;tZ<Zmax;tZ++) {
    for(int tN=tZ;tN<tZ*3 && tN<Nmax;tN++) {
      
      // Add parent nucleus if missing
      size_t ip;
      add_missing(m.dist,tZ,tN,ip);

      // Copy m.dist distribution, neutron, and proton densities
      m_new.dist=m.dist;
      m_new.n->n=m.n->n;
      m_new.p->n=m.p->n;
      
      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<m.dist.size();j++) {
	ntot+=m.dist[j].n;
      }

      // New nucleus after electron capture
      int newZ=m.dist[ip].Z-1;
      int newN=m.dist[ip].N+1;

      // Find parent in m_new and add daughter nucleus if missing
      size_t ip2, id;
      add_missing(m_new.dist,tZ,tN,ip2);
      add_missing(m_new.dist,newZ,newN,id);
      
      // Update densities
      double nshift=ntot*delta_n/10.0;
      m.dist[ip].n+=nshift;
      m_new.dist[id].n+=nshift;

      // Fix pressure of new distribution
      int ret=0;
      if (m.n->n>0.0) {
	a->dt.gibbs_energy_dist(m,T);
	double phi_old=a->dt.get_phi(m,T);
	ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),m_new,T,0.0);
      } else {
	a->dt.gibbs_energy_dist(m,T);
	a->dt.gibbs_fixp(m.pr,m_new,T);
      }
      
      if (ret==0) {
	
	double gpb1, gpb2;
	a->dt.baryon_density(m,T);
	gpb1=m.gb/m.nb;
	a->dt.baryon_density(m_new,T);
	gpb2=m_new.gb/m_new.nb;
	
	gpb_store(tZ,tN)=(gpb2-gpb1)*o2scl_const::hc_mev_fm;
	double mtot1=m_new.dist[ip2].m;
	double mtot2=m_new.dist[id].m;
	m_store(tZ,tN)=(mtot2-mtot1)*o2scl_const::hc_mev_fm;

	// End of loop for 'if (ret==0)'
      }

      // Remove extra nuclei from distributions
      m.dist[ip].n-=nshift;
      m_new.dist[id].n-=nshift;
      a->prune_distribution(m.dist);
      a->prune_distribution(m_new.dist);
      
      // End of loop for 'if (nshift>0)'
    }
    
    // End of loop if 'm.dist[i]>4'
  }
  
  return 0;
}

int rxns::elec_capture(crust_driver *a, matter &m, matter &m_new, 
		       double T, int &cnt, double &heat) {
  
  bool restart=true;
  while (restart==true) {

    restart=false;
    // Try an electron capture for each nucleus in the distribution
    for(int ii=m.dist.size()-1;ii>0;ii--) {
      size_t i=(size_t)ii;
      //for(size_t i=0;i<m.dist.size();i++) {
      
      // Ensure electron captures don't decrease Z too much
      if (m.dist[i].Z>4) {
	
	// Compute total number of nuclei per unit volume
	double ntot=0.0;
	for(size_t j=0;j<m.dist.size();j++) {
	  ntot+=m.dist[j].n;
	}

	// Copy m.dist distribution, neutron, and proton densities
	m_new.dist=m.dist;
	m_new.n->n=m.n->n;
	m_new.p->n=m.p->n;

	// New nucleus after electron capture
	int newZ=m.dist[i].Z-1;
	int newN=m.dist[i].N+1;
      
	// Find new nucleus, and add if missing
	size_t ix;
	add_missing(m_new.dist,newZ,newN,ix);

	// Update density
	double nshift=ntot*delta_n;

	// If there aren't enough parent nuclei for the requested
	// number of electron captures, then just capture all parent
	// nuclei
	if (nshift>m_new.dist[i].n) nshift=m_new.dist[i].n;
	
	m_new.dist[i].n-=nshift;
	m_new.dist[ix].n+=nshift;
	
	if (nshift>0.0) {
          
          // Fix pressure of new distribution
          int ret=0;
          if (m.n->n>0.0) {
            a->dt.gibbs_energy_dist(m,T);
            double phi_old=a->dt.get_phi(m,T);
            ret=a->dt.gibbs_fixp_neutron
	      (m.pr,m.n->n*(1.0-phi_old),m_new,T,0.0);
          } else {
            a->dt.gibbs_energy_dist(m,T);
            a->dt.gibbs_fixp(m.pr,m_new,T);
          }
          
          if (ret==0) {
            
            double gpb1, gpb2;
            
            if (m.pr>7.0e-6 && ec_one_cell) {
              a->dt.gibbs_energy_per_baryon_cell(m,T,m.pr,i,gpb1);
              a->dt.gibbs_energy_per_baryon_cell(m_new,T,m.pr,ix,gpb2);
              // Also compute baryon densities for later use
              a->dt.baryon_density(m,T);
              a->dt.baryon_density(m_new,T);
            } else {
              a->dt.baryon_density(m,T);
              gpb1=m.gb/m.nb;
              a->dt.baryon_density(m_new,T);
              gpb2=m_new.gb/m_new.nb;
            }
            
            if (gpb2<gpb1) {
	      // It's important that update_flow() is before
	      // dist_switch_gb() because dist_switch_gb() 
	      // rearranges the distribution
	      a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"ec",
			     (m.gb/m.nb-m_new.gb/m_new.nb)*
			     hc_mev_fm*a->ec_heating,m.rho);
	      double heat2=0.0;
              a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"ec");
              heat2*=a->ec_heating;
              heat+=heat2;
            }
            
          }
          
          // End of loop for 'if (nshift>0)'
        }

	// End of loop if 'm.dist[i]>4'
      }

      // End of loop over distribution
    }

    // End of restart loop
  }
    
  return 0;
}

int rxns::beta_decay(crust_driver *a, matter &m, matter &m_new, 
		     double T, int &cnt, double &heat) {
  
  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;

  bool restart=true;
  while (restart==true) {

    restart=false;
    // Try a beta decay for each nucleus in the distribution
    for(size_t i=0;i<current.size();i++) {

      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<current.size();j++) {
	ntot+=current[j].n;
      }

      // Copy current distribution
      trial=current;
      m_new.n->n=m.n->n;
      m_new.p->n=m.p->n;

      // New nucleus after beta decay
      int newZ=current[i].Z+1;
      int newN=current[i].N-1;
      
      // Find new nucleus, and add if missing
      size_t ix;
      add_missing(trial,newZ,newN,ix);

      // Update density
      double nshift=ntot*delta_n;

      // If there aren't enough parent nuclei for the requested
      // number of beta decays, then just allow all parent
      // nuclei to decay
      if (nshift>trial[i].n) nshift=trial[i].n;
      
      trial[i].n-=nshift;
      trial[ix].n+=nshift;
      
      if (nshift>0.0) {
	
	// Fix pressure of new distribution
	int ret=0;
	if (m.n->n>0.0) {
	  a->dt.gibbs_energy_dist(m,T);
	  double phi_old=a->dt.get_phi(m,T);
	  ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),
				       m_new,T,0.0);
	} else {
	  a->dt.gibbs_energy_dist(m,T);
	  a->dt.gibbs_fixp(m.pr,m_new,T);
	}
	
	if (ret==0) {
	  
	  a->dt.baryon_density(m,T);
	  double gpb1=m.gb/m.nb;
	  a->dt.baryon_density(m_new,T);
	  double gpb2=m_new.gb/m_new.nb;
	  
	  if (gpb2<gpb1) {
	    a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"bd",
			   (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,m.rho);
	    double heat2=0.0;
	    a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"bd");
	    heat+=heat2;
	  }
	  
	}
	
      }
      
    }
  }
  
  return 0;
}

int rxns::emit_neutron(crust_driver *a, matter &m, matter &m_new, 
		       double T, int &cnt, double &heat) {

  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;
  
  bool restart=true;
  while (restart==true) {

    restart=false;
    // Try a neutron emission for each nucleus in the distribution
    for(size_t i=0;i<current.size();i++) {

      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<current.size();j++) {
	ntot+=current[j].n;
      }

      // Copy current distribution
      trial=current;
      m_new.p->n=m.p->n;

      // New nucleus after neutron capture
      int newZ=current[i].Z;
      int newN=current[i].N-1;

      // Find new nucleus, and add if missing
      size_t ix;
      add_missing(trial,newZ,newN,ix);

      // Update densities
      double nshift=ntot*delta_n;

      // If there aren't enough parent nuclei for the requested
      // number of neutron emissions, then just emit a neutron for
      // all parent nuclei
      if (nshift>trial[i].n) nshift=trial[i].n;

      trial[i].n-=nshift;
      trial[ix].n+=nshift;
      m_new.n->n=m.n->n+nshift;

      if (nshift>0.0) {
	
	// Fix pressure of new distribution
	//gibbs_fixp(m,m_new,T);
	a->dt.gibbs_energy_dist(m,T);
	double phi_old=a->dt.get_phi(m,T);
      
	int ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old)+nshift,
					 m_new,T,nshift);
	
	if (ret==0) {
	  a->dt.baryon_density(m,T);
	  double gpb1=m.gb/m.nb;
	  a->dt.baryon_density(m_new,T);
	  double gpb2=m_new.gb/m_new.nb;

	  if (gpb2<gpb1) {
	    a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"en",
			   (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,m.rho);
	    double heat2=0.0;
	    a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"en");
	    heat+=heat2;
	  }
	}
	
      }

    }
  }
    
  return 0;
}

int rxns::en_summary(crust_driver *a, matter &m, matter &m_new, 
		     double T, ubmatrix &gpb_store,
		     ubmatrix &m_store) {

  int Zmax=gpb_store.size1();
  int Nmax=gpb_store.size2();
  int Zmin=6;

  for(int tZ=Zmin;tZ<Zmax;tZ++) {
    for(int tN=tZ;tN<tZ*3 && tN<Nmax;tN++) {
      
      // Add parent nucleus if missing
      size_t ip;
      add_missing(m.dist,tZ,tN,ip);

      // Copy m.dist distribution, neutron, and proton densities
      m_new.dist=m.dist;
      m_new.n->n=m.n->n;
      m_new.p->n=m.p->n;
      
      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<m.dist.size();j++) {
	ntot+=m.dist[j].n;
      }

      // New nucleus after electron capture
      int newZ=m.dist[ip].Z;
      int newN=m.dist[ip].N-1;

      // Find parent in m_new and add daughter nucleus if missing
      size_t ip2, id;
      add_missing(m_new.dist,tZ,tN,ip2);
      add_missing(m_new.dist,newZ,newN,id);
      
      // Update densities
      double nshift=ntot*delta_n/10.0;
      m.dist[ip].n+=nshift;
      m.dist[id].n+=nshift;
      m_new.n->n+=nshift;

      // Fix pressure of new distribution
      int ret=0;
      if (m.n->n>0.0) {
	a->dt.gibbs_energy_dist(m,T);
	double phi_old=a->dt.get_phi(m,T);
	ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),m_new,T,0.0);
      } else {
	a->dt.gibbs_energy_dist(m,T);
	a->dt.gibbs_fixp(m.pr,m_new,T);
      }
      
      if (ret==0) {
	
	double gpb1, gpb2;
	a->dt.baryon_density(m,T);
	gpb1=m.gb/m.nb;
	a->dt.baryon_density(m_new,T);
	gpb2=m_new.gb/m_new.nb;
	
	gpb_store(tZ,tN)=gpb2-gpb1;
	double mtot1=m_new.dist[ip2].m;
	double mtot2=m_new.dist[id].m;
	m_store(tZ,tN)=mtot2-mtot1;

	// End of loop for 'if (ret==0)'
      }

      // Remove extra nuclei and extra neutrons from distributions
      m.dist[ip].n-=nshift;
      m_new.dist[id].n-=nshift;
      m_new.n->n-=nshift;
      a->prune_distribution(m.dist);
      a->prune_distribution(m_new.dist);
      
      // End of loop for 'if (nshift>0)'
    }
    
    // End of loop if 'm.dist[i]>4'
  }
  
  return 0;
}

int rxns::neut_capture(crust_driver *a, matter &m, matter &m_new, 
		       double T, int &cnt, double &heat) {
  
  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;
  
  bool restart=true;
  while (restart==true) {

    restart=false;
    // Try an neutron capture for each nucleus in the distribution
    for(size_t i=0;i<current.size();i++) {

      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<current.size();j++) {
	ntot+=current[j].n;
      }

      // Copy current distribution
      trial=current;
      m_new.p->n=m.p->n;

      // New nucleus after neutron capture
      int newZ=current[i].Z;
      int newN=current[i].N+1;
      
      // Find new nucleus, and add if missing
      size_t ix;
      add_missing(trial,newZ,newN,ix);

      // Update density
      double nshift=ntot*delta_n;

      // If there aren't enough parent nuclei for the requested
      // number of neutron captures, then just capture all parent
      // nuclei
      if (nshift>trial[i].n) nshift=trial[i].n;

      // If there aren't enough neutrons, adjust accordingly
      if (nshift>m.n->n) nshift=m.n->n;
	
      trial[i].n-=nshift;
      trial[ix].n+=nshift;
      m_new.n->n=m.n->n-nshift;

      if (nshift>0.0) {
	
	// Fix pressure of new distribution
	//gibbs_fixp(m,m_new,T);
	a->dt.gibbs_energy_dist(m,T);
	double phi_old=a->dt.get_phi(m,T);

	int ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old)-nshift,
					 m_new,T,-nshift);
	
	if (ret==0) {
	  a->dt.baryon_density(m,T);
	  double gpb1=m.gb/m.nb;
	  a->dt.baryon_density(m_new,T);
	  double gpb2=m_new.gb/m_new.nb;
	  
	  if (gpb2<gpb1) {
	    a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"nc",
			   (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,m.rho);
	    double heat2=0.0;
	    a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"nc");
	    heat+=heat2;
	  }
	}
      }

    }
  }
    
  return 0;
}

int rxns::gen_reaction(crust_driver *a, matter &m, matter &m_new, 
		       double T, int &cnt) {

  a->dt.baryon_density(m,T);  
  double nb_base=m.nb;
  
  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;

  // Try an generic reaction for each nucleus in the distribution
  for(size_t i=0;i<current.size();i++) {
	
    // Compute loop sizes
    int oldZ=current[i].Z;
    int oldN=current[i].N;
    int deltaZ_lo, deltaZ_hi, deltaN_lo, deltaN_hi;
    a->delta_ZN(oldZ,oldN,deltaZ_lo,deltaZ_hi,deltaN_lo,deltaN_hi);

    for(int newZ=oldZ-deltaZ_lo;newZ<=oldZ+deltaZ_hi;newZ++) {
      for(int newN=oldN-deltaN_lo;newN<=oldN+deltaN_hi;newN++) {
	
	if (newZ!=oldZ || newN!=oldN) {
	  
	  if (newZ<=10) newZ=10;
	  if (newN<=10) newN=10;

	  bool restart=true;
	  while (restart==true) {

	    restart=false;

	    // Compute total number of nuclei per unit volume
	    double ntot=0.0;
	    for(size_t j=0;j<current.size();j++) {
	      ntot+=current[j].n;
	    }
	    
	    // Copy current distribution
	    trial=current;
	    m_new.n->n=m.n->n;
	    m_new.p->n=m.p->n;
	    
	    // Find new nucleus, and add if missing
	    size_t ix;
	    add_missing(trial,newZ,newN,ix);
	    
	    // Update densities
	    double nshift=ntot*delta_n;
	    
	    // If there aren't enough parent nuclei for the requested
	    // number of neutron emissions, then just emit a neutron for
	    // all parent nuclei
	    if (nshift>trial[i].n) nshift=trial[i].n;
	    
	    double oldA=oldZ+oldN;
	    double newA=newZ+newN;

	    trial[i].n-=nshift;
	    trial[ix].n+=nshift*oldA/newA;
	    
	  }
	}
      }
    }

  }
    
  return 0;
}

int rxns::pyc_fusion(crust_driver *a, matter &m, matter &m_new, 
		     double T, int &cnt, double &heat) {
  
  vector<nucleus> &current=m.dist;
  vector<nucleus> &trial=m_new.dist;

  bool restart=true;
  while (restart==true) {
    
    restart=false;

    // Try an fusion for each pair of nuclei in the distribution
    for(size_t i=0;i<current.size();i++) {
      for(size_t k=0;k<=i;k++) {

	// Ensure we don't fuse beyond shell gap limits
	if (current[i].N+current[k].N<406 && 
	    current[i].Z+current[k].Z<120) {

	  // Compute total number of nuclei per unit volume
	  double ntot=0.0;
	  for(size_t j=0;j<current.size();j++) {
	    ntot+=current[j].n;
	  }
	
	  // Determine if fusion is allowed
	  bool allowed=false;

	  if (a->simple_pyc) {

	    // Compute mass density
	    a->dt.mass_density(m,T);
	    double rho=m.rho;
	    
	    if (a->pycno_allowed(current[i].Z,current[k].Z,rho)) {
	      allowed=true;
	    }

	  } else {
	  
	    // Compute pressure (also computes rest mass energy density
	    // and baryon density)
	    a->dt.gibbs_energy_dist(m,T);

	    // Compute the mass fraction in nuclei
	    double phi=a->dt.get_phi(m,T);
	    // Subtract out contribution from electrons and quasi-free
	    // neutrons
	    double XN=(m.rho-m.e.n*m.e.m-m.n->n*(1.0-phi)*m.n->m)/m.rho;
	  
	    // Compute average Z
	    double avgZ=m.average_Z();
	  
	    // Compute average A
	    double avgA=m.average_A();
	    
	    if (a->pyc.is_allowed(current[i].Z,current[k].Z,
				  current[i].Z+current[i].N,
				  current[k].Z+current[k].N,
				  current[i].m,current[k].m,
				  current[i].n,current[k].n,
				  ntot,avgZ,avgA,m.rho,XN,mdot,m.pr)) {
	      allowed=true;
	    }
	  }

	  if (allowed==true) {

	    // Copy current distribution
	    trial=current;

	    // New nucleus after fusion
	    int newZ=current[i].Z+current[k].Z;
	    int newN=current[i].N+current[k].N;
      
	    // Find new nucleus, and add if missing
	    size_t ix;
	    add_missing(trial,newZ,newN,ix);

	    // Update density
	    double nshift=ntot*delta_n;

	    // If there aren't enough parent nuclei for the requested
	    // number of fusion reactions, then just fuse all nuclei
	    // present
	    if (i==k) {
	      if (nshift>trial[i].n/2.0) nshift=trial[i].n/2.0;
	    
	      trial[i].n-=2*nshift;
	      trial[ix].n+=nshift;
	    } else {
	      if (nshift>trial[i].n) nshift=trial[i].n;
	      if (nshift>trial[k].n) nshift=trial[k].n;
	    
	      trial[i].n-=nshift;
	      trial[k].n-=nshift;
	      trial[ix].n+=nshift;
	    }

	    if (trial[i].n<0.0 || trial[k].n<0.0 || trial[ix].n<0.0) {
	      cout << "Density negative in pyc_fusion." 
		   << endl;
	      cout << i << " " << k << " " << ix << endl;
	      cout << trial[i].n << " " << trial[k].n << " "
		   << trial[ix].n << " " << nshift << endl;
	      for(size_t j=0;j<m_new.dist.size();j++) {
		cout << j << " " << m_new.dist[j].Z << " " 
		     << m_new.dist[j].N << " "
		     << m_new.dist[j].n << endl;
	      }
	      cout << m_new.dist.size() << " " << m_new.dist[0].Z << " "
		   << m_new.dist[0].N << " " << m_new.n->n << endl;
	      cout << m_new.p->n << " " << m_new.e.n << " " << T << endl;
	      exit(-1);
	    }
	  
	    if (nshift>0.0) {
	    
	      // Fix pressure of new distribution
	      int ret=0;
	      if (m.n->n>0.0) {
		a->dt.gibbs_energy_dist(m,T);
		double phi_old=a->dt.get_phi(m,T);
		ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),
					     m_new,T,0.0);
	      } else {
		a->dt.gibbs_energy_dist(m,T);
		a->dt.gibbs_fixp(m.pr,m_new,T);
	      }

	      if (ret!=0 && m.rho>1.5e12) {
		cout << "Value ret nonzero at high density" << endl;
		exit(-1);
	      }
	      
	      if (ret==0) {
		
		a->dt.baryon_density(m,T);
		double gpb1=m.gb/m.nb;
		a->dt.baryon_density(m_new,T);
		double gpb2=m_new.gb/m_new.nb;
		
		if (gpb2<gpb1) {
		  
		  if (current[i].Z==6 && current[k].Z==6 && 
		      fused_6==false) {
		    cout << "Fusing Carbon." << endl;
		    fused_6=true;
		  }
		  if (current[i].Z==8 && current[k].Z==8 && 
		      fused_8==false) {
		    cout << "Fusing Oxygen." << endl;
		    fused_8=true;
		  }
		  if (current[i].Z==10 && current[k].Z==10 && 
		      fused_10==false) {
		    cout << "Fusing Neon." << endl;
		    fused_10=true;
		  }
		  if (current[i].Z==12 && current[k].Z==12 && 
		      fused_12==false) {
		    cout << "Fusing Magnesium." << endl;
		    fused_12=true;
		  }
		  
		  a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"pf",
				 (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,
				 m.rho);
		  double heat2=0.0;
		  a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"pf");
		  heat+=heat2;
		  
		  // End of (gpb2<gpb1) if statement
		}
		// End of (ret==0) if statement
	      }
	      // End of (nshift>0) if statement
	    }
	    // End of (allowed==true) if statement
	  }
	  // End of check to ensure N1+N2<406
	}
	// End of k loop
      }
      // End of i loop
    }
    // End of while (restart==true) loop
  }

  //char ch;
  //cin >> ch;
    
  return 0;
}

int rxns::emit_fragment(crust_driver *a, matter &m, matter &m_new, 
			double T, int &cnt, double &heat) {
  
  bool restart=true;
  while (restart==true) {

    restart=false;

    if (false) {

      // Try to emit a proton for each nucleus in the distribution
      for(size_t i=0;i<m.dist.size();i++) {

	// Compute total number of nuclei per unit volume
	double ntot=0.0;
	for(size_t j=0;j<m.dist.size();j++) {
	  ntot+=m.dist[j].n;
	}

	// Copy m.dist distribution, neutron, and proton densities
	m_new.dist=m.dist;
	m_new.n->n=m.n->n;
	m_new.p->n=m.p->n;

	// New nuclei
	int newZ=m.dist[i].Z-1;
	int newN=m.dist[i].N;
	int newZ2=1;
	int newN2=0;
      
	// Find new nucleus, and add if missing
	size_t ix;
	add_missing(m_new.dist,newZ,newN,ix);
	size_t ix2;
	add_missing(m_new.dist,newZ2,newN2,ix2);

	// Update density
	double nshift=ntot*delta_n;

	// If there aren't enough parent nuclei for the requested
	// number of electron captures, then just capture all parent
	// nuclei
	if (nshift>m_new.dist[i].n) nshift=m_new.dist[i].n;
      
	m_new.dist[i].n-=nshift;
	m_new.dist[ix].n+=nshift;
	m_new.dist[ix2].n+=nshift;

	if (nshift>0.0) {

	  // Fix pressure of new distribution
	  int ret=0;
	  if (m.n->n>0.0) {
	    a->dt.gibbs_energy_dist(m,T);
	    double phi_old=a->dt.get_phi(m,T);
	    ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),
					 m_new,T,0.0);
	  } else {
	    a->dt.gibbs_energy_dist(m,T);
	    a->dt.gibbs_fixp(m.pr,m_new,T);
	  }
	  
	  if (ret==0) {
	    
	    a->dt.baryon_density(m,T);
	    double gpb1=m.gb/m.nb;
	    a->dt.baryon_density(m_new,T);
	    double gpb2=m_new.gb/m_new.nb;
	    
	    if (gpb2<gpb1) {
	      cout << "Emitting proton" << endl;
	      exit(-1);
	      a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"ef",
			     (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,m.rho);
	      double heat2=0.0;
	      a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"ef");
	      heat+=heat2;
	    }
	    
	  }
	  
	}
	
      }

    }

    // Try to emit an alpha particle for each nucleus in the distribution
    for(size_t i=0;i<m.dist.size();i++) {

      // Compute total number of nuclei per unit volume
      double ntot=0.0;
      for(size_t j=0;j<m.dist.size();j++) {
	ntot+=m.dist[j].n;
      }

      // Copy m.dist distribution, neutron, and proton densities
      m_new.dist=m.dist;
      m_new.n->n=m.n->n;
      m_new.p->n=m.p->n;

      // New nuclei
      int newZ=m.dist[i].Z-2;
      int newN=m.dist[i].N-2;
      int newZ2=2;
      int newN2=2;
      
      // Find new nucleus, and add if missing
      size_t ix;
      add_missing(m_new.dist,newZ,newN,ix);
      size_t ix2;
      add_missing(m_new.dist,newZ2,newN2,ix2);

      // Update density
      double nshift=ntot*delta_n;

      // If there aren't enough parent nuclei for the requested
      // number of electron captures, then just capture all parent
      // nuclei
      if (nshift>m_new.dist[i].n) nshift=m_new.dist[i].n;
      
      m_new.dist[i].n-=nshift;
      m_new.dist[ix].n+=nshift;
      m_new.dist[ix2].n+=nshift;
      
      if (nshift>0.0) {
	
	// Fix pressure of new distribution
	int ret=0;
	if (m.n->n>0.0) {
	  a->dt.gibbs_energy_dist(m,T);
	  double phi_old=a->dt.get_phi(m,T);
	  ret=a->dt.gibbs_fixp_neutron(m.pr,m.n->n*(1.0-phi_old),
				       m_new,T,0.0);
	} else {
	  a->dt.gibbs_energy_dist(m,T);
	  a->dt.gibbs_fixp(m.pr,m_new,T);
	}
	
	if (ret==0) {
	  
	  a->dt.baryon_density(m,T);
	  double gpb1=m.gb/m.nb;
	  a->dt.baryon_density(m_new,T);
	  double gpb2=m_new.gb/m_new.nb;
	  
	  if (gpb2<gpb1) {
	    cout << "Emitting alpha" << endl;
	    exit(-1);
	    a->update_flow(m.dist[i].Z,m.dist[i].N,newZ,newN,"ef",
			   (m.gb/m.nb-m_new.gb/m_new.nb)*hc_mev_fm,m.rho);
	    double heat2=0.0;
	    a->dist_switch_gb(m,m_new,cnt,restart,heat2,T,"ef");
	    heat+=heat2;
	  }
	  
	}
	
	
      }
      
    }

  }
    
  return 0;
}
  
