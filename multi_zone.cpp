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
#include "multi_zone.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

multi_zone::multi_zone() {

  cng.units_cmd_string="units -f indata/units_hck.dat ";
  //conv=cng.convert("1/fm^4","Msun/km^3",1.0);
  conv=1.7684743e-4;
  
  schwarz_km=o2scl_mks::schwarzchild_radius/1.0e3;

  sdebug=0;
}

double multi_zone::min_fun(size_t nv, const ubvector &x) {

  // Set neutron densities for all zones but the bottom zone
  for(size_t i=0;i<n_zones-1;i++) {
    zones[i]->n->n=x[i]*min_scale[i];
    if (zones[i]->n->n<0.0) return 1.0e3;
  }
  
#ifdef O2SCL_NEVER_DEFINED
  
  mm_funct_mfptr<multi_zone> mfm(this,&multi_zone::solve_fun);

  // Ensure the initial guess is the same for all calls to 
  // min_fun()
  x_solve.set_all(1.0);

  gmr.msolve(n_solve,x_solve,mfm);

  // Check if pressure is non-monotonic
  for(size_t i=0;i<n_zones-1;i++) {
    if (P[i+1]<P[i]) gmtot+=100.0;
  }
  if (P[n_zones-1]>Pbottom) gmtot+=100.0;
  if (P[0]<Ptop) gmtot+=100.0;
  
#endif

  // Minimize the gravitational mass
  return gmtot;
}

int multi_zone::solve_fun(size_t nv, const ubvector &x, 
			  ubvector &y) {
  
  if (sdebug>1) {
    cout << "nv,n_solve: " << nv << " " << n_solve << endl;
  }

  // Set quasi-free neutron density for the bottom zone
  size_t ix=0;
  zones[n_zones-1]->n->n=x[ix]*solve_scale[ix];
  if (zones[n_zones-1]->n->n<0.0) return 1;
  ix++;

  // Set the nuclear densities
  for(size_t i=0;i<n_zones;i++) {
    for(size_t j=0;j<zones[i]->dist.size();j++) {
      zones[i]->dist[j].n=x[ix]*solve_scale[ix];
      if (zones[i]->dist[j].n<0.0) return 2;
      ix++;
    }
  }

  // Determine pressure and energy densities
  for (size_t i=0;i<n_zones;i++) {
    dtp->gibbs_energy_dist(*zones[i],Tptr);
    // Convert to Msun/km^4
    P[i]=zones[i]->pr*conv;
    eps[i]=zones[i]->fr*conv;
    if (sdebug>1) {
      cout << "Zone " << i << " " << zones[i]->n->n << " " 
	   << zones[i]->pr << " " << zones[i]->fr << " " 
	   << P[i] << " " << eps[i] << " "
	   << zones[i]->gb/zones[i]->nb << endl;
	   
    }
  }

  if (first_time) {
    Pbottom=zones[n_zones-1]->pr*
      sqrt(zones[n_zones-1]->pr/zones[n_zones-2]->pr)*conv;
    Ptop=zones[0]->pr/sqrt(zones[1]->pr/zones[0]->pr)*conv;
    if (sdebug>1) {
      cout << "Pbottom, Ptop: " << Pbottom << " " << Ptop << endl;
      cout << endl;
    }
  }

  // Determine radii and total nuclear numbers from thermodynamics
  ix=0;
  nntot=0.0;
  for (int i=n_zones-1;i>=0;i--) {

    // Construct the pressures on the radial grid points

    double Plo, Phi, dP;
    if (i==0) {
      Phi=sqrt(P[1]*P[0]);
      Plo=Ptop;
    } else if (i==((int)(n_zones-1))) {
      Phi=Pbottom;
      Plo=sqrt(P[n_zones-2]*P[n_zones-1]);
    } else {
      Phi=sqrt(P[i+1]*P[i]);
      Plo=sqrt(P[i-1]*P[i]);
    }

    if (Plo>Phi) return 4;

    if (sdebug>1) cout << "Plo,Phi: " << Plo << " " << Phi << endl;

    // Compute radii and gravitational masses

    double rlast, gmlast;
    if (i==((int)(n_zones-1))) {
      gmlast=1.4;
      rlast=11.0;
    } else {
      rlast=r[i+1];
      gmlast=gm[i+1];
    }
    double dPdr=-schwarz_km/2.0*gmlast*eps[i]/rlast/rlast/
      (1.0-schwarz_km*gmlast/rlast);
    r[i]=rlast-(Phi-Plo)/dPdr;
    gm[i]=gmlast+4.0*pi*rlast*rlast*eps[i];

    if (sdebug>1) cout << "r_i, gm_i: " << r[i] << " " << gm[i] << endl;

    // Compute numbers of neutrons and nuclei

    // 1.0e54 to convert from fm^{-3} to km^{-3}
    nn[i]=4.0*pi*pow(rlast,2.0)*zones[i]->n->n*1.0e54*
      sqrt(1.0-schwarz_km*gmlast/rlast)*(r[i]-rlast);
    nntot+=nn[i];

    for(size_t j=0;j<zones[i]->dist.size();j++) {
      
      // 1.0e54 to convert from fm^{-3} to km^{-3}
      nN[ix]=4.0*pi*pow(rlast,2.0)*zones[i]->dist[j].n*1.0e54*
	sqrt(1.0-schwarz_km*gmlast/rlast)*(r[i]-rlast);

      if (first_time) nN_base[ix]=nN[ix];
      
      if (sdebug>1) {
	cout << "nN: " << ix << " " << zones[i]->dist[j].n << " " 
	     << (r[i]-rlast) << " " << nN[ix] << " " << nN_base[ix] << endl;
      }
      
      y[ix]=(nN[ix]-nN_base[ix])/nN_base[ix];
      ix++;
      
    }

    // Total gravitational mass for minimization
    if (i==0) gmtot=gm[i];
    
    if (sdebug>1) cout << endl;
  }

  // Solve for equal number of neutrons
  if (first_time) {
    nn_base=nntot;
  }
  if (sdebug>1) {
    cout << "nn: " << nntot << " " << nn_base << endl;
  }
  y[ix]=(nntot-nn_base)/nn_base;
  ix++;

  // Check gmtot
  if (sdebug>1) {
    cout << "gmtot: " << gmtot << " " << ix << " " << nv << " " 
	 << n_solve << endl;
  }
  if (gmtot<0.0) return 3;

  return 0;
}

void multi_zone::test_stability(double T, dist_thermo &dt, 
				size_t nz, table_units<> &t, size_t fz) {
  
#ifdef O2SCL_NEVER_DEFINED

  cout.precision(6);

  Tptr=T;
  dtp=&dt;

  // ----------------------------------------------------------
  // Read table into zones object

  read_zones(nz,t,fz);
  
  // ----------------------------------------------------------
  // First compute the total number of nuclei
  
  ubvector x_min, y_solve;
  double y;

  for(size_t i=0;i<n_zones-1;i++) {
    x_min.push_back(zones[i]->n->n);
  }
  n_min=n_zones-1;

  x_solve.free();
  x_solve.push_back(zones[n_zones-1]->n->n);
  for(size_t i=0;i<n_zones;i++) {
    for(size_t j=0;j<zones[i]->dist.size();j++) {
      x_solve.push_back(zones[i]->dist[j].n);
    }
  }
  n_solve=x_solve.size();

  if (false) {
    cout << "n_zones, n_solve, n_min: " << n_zones << " " << n_solve << " " 
	 << n_min << endl;
  }

  // Resize all the vectors
  solve_scale.resize(n_solve);
  min_scale.resize(n_min);
  nn.resize(n_zones);
  nN.resize(n_solve-1);
  nN_base.resize(n_solve-1);
  y_solve.resize(n_solve);

  // Renormalize and set scales
  for(size_t i=0;i<n_min;i++) {
    min_scale[i]=x_min[i];
    x_min[i]=1.0;
  }
  for(size_t i=0;i<n_solve;i++) {
    solve_scale[i]=x_solve[i];
    x_solve[i]=1.0;
  }

  first_time=true;
  sdebug=0;
  y=solve_fun(n_solve,x_solve,y_solve);
  sdebug=0;
  first_time=false;

  {
    cout << 0.0 << " " << Ptop/conv << endl;
    for(size_t i=0;i<nz;i++) {
      cout << zones[i]->n->n << " " << zones[i]->pr << " " 
	   << zones[i]->gb/zones[i]->nb << " "
	   << zones[i]->rho 
	   << endl;
    }
    cout << 0.0 << " " << Pbottom/conv << endl;
    cout << gmtot << endl;
    cout << endl;
  }

  min_fun(n_min,x_min);

  // ----------------------------------------------------------
  // Now minimize

  for(size_t iii=0;iii<100;iii++) {
    
    ubvector step(1);
    step[0]=1.0e-4;
    gms.set_step(1,step);
    
    multi_funct_mfptr<multi_zone> mfm(this,&multi_zone::min_fun);
    gms.ntrial=100000;
    gms.verbose=0;
    gms.tol_abs=1.0e-9;
    gms.mmin(n_min,x_min,y,mfm);

    gmtot=min_fun(n_min,x_min);
    
    if (false) {
      sdebug=2;
      y=solve_fun(n_solve,x_solve,y_solve);
      sdebug=0;
    }
    
    {
      cout << 0.0 << " " << Ptop/conv << endl;
      for(size_t i=0;i<nz;i++) {
	cout << zones[i]->n->n << " " << zones[i]->pr << " " 
	     << zones[i]->gb/zones[i]->nb << " "
	     << zones[i]->rho 
	     << endl;
      }
      cout << 0.0 << " " << Pbottom/conv << endl;
      cout << gmtot << endl;
      cout << endl;
    }

  }

#endif
  
  return;
}

void multi_zone::read_zones(size_t nz, table_units<> &tab, size_t fz) {

  first_zone=fz;
  n_zones=nz;

  P.resize(n_zones);
  eps.resize(n_zones);
  r.resize(n_zones);
  gm.resize(n_zones);

  cout << "Reading table to get zone information." << endl;

  for(size_t i=0;i<n_zones;i++) {

    cout << "Zone " << i << "." << endl;

    // Create the matter object
    matter *mx=new matter(false);
    zones.push_back(mx);
    matter &m=*zones[i];

    // Set the quasi-free neutron density from the table
    m.n->n=tab.get("nn",fz+i);
    cout << "Neutron density: " << m.n->n << endl;

    double nptot=0.0;

    size_t firstcol=((size_t)(tab.get_constant("first_col")*1.001));
    cout << "Firstcol: " << firstcol << endl;

    for(size_t j=firstcol;j<tab.get_ncolumns();j++) {
      if (tab.get(j,fz+i)>0.0) {
	int Z, N, A;
	nmi.parse_elstring(tab.get_column_name(j),Z,N,A);
	nucleus nuc;
	nuc.Z=Z;
	nuc.N=N;
	nuc.n=tab.get(j,fz+i);
	m.dist.push_back(nuc);
	cout << "Found nucleus (" << Z << "," << N << ") with density "
	     << nuc.n << endl;
	nptot+=Z*nuc.n;
      }
    }

    // Set electron density
    m.e.n=nptot;
    cout << "Electron density: " << m.e.n << endl;
  }
  cout << endl;

  return;
}
