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
#include "pyc_rates.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

pyc_rates::pyc_rates() : cng(o2scl_settings.get_convert_units()) {
  Cpyc=3.9;
  Cexp=2.638;
  Cpl=1.25;
  debug=false;
  use_fit=false;

  set_fit();

  allow_highZ=true;
  loaded=false;
}

int pyc_rates::fusion_times(dist_thermo &dt, matter &m, double T, 
			    double &t_acc, double &t_c_fusion,
			    double &t_o_fusion, double &t_ne_fusion,
			    double &t_mg_fusion) {
  
  // Compute total number of nuclei per unit volume
  double ntot=0.0;
  for(size_t j=0;j<m.dist.size();j++) {
    ntot+=m.dist[j].n;
  }
  
  dt.gibbs_energy_dist(m,T);
  
  // Compute the mass fraction in nuclei
  double phi=dt.get_phi(m,T);
  // Subtract out contribution from electrons and quasi-free
  // neutrons
  double rho_fm4=cng.convert("g/cm^3","1/fm^4",m.rho);
  double XN=(rho_fm4-m.e.n*m.e.m-m.n->n*(1.0-phi)*m.n->m)/rho_fm4;
  
  // Estimate the mass of Carbon
  double m_carbon=12.0*931.0/197.33;
  double m_oxygen=16.0*931.0/197.33;
  double m_neon=20.0*931.0/197.33;
  double m_magnesium=24.0*931.0/197.33;

  get_times(6,6,12,12,m_carbon,m_carbon,ntot,ntot,ntot,
	    m.average_Z(),m.average_A(),m.rho,XN,1.0e-10,m.pr,
	    t_c_fusion,t_acc);
  get_times(8,8,16,16,m_oxygen,m_oxygen,ntot,ntot,ntot,
	    m.average_Z(),m.average_A(),m.rho,XN,1.0e-10,m.pr,
	    t_o_fusion,t_acc);
  get_times(10,10,20,20,m_neon,m_neon,ntot,ntot,ntot,
	    m.average_Z(),m.average_A(),m.rho,XN,1.0e-10,m.pr,
	    t_ne_fusion,t_acc);
  get_times(12,12,24,24,m_magnesium,m_magnesium,ntot,ntot,ntot,
	    m.average_Z(),m.average_A(),m.rho,XN,1.0e-10,m.pr,
	    t_mg_fusion,t_acc);

  if (false && m.rho>1.0e12) {

    debug=true;
    get_times(6,6,12,12,m_carbon,m_carbon,ntot,ntot,ntot,
	      m.average_Z(),m.average_A(),m.rho,XN,1.0e-10,m.pr,
	      t_c_fusion,t_acc);
    debug=false;

    cout << "m_carbon: " << m_carbon << endl;
    cout << "ntot: " << ntot << endl;
    cout << "XN: " << XN << endl;
    cout << "pr: " << m.pr << endl;
    cout << "t_acc: " << t_acc << endl;
    cout << "t_C_fusion: " << t_c_fusion << endl;
    exit(-1);
  }

  return 0;
}

int pyc_rates::load_data(string dir) {
  hdf_file hf;
  string fn="indata/beard10.o2";
  string name;
  hf.open(fn);
  hdf_input(hf,b10,name);
  hf.close();
  loaded=true;
  return 0;
}

int pyc_rates::set_fit() {
  double temp[10][13]=
    {{6,6,7.4836,0.1759,0.1759,0.0040,0.0040,0.0400,
      1.3736,3.5499,0.2658,0.35,0.11},
     {6,8,7.8671,0.1740,0.1280,-0.0045,-0.0310,0.0412,
      1.5438,5.2576,0.2306,0.47,0.12},
     {6,10,7.9387,0.1720,0.1206,-0.0171,-0.0035,0.0400,
      1.9478,3.7661,0.2328,0.62,0.16},
     {6,12,8.0513,0.1705,0.1014,-0.0210,-0.0186,0.0386,
      2.4327,4.0059,0.1844,0.56,0.17},
     {8,8,8.0641,0.1266,0.1266,-0.0377,-0.0377,0.0388,
      2.1998,6.0147,0.1547,0.50,0.14},
     {8,10,8.1191,0.1257,0.1183,-0.0461,-0.0068,0.0371,
      2.9486,3.5127,0.1702,0.65,0.19},
     {8,12,8.2404,0.1246,0.0994,-0.0500,-0.0216,0.0357,
      3.7433,2.8303,0.1417,0.65,0.20},
     {10,10,8.1419,0.1175,0.1175,-0.0107,-0.0107,0.0348,
      4.2215,0.1225,0.1717,1.00,0.27},
     {10,12,8.2880,0.1160,0.0987,-0.0157,-0.0273,0.0339,
      5.2525,0.2141,0.1342,0.86,0.28},
     {12,12,8.4509,0.0976,0.0976,-0.0288,-0.0288,0.0332,
      5.9785,0.5263,0.1393,0.82,0.28}};
  for(size_t i=0;i<10;i++) {
    for(size_t j=0;j<13;j++) {
      fit[i][j]=temp[i][j];
    }
  }
  return 0;
}

int pyc_rates::test_rate2() {
  
  // Reproducing Fig. 3 in Yakovlev06
  double m_carbon=12.0*931.0/197.33;
  double m_carbon_g=cng.convert("1/fm","g",m_carbon);
  double m_oxygen=16.0*931.0/197.33;
  double m_oxygen_g=cng.convert("1/fm","g",m_oxygen);
  double m_neon=20.0*931.0/197.33;
  double m_neon_g=cng.convert("1/fm","g",m_oxygen);
  double m_magnesium=24.0*931.0/197.33;
  double m_magnesium_g=cng.convert("1/fm","g",m_oxygen);

  table_units<> tab;
  tab.line_of_names("rho CC CO OO");
  tab.set_unit("rho","g/cm^3");
  tab.set_unit("CC","1/cm^3/s");
  tab.set_unit("CO","1/cm^3/s");
  tab.set_unit("OO","1/cm^3/s");
  
  // CC and OO seem to work, but CO doesn't
  use_fit=false;
  for(double lrho=8.5;lrho<=10.5001;lrho+=0.1) {
    double rho=pow(10.0,lrho);
    double nb_CC=rho/m_carbon_g/1.0e39;
    double R_CC=rate(6,6,12,12,m_carbon,m_carbon,
		     nb_CC,nb_CC,nb_CC,6,12,rho,1.0);
    double nb_CO=rho*2.0/(m_carbon_g+m_oxygen_g)/1.0e39;
    double R_CO=rate(6,8,12,16,m_carbon,m_oxygen,
		     nb_CO/2.0,nb_CO/2.0,nb_CO,7,14,rho,1.0);
    double nb_OO=rho/m_oxygen_g/1.0e39;
    double R_OO=rate(8,8,16,16,m_oxygen,m_oxygen,
		     nb_OO,nb_OO,nb_OO,8,16,rho,1.0);
    cout << rho << " " << R_CC*1.0e39 << " "
	 << R_CO*1.0e39 << " "
	 << R_OO*1.0e39 << endl;
    double line[4]={rho,R_CC*1.0e39,R_CO*1.0e39,R_OO*1.0e39};
    tab.line_of_data(4,line);
  }
  cout << endl;

  hdf_file hf;
  hf.open("y06fig3.o2");
  hdf_output(hf,tab,"y06");
  hf.close();

  // use_fit=true doesn't work yet
  use_fit=true;
  for(double lrho=8.5;lrho<=10.5001;lrho+=0.1) {
    double rho=pow(10.0,lrho);
    double nb_CC=rho/m_carbon_g/1.0e39;
    double R_CC=rate(6,6,12,12,m_carbon,m_carbon,
		     nb_CC,nb_CC,nb_CC,6,12,rho,1.0);
    double nb_CO=rho*2.0/(m_carbon_g+m_oxygen_g)/1.0e39;
    double R_CO=rate(6,8,12,16,m_carbon,m_oxygen,
		     nb_CO/2.0,nb_CO/2.0,nb_CO,7,14,rho,1.0);
    double nb_OO=rho/m_oxygen_g/1.0e39;
    double R_OO=rate(8,8,16,16,m_oxygen,m_oxygen,
		     nb_OO,nb_OO,nb_OO,8,16,rho,1.0);
    cout << rho << " " << R_CC*1.0e39 << " "
	 << R_CO*1.0e39 << " "
	 << R_OO*1.0e39 << endl;
  }
  cout << endl;
  
  // Now compute the rate at 10^12 g/cm^3 for pure carbon
  // and the corresponding fusion timescale
  use_fit=false;
  for(double rho=1.0e10;rho<1.1e14;rho*=sqrt(10.0)) {
    double nb_CC=rho/m_carbon_g/1.0e39;
    double R_CC=rate(6,6,12,12,m_carbon,m_carbon,
		     nb_CC,nb_CC,nb_CC,6,12,rho,1.0);
    cout << "At higher densities: " << rho << " " << R_CC*1.0e39 << " "
	 << nb_CC/R_CC << endl;
  }
  cout << endl;
  
  for(double rho=1.0e10;rho<1.1e14;rho*=sqrt(10.0)) {
    double nb_NN=rho/m_neon_g/1.0e39;
    double R_NN=rate(10,10,20,20,m_neon,m_neon,
		     nb_NN,nb_NN,nb_NN,10,20,rho,1.0);
    cout << "At higher densities: " << rho << " " << R_NN*1.0e39 << " "
	 << nb_NN/R_NN << endl;
  }
  cout << endl;
  
  for(double rho=1.0e10;rho<1.1e14;rho*=sqrt(10.0)) {
    double nb_MM=rho/m_magnesium_g/1.0e39;
    double R_MM=rate(12,12,24,24,m_magnesium,m_magnesium,
		     nb_MM,nb_MM,nb_MM,12,24,rho,1.0);
    cout << "At higher densities: " << rho << " " << R_MM*1.0e39 << " "
	 << nb_MM/R_MM << endl;
  }
  
  return 0;
}

int pyc_rates::test_Sfactor() {
  table<> t;
  t.line_of_names("E C10 C12 C16 C20 C24 C10a C12a C16a C20a C24a");

  for(double E=0.0;E<20.01;E+=0.25) {

    double line[11]={E,log10(Sfactor(6,6,10,10,E)),
		     log10(Sfactor(6,6,12,12,E)),
		     log10(Sfactor(6,6,16,16,E)),
		     log10(Sfactor(6,6,20,20,E)),
		     log10(Sfactor(6,6,24,24,E)),
		     0.0,0.0,0.0,0.0,0.0};
    use_fit=true;
    line[6]=log10(Sfactor(6,6,10,10,E));
    line[7]=log10(Sfactor(6,6,12,12,E));
    line[8]=log10(Sfactor(6,6,16,16,E));
    line[9]=log10(Sfactor(6,6,20,20,E));
    line[10]=log10(Sfactor(6,6,24,24,E));
    use_fit=false;

    cout << line[0] << " " << line[1] << " " << line[2] << " "
	 << line[3] << " " << line[4] << " " << line[5] << endl;
    t.line_of_data(11,line);
  }
  hdf_file hf;
  hf.open_or_create("Stest.o2");
  hdf_output(hf,t,"Stest");
  hf.close();
  return 0;
}

double pyc_rates::Sfactor(size_t Z1, size_t Z2, size_t A1, size_t A2,
			  double E) {

  if (loaded==false) {
    load_data("indata");
  }
  
  if (use_fit) {
    if (Z1>=6 && Z1<=12 && Z2>=6 && Z2<=12) {
      double A10=2*((double)Z1);
      double A20=2*((double)Z2);
      size_t row=100;
      for(size_t i=0;i<10;i++) {
	if (fit[i][0]==Z1 && fit[i][1]==Z2) row=i;
      }
      if (row==100) {
	cout << Z1 << " " << Z2 << " " << A1 << " " << A2 << endl;
	O2SCL_ERR("Row not found in Yakovlev's table.",o2scl::exc_efailed);
      }
      double R=fit[row][2], dR1, dR2;
      if (A1>=((size_t)(A10+1.0e-10))) {
	dR1=fit[row][3];
      } else {
	dR1=fit[row][5];
      }
      if (A2>=((size_t)(A20+1.0e-10))) {
	dR2=fit[row][4];
      } else {
	dR2=fit[row][6];
      }
      double delta=fit[row][7];
      double S0=fit[row][8];
      double xi0=fit[row][9];
      double xi1=fit[row][10];
      double alpha=Z1*Z2*o2scl_const::fine_structure;
      double RC=R+dR1*fabs(((double)A1)-A10)+
	dR2*fabs(((double)A2)-A20);
      double EC=alpha/RC;
      double xi=xi0+xi1*(A1+A2);
      double mu=(A1+A2)/A1/A2*931/hc_mev_fm;
      double ER=alpha*alpha*mu/2.0;
      if (E<EC) {
	double EC1=EC*(2.0+2.0*delta)/(2.0+3.0*delta);
	double RC1=RC*(1.0+delta);
	double beta=1.0/delta/(2.0+3.0*delta);
	double gamma=pow(2.0+3.0*delta,1.5)*sqrt(delta)/pow(1.0+delta,2.0);
	double sqt=sqrt(ER/EC);
	double sqt1=sqrt(ER/EC1);
	double g0r=8.0*sqt1;
	double g1r=-4.0/3.0/EC1*sqt1;
	double g2r=-0.2/EC1/EC1*sqt1;
	double xlo=sqrt(delta/(2.0+3.0*delta));
	double g0l=-gamma*sqt*(pi/2.0+asin(xlo)+xlo*sqrt(1.0-xlo*xlo));
	double g1l=gamma/EC*sqrt(pi/2.0+asin(xlo));
	double g2l=gamma/4.0/EC/EC*sqt*xlo/sqrt(1.0-xlo*xlo);
	return S0*exp((g0r+g0l)+E*(g1r+g1l)+E*E*(g1r*g1l));
      } 
      double eta=sqrt(ER/E);
      return S0*exp(2.0*pi*eta)*sqrt(E/EC)*(1.0+xi*(E-EC)/E);
    }
  }
    
  bool row_found=false;
  size_t row=0;
  for(size_t i=0;i<b10.get_nlines();i++) {
    if (((size_t)(b10.get("Z1",i)+1.0e-10))==Z1) {
      if (((size_t)(b10.get("Z2",i)+1.0e-10))==Z2) {
	if (((size_t)(b10.get("A1",i)+1.0e-10))==A1) {
	  if (((size_t)(b10.get("A2",i)+1.0e-10))==A2) {
	    row_found=true;
	    row=i;
	  }
	}
      }
    }
    if (row_found==true) i=b10.get_nlines();
  }

  if (row_found) {
    double B1=b10.get("B1",row);
    double B2=b10.get("B2",row);
    double B3=b10.get("B3",row);
    double C1=b10.get("C1",row);
    double C2=b10.get("C2",row);
    double C3=b10.get("C3",row);
    double C4=b10.get("C4",row);
    double EC=b10.get("EC",row);
    double D=b10.get("D",row);
    return exp(B1+B2*E+B3*E*E+(C1+C2*E+C3*E*E+C4*E*E*E)/
	       (1.0+exp((EC-E)/D)));
  }
    
  return 0.0;
    
}

double pyc_rates::Tp(size_t Z1, size_t Z2, double n12, double mu12) {
  return (sqrt(4.0*pi*o2scl_const::fine_structure*Z1*Z2*n12/2.0/mu12));
}
  
double pyc_rates::rB(size_t Z1, size_t Z2, double m1, double m2,
		     double mu12) {
  return (0.5/mu12/Z1/Z2/o2scl_const::fine_structure);
}
  
int pyc_rates::lambda_omega(size_t Z1, size_t Z2, 
			    double m1, double m2, double n1, double n2,
			    double &lambda, double &omega) {
  double a1=cbrt(3.0/n1/4.0/pi);
  double a2=cbrt(3.0/n2/4.0/pi);
  double a12=(a1+a2)/2.0;
  double n12=3.0/4.0/pi/a12/a12/a12;
  double mu12=m1*m2/(m1+m2);
  lambda=rB(Z1,Z2,m1,m2,mu12)*cbrt(n12/2.0);
  omega=Tp(Z1,Z2,n12,mu12);
  return 0;
}

int pyc_rates::lambda_omega2(size_t Z1, size_t Z2, size_t A1, size_t A2,
			     double m1, double m2, double n1, double n2,
			     double rho, double Xn, double aveZ, double aveA,
			     double &lambda, double &omega) {
  // This one doesn't work currently
  double a1=cbrt(3.0/n1/4.0/pi);
  double a2=cbrt(3.0/n2/4.0/pi);
  double a12=(a1+a2)/2.0;
  double n12=3.0/4.0/pi/a12/a12/a12;
  double mu12=m1*m2/(m1+m2);
  omega=Tp(Z1,Z2,n12,mu12);
  lambda=(((double)A1)+A2)/A1/A2/Z1/Z2/
    (cbrt(((double)Z1))+cbrt(((double)Z2)))*
    cbrt(rho*(1.0-Xn)*aveZ/aveA/1.3574e11);
  return 0;
}

int pyc_rates::test_rate() {

  // Mass and number density of Neon from Yakovlev's paper
  double mNE=31695.4/197.33;
  double nNE=1.809e-5;
  double rho=1.688e12;
  double Xn=0.39;
  double r=rate(10,10,34,34,mNE,mNE,nNE,nNE,nNE,10,34,rho,Xn);
  cout << "rate: " << r << endl;
  cout << "time: " << nNE/3.0/r << endl;
  for(rho=1.0e10;rho<=1.0e13;rho*=1.2) {
    double nB=rho/o2scl_cgs::mass_neutron/1.0e39;
    nNE=nB*(1.0-Xn)/34;
    double ne=10.0*nNE;
    double mue=cbrt(3.0*pi*pi*ne);
    cout << rho << " "
	 << nNE << " " 
	 << is_allowed(10,10,34,34,mNE,mNE,nNE,nNE,nNE,10,34,rho,Xn,
		       1.0e-10,mue) << endl;
  }
  cout << endl;

  cout << "N,Z dependence: " << endl;
  cout << "Vary N with fixed Z: " << endl;
  for(int Zx=6;Zx<=10;Zx+=2) {
    for (int Nx=4;Nx<=12;Nx+=2) {
      int Z1=Zx, N1=Nx;
      int Z2=Zx, N2=Nx;
      double m1=(Z1*938.3+N1*939.6-7*(N1+Z1))/197.33;
      double m2=(Z2*938.3+N2*939.6-7*(N2+Z2))/197.33;
      rho=1.0e12;
      Xn=0.3;
      double n1=1.0e-5, n2=1.0e-5;
      r=rate(Z1,Z2,N1,N2,m1,m2,n1,n2,n1,
	     (Z1+Z2)/2.0,(Z1+Z2+N1+N2)/2.0,rho,Xn);
      cout.width(2);
      cout << Z1 << " ";
      cout.width(2);
      cout << N1 << " ";
      cout.width(2);
      cout << Z2 << " ";
      cout.width(2);
      cout << N2 << " " << r << endl;
    }
    cout << endl;
  }
  cout << endl;

  cout << "Vary Z with fixed N: " << endl;
  for (int Nx=4;Nx<=12;Nx+=2) {
    for(int Zx=6;Zx<=12;Zx+=2) {
      for(int dn=0;dn<=8;dn+=2) {
	int Z1=Zx, N1=Nx;
	int Z2=Zx, N2=Nx+dn;
	double m1=(Z1*938.3+N1*939.6-7*(N1+Z1))/197.33;
	double m2=(Z2*938.3+N2*939.6-7*(N2+Z2))/197.33;
	rho=1.0e12;
	Xn=0.3;
	double n1=1.0e-5, n2=1.0e-5;
	r=rate(Z1,Z2,N1,N2,m1,m2,n1,n2,n1,
	       (Z1+Z2)/2.0,(Z1+Z2+N1+N2)/2.0,rho,Xn);
	if (r>0.0) {
	  cout.width(2);
	  cout << Z1 << " ";
	  cout.width(2);
	  cout << N1 << " ";
	  cout.width(2);
	  cout << Z2 << " ";
	  cout.width(2);
	  cout << N2 << " " << log10(r) << endl;
	}
      }
    }
    cout << endl;
  }
  cout << endl;

  size_t ix=0;
  for (int Nx=4;Nx<=12;Nx+=2) {
    for(int dn=0;dn<=8;dn+=2) {
      for(int dz=0;dz<=4;dz+=2) {
	string fname="sf_test/prt_"+
	  itos(Nx)+"_"+itos(Nx+dn)+"_"+itos(dz)+".out";
	ofstream fout(fname.c_str());
	  
	for(int Zx=6;Zx<=12;Zx+=2) {
	  int Z1=Zx, N1=Nx;
	  int Z2=Zx+dz, N2=Nx+dn;
	  double m1=(Z1*938.3+N1*939.6-7*(N1+Z1))/197.33;
	  double m2=(Z2*938.3+N2*939.6-7*(N2+Z2))/197.33;
	  rho=1.0e12;
	  Xn=0.3;
	  double n1=1.0e-5, n2=1.0e-5;
	    
	  double sf=Sfactor(Z1,Z2,N1+Z1,N2+Z2,2.0);
	    
	  if (sf>0.0) {
	    fout.width(2);
	    fout << Z1 << " ";
	    fout.width(2);
	    fout << N1 << " ";
	    fout.width(2);
	    fout << Z2 << " ";
	    fout.width(2);
	    fout << N2 << " " << log10(sf) << endl;
	  }
	}
	fout.close();
	ix++;
      }
    }
  }

  return 0;
}

double pyc_rates::rate(size_t Z1, size_t Z2, size_t A1, size_t A2,
		       double m1, double m2, double n1, double n2,
		       double ntot, double avgZ, double avgA, double rho, 
		       double XN) {

  double lambda, omega;
  double x1=n1/ntot;
  double x2=n2/ntot;
  if (true) {
    lambda_omega(Z1,Z2,m1,m2,n1,n2,lambda,omega);
  } else {
    // This one doesn't work currently
    lambda_omega2(Z1,Z2,A1,A2,m1,m2,n1,n2,rho,XN,avgZ,avgA,lambda,omega);
  }
  //x1=1.0;
  //x2=1.0;

  if (Z1>12) Z1=12;
  if (Z2>12) Z2=12;
    
  double Sf=Sfactor(Z1,Z2,A1,A2,omega);

  double r=8.0e7*Cpyc*rho*XN*x1*x2*A1*A2*avgA*Z1*Z1*Z2*Z2/
    pow(A1+A2,2.0)*Sf*pow(lambda,3.0-Cpl)*exp(-Cexp/sqrt(lambda));
  // Factor for identical nuclei
  if (Z1==Z2 && A1==A2) r/=2.0;

  if (debug) {
    cout << "rate() omega,rho,fact,S: " << omega << " " << rho << " " 
	 << XN*x1*x2*A1*A2*avgA*Z1*Z1*Z2*Z2/pow(A1+A2,2.0) << " " 
	 << Sfactor(Z1,Z2,A1,A2,omega) << endl;
    cout << "lambda,n1,n2,fac,r: " 
	 << lambda << " " << n1 << " " << n2 << " "
	 << pow(lambda,3.0-Cpl)*exp(-Cexp/sqrt(lambda)) << " "
	 << r << endl;
  }
  return r;
}

int pyc_rates::test_mixture() {

  size_t Z1=6, Z2=6, N1=6, N2=8, A1=12, A2=14;
  double m1=939.0*A1/197.33, m2=m1;//939.0*A2/197.33;
  
  cout.setf(ios::left);
  cout.width(10);
  cout << "x";
  cout.width(10);
  cout << "lambda";
  cout.width(10);
  cout << "omega";
  cout.width(10);
  cout << "Sf";
  cout.width(10);
  cout << "alpha";
  cout.width(10);
  cout << "rate";
  cout.width(10);
  cout << "dens_rat";
  cout.width(10);
  cout << "ratio" << endl;
  cout.unsetf(ios::left);

  for(double x=0.9;x>=0.099;x-=0.1) {

    // Set densities
    double ntot=1.0e-5;
    double n1=x*ntot;
    double n2=(1.0-x)*ntot;

    // Mass density
    double rho=(m1*n1+m2*n2)*3.5176722e14;

    double XN=0.1;

    double avgA=(A1*n1+A2*n2)/ntot;
    double avgZ=(Z1*n1+Z2*n2)/ntot;

    double lambda, omega;
    lambda_omega2(Z1,Z2,A1,A2,m1,m2,n1,n2,rho,XN,avgZ,avgA,lambda,omega);

    double Sf=Sfactor(Z1,Z2,A1,A2,omega);
    double alpha=pow(lambda,3.0-Cpl)*exp(-Cexp/sqrt(lambda));

    double r=rate(Z1,Z2,A1,A2,m1,m2,n1,n2,ntot,avgZ,avgA,rho,XN);

    double ratio=ntot/3.0/r;

    cout.precision(3);
    cout << x << " " << lambda << " " << omega << " " 
	 << Sf << " " << alpha << " "
	 << r << " " << pow(n1*n2/ntot/ntot,2.0) << " " 
	 << ratio << endl;
  }

  return 0;
}

void pyc_rates::get_times(size_t Z1, size_t Z2, size_t A1, size_t A2,
			  double m1, double m2, double n1, double n2,
			  double ntot, double avgZ, double avgA, double rho,
			  double XN, double Mdot, double P, double &t_fusion,
			  double &t_acc) {

  // Convert from fm^{-4} to g/cm/s^2
  P*=3.1615e35;
  
  // Compute column depth in g/cm^2 assuming a fixed surface gravity
  // of 2.43e14 cm/s^2
  double column=P/2.43e14;

  // Compute mass accretion rate per unit area assuming a radius of 12
  // km in g/cm^2/s. One year is approximately 3.1557e7 seconds.
  double sigma_dot=Mdot*o2scl_cgs::solar_mass/3.1557e7/
    (4.0*pi*1.2e6*1.2e6);
  
  // Compute the timescale in s
  t_acc=column/sigma_dot;

  double r=rate(Z1,Z2,A1,A2,m1,m2,n1,n2,ntot,avgZ,avgA,rho,XN);

  // Give a floor to the rate to avoid dividing by zero
  if (r<1.0e-100) r=1.0e-100;

  t_fusion=ntot/3.0/r;
  
  return;
}

bool pyc_rates::is_allowed(size_t Z1, size_t Z2, size_t A1, size_t A2,
			   double m1, double m2, double n1, double n2,
			   double ntot, double avgZ, double avgA, double rho,
			   double XN, double Mdot, double P) {

  if (allow_highZ==false && (Z1>12 || Z2>12)) return false;
    
  // Ensure Z1<Z2 and, if Z1=Z2, then also A1<A2
  if (Z2<Z1 || (Z1==Z2 && A2<A1) ) {
    int temp=Z1;
    Z1=Z2;
    Z2=temp;
    temp=A1;
    A1=A2;
    A2=temp;
  }

  // Ensure Z's and A's are even by increasing Z and decreasing A.
  // Rates tend to decrease with Z and increase with A, so this
  // choice minimizes fusion.
  if (Z1%2==1) Z1++;
  if (Z2%2==1) Z2++;
  if (A1%2==1) A1--;
  if (A2%2==1) A2--;

  if (true) {
    // Just assume Z<=4 always fuses
    if (Z1<=4 || Z2<=4) return true;
  } else {
    // If Z is 4, replace them with 6
    if (Z1==4) Z1==6;
    if (Z2==4) Z2==6;
  }

  if (!use_fit) {
    // If A>Amax, set A=Amax. This will increase fusion slightly for
    // heavier isotopes.
    if (Z1==6 && A1>24) A1=24;
    if (Z1==8 && A1>28) A1=28;
    if (Z1==10 && A1>40) A1=40;
    if (Z1==12 && A1>46) A1=46;
    if (Z2==6 && A2>24) A2=24;
    if (Z2==8 && A2>28) A2=28;
    if (Z2==10 && A2>40) A2=40;
    if (Z2==12 && A2>46) A2=46;
  }

  double t_fusion, t_acc;
  get_times(Z1,Z2,A1,A2,m1,m2,n1,n2,ntot,avgZ,avgA,rho,XN,Mdot,P,
	    t_fusion,t_acc);
  
  if (t_fusion<t_acc) {
    return true;
  }
  return false;
}
