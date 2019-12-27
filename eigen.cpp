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
#include "eigen.h"

using namespace std;
using namespace o2scl;
using namespace crust;
using namespace o2scl_const;
using namespace o2scl_hdf;

shear_eigen::shear_eigen() {
  ntab=400;
  r_grid.resize(ntab);
  xi_start.resize(2);
  xi.resize(ntab,2);
  dxi.resize(ntab,2);
  xi_err.resize(ntab,2);
}
  
int shear_eigen::derivs(double r, size_t nv, const ubvector &y, 
			ubvector &dydx) {

  double A, B;
  if (freq_gr) {
    A=shtab->interp("r",r,"Adsb");
    B=shtab->interp("r",r,"Bdsb");
  } else {
    A=shtab->interp("r",r,"Asw");
    B=shtab->interp("r",r,"Bsw");
  }

  dydx[0]=y[1];
  dydx[1]=-A*y[1]-B*y[0];

  return 0;
}

int shear_eigen::eigen_solve(size_t nv, const ubvector &x, 
			     ubvector &y) {

  shtab->set_constant("omega",x[0]);

  // Last coefficient
  shtab->function_column("-vs^2*(ell-1)*(ell+2)/R^2/(vs^2+vA^2)","B2sw");
  shtab->function_column("-e2l*vs^2*(ell-1)*(ell+2)/r^2/(vs^2+vA^2)","B2dsb");
  
  shtab->function_column("(1+vA^2)/(vs^2+vA^2)","B1sw");
  shtab->function_column("e2l/e2n*(1+vA^2)/(vs^2+vA^2)","B1dsb");

  shtab->function_column("(omega/c_km_s)^2*B1sw+B2sw","Bsw");
  shtab->function_column("(omega/c_km_s)^2*B1dsb+B2dsb","Bdsb");
  
  // Set units for last coefficient
  shtab->set_unit("B2sw","1/km^2");
  shtab->set_unit("Bsw","1/km^2");
  shtab->set_unit("B2dsb","1/km^2");
  shtab->set_unit("Bdsb","1/km^2");

  std::cout << "Missing 20." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED  
  
  ode_funct_mfptr<shear_eigen> ofm(this,&shear_eigen::derivs);

  // Initialize grid
  for(size_t i=0;i<ntab;i++) {
    r_grid[i]=Radc+(Rad-Radc)/((double)(ntab-1))*((double)i);
  }

  // Starting values of xi
  xi[0][0]=1.0;
  xi[0][1]=0.0;
  ois.ntrial=10000;
  ois.solve_grid<omatrix,omatrix_row>((Rad-Radc)/1.0e1,2,ntab,r_grid,
				      xi,xi_err,dxi,ofm);

  for(size_t i=0;i<ntab-1;i++) {
    nodes=0;
    if (xi[i][0]*xi[i+1][0]<0.0) {
      nodes++;
    }
  }
    
  y[0]=xi[ntab-1][1];
  
#endif
    
  return 0;
}

void shear_eigen::store_table(std::string fname) {
  
  ofstream fout(fname.c_str());
  table_units<> xitab;
  xitab.line_of_names("r xi xip");
  xitab.set_unit("r","km");
  xitab.set_unit("xip","1/km");
  for(size_t i=0;i<ntab;i++) {
    double line[3]={r_grid[i],xi(i,0),xi(i,1)};
    xitab.line_of_data(3,line);
  }

  hdf_file hf;
  hf.open(fname);
  hdf_output(hf,xitab,"xi");
  hf.close();
  
  return;
}

int shear_eigen::compute_freq(table_units<> *shear, double Rc, double R, 
			      double ell, double &freq_lo,
			      double freq_hi, int verbose) {
  
  shtab=shear;
  Radc=Rc;
  Rad=R;
    
  shtab->add_constant("omega",0.0);
  shtab->set_constant("ell",ell);

  shtab->sort_table("r");
  if (!shtab->is_column("sqrtB")) shtab->new_column("sqrtB");
  
  ubvector v_xip, v_omega;
    
  cout << "f omega xip integ nodes" << endl;

  ubvector x(1), y(1);
  double domega=pow(freq_hi/freq_lo,0.002);
  for(double omega=2.0*pi*freq_lo;omega<2.0*pi*freq_hi;omega*=domega) {
    
    x[0]=omega;
    eigen_solve(1,x,y);

    for(size_t i=0;i<shtab->get_nlines();i++) {
      if (freq_gr) {
	if (shtab->get("Bdsb",i)>0.0) {
	  shtab->set("sqrtB",i,sqrt(shtab->get("Bdsb",i)));
	} else {
	  shtab->set("sqrtB",i,0.0);
	}
      } else {
	if (shtab->get("Bsw",i)>0.0) {
	  shtab->set("sqrtB",i,sqrt(shtab->get("Bsw",i)));
	} else {
	  shtab->set("sqrtB",i,0.0);
	}
      }
	
    }

    std::cout << "Missing 21." << std::endl;
    exit(-1);
#ifdef O2SCL_NEVER_DEFINED
    v_omega.push_back(x[0]);
    v_xip.push_back(y[0]);
    
    if (verbose>0) {
      cout << x[0]/2.0/pi << " " << x[0] << " ";
      cout.setf(ios::showpos);
      cout << y << " ";
      cout << shtab->integ("r",shtab->get_constant("Rc"),
			   shtab->get_constant("R_solid"),"sqrtB") << " ";
      cout.unsetf(ios::showpos);
      cout << nodes << " ";
      cout << endl;
    }
#endif
  }

  // Look for the first row where the derivative changes sign
  size_t row=0;
  bool found=false;
  for(size_t i=0;i<v_omega.size()-1 && found==false;i++) {
    if (v_xip[i]*v_xip[i+1]<0.0) {
      row=i;
      found=true;
      cout << "Found: " << v_xip[i] << " " << v_xip[i+1] << endl;
    }
  }

  std::cout << "Missing 22." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED
  if (found) {
    int start=row-5, end=row+5;
    if (start<0) start=0;
    if (end>=((int)v_omega.size())) end=v_omega.size()-1;
    ubvector w_omega, w_xip;
    for(int i=start;i<end;i++) {
      w_omega.push_back(v_omega[i]);
      w_xip.push_back(v_xip[i]);
    }
    if (verbose>0) {
      for(size_t i=0;i<w_omega.size();i++) {
	cout << i << " " << w_omega[i] << " " << w_xip[i] << endl;
      }
    }
    freq_lo=si.interp(0.0,w_omega.size(),w_xip,w_omega)/2.0/pi;
    cout << "Final frequency: " << freq_lo << endl;

    // Now use the frequency to get the final values of xi
    x[0]=2.0*pi*freq_lo;
    eigen_solve(1,x,y);
    
    // Compute B^{1/2}
    for(size_t i=0;i<shtab->get_nlines();i++) {
      if (freq_gr) {
	if (shtab->get("Bdsb",i)>0.0) {
	  shtab->set("sqrtB",i,sqrt(shtab->get("Bdsb",i)));
	} else {
	  shtab->set("sqrtB",i,0.0);
	}
      } else {
	if (shtab->get("Bsw",i)>0.0) {
	  shtab->set("sqrtB",i,sqrt(shtab->get("Bsw",i)));
	} else {
	  shtab->set("sqrtB",i,0.0);
	}
      }
    }
    
    
  } else {
    freq_lo=0.0;
  }
#endif

  return 0;
}

int tov_shear::tov(std::vector<std::string> &sv, bool itive_com) {

  if (sv.size()<3) {
    cout << "No crust filename or no output prefix given." << endl;
    return o2scl::exc_efailed;
  }
    
  string prefix=sv[2];

  // Load crust table
  table_units<> feq;
  hdf_file hf;
  hf.open(sv[1]);
  hdf_input(hf,feq,"feq");
  hf.close();
  cout << "Loaded crust table '" << sv[1] << "'." << endl;

  // Load core table
  cout << "In_dir: " << in_dir << endl;
  table_units<> slb;
  hf.open(in_dir+((string)"/slb11.o2"));
  std::string name;
  hdf_input(hf,slb,name);
  hf.close();
  cout << "Loaded core table 'slb11.o2'." << endl;

  // Create new table for use in TOV
  table_units<> te;
  te.line_of_names("ed pr");
  te.set_unit("ed","1/fm^4");
  te.set_unit("pr","1/fm^4");

  double ed_last=0.0, pr_last=0.0;

  string slb_col="pr";
  if (hd_flag==1) slb_col="pr_p1";
  else if (hd_flag==3) slb_col="pr_m1";
  else if (hd_flag==2) slb_col="pr_p2";
  else if (hd_flag==4) slb_col="pr_m2";
  cout << "Using column: " << slb_col 
       << " in high-density EOS file." << endl;

  // Add crust rows
  for(size_t i=0;i<feq.get_nlines();i++) {
    double ed=feq.get("fr",i);
    double pr_crust=feq.get("pr",i);
    double pr_core=slb.interp("ed",ed*hc_mev_fm,
			      slb_col)/hc_mev_fm;
    if (hd_flag==0) {
      pr_core=(slb.interp("ed",ed*hc_mev_fm,"pr_m1")+
	       slb.interp("ed",ed*hc_mev_fm,"pr_p1"))/2.0/hc_mev_fm;
      
    }
    
    if (feq.get("nb",i)<0.1) {
      double line[2]={ed,pr_crust};
      ed_last=ed;
      pr_last=pr_crust;
      te.line_of_data(2,line);
      cout << feq.get("nb",i) << " "
	   << line[0] << " " << line[1] << " " << pr_core << endl;
    }
  }
  cout << "Switch." << endl;

  for(double ed=ed_last*1.1;ed<10.0;ed*=1.1) {
    double pr=slb.interp("ed",ed*hc_mev_fm,slb_col)/hc_mev_fm;
    if (hd_flag==0) {
      pr=(slb.interp("ed",ed*hc_mev_fm,"pr_m1")+
	  slb.interp("ed",ed*hc_mev_fm,"pr_p1"))/2.0/hc_mev_fm;
    }
    if (pr>pr_last) {
      double line[2]={ed,pr};
      te.line_of_data(2,line);
      cout << line[0] << " " << line[1] << endl;
    }
  }

  if (true) {
    if (mag_field>0.0) {
      for(size_t i=0;i<te.get_nlines();i++) {
	double corr=o2scl_const::gauss2_fm4*mag_field*mag_field/8.0/pi;
	if (corr>te["ed"][i]) corr=te["ed"][i];
	if (corr>te["pr"][i]) corr=te["pr"][i];
	te.set("ed",i,te.get("ed",i)+corr);
	te.set("pr",i,te.get("pr",i)+corr);
      }
    }
  }

  string fname1=prefix+"_tov_in.o2";
  hf.open(fname1);
  hdf_output(hf,te,"tov");
  hf.close();
  cout << "Wrote EOS for TOV in '" << fname1 << "'." << endl;

  eos_tov_interp teos;
  //teos.set_units("1/fm^4","1/fm^4","");
  teos.read_table(te,"ed","pr","");

  tov_solve ts;
  ts.verbose=0;
  ts.set_units("1/fm^4","1/fm^4","");
  ts.set_eos(teos);

  {
    cout << "Going to mvsr." << endl;
    ts.mvsr();
    cout << "Done with mvsr." << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname2=prefix+"_mvsr_out.o2";
    hf.open(fname2);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.0" << endl;
    ts.fixed(1.0);
    cout << "Done with 1.0" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();

    string fname3=prefix+"_m10.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.1" << endl;
    ts.fixed(1.1);
    cout << "Done with 1.1" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m11.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.2" << endl;
    ts.fixed(1.2);
    cout << "Done with 1.2" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m12.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.3" << endl;
    ts.fixed(1.3);
    cout << "Done with 1.3" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m13.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.4" << endl;
    ts.fixed(1.4);
    cout << "Done with 1.4" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
    cout << "Rad 1.4: " << mr->max("r") << endl;
      
    string fname3=prefix+"_m14.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.5" << endl;
    ts.fixed(1.5);
    cout << "Done with 1.5" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m15.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.6" << endl;
    ts.fixed(1.6);
    cout << "Done with 1.6" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m16.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.7" << endl;
    ts.fixed(1.7);
    cout << "Done with 1.7" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m17.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  {
    cout << "Going to 1.8" << endl;
    ts.fixed(1.8);
    cout << "Done with 1.8" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m18.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }
  if (false) {
    cout << "Going to 1.9" << endl;
    ts.fixed(1.9);
    cout << "Done with 1.9" << endl;
    shared_ptr<table_units<> > mr=ts.get_results();
      
    string fname3=prefix+"_m19.o2";
    hf.open(fname3);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }

  return 0;
}

int tov_shear::shear(std::vector<std::string> &sv, bool itive_com) {
    
  if (sv.size()<3) {
    cout << "Args: [crust filename] [prefix]" << endl;
    return o2scl::exc_efailed;
  }
  string crust_fn=sv[1];
  string prefix=sv[2];
    
  // -------------------------------------------
  // Load crust table

  table_units<> feq;
  hdf_file hf;
  hf.open(crust_fn);
  hdf_input(hf,feq,"feq");
  hf.close();
  if (verbose>0) cout << "Loaded crust table '" << sv[1] << "'." << endl;

  std::cout << "Missing 25." << std::endl;
  exit(-1);
#ifdef O2SCL_NEVER_DEFINED  
  if (false) {

    if (sv[1]=="prc_feq/SLy4_cs01_feq.o2") {
      cout << "Hack to remove rows with decreasing pressure in SLy4." << endl;
      feq.delete_row(feq.get_nlines()-1);
    }
    
    //feq.set_interp_type(itp_linear);
    
    // Add a row exactly at the transition
    ubvector line_trans(feq.get_ncolumns());
    for(size_t i=0;i<feq.get_ncolumns();i++) {
      if (feq.get_column_name(i)==((string)"nb")) {
	line_trans[i]=0.08;
      } else {
	line_trans[i]=feq.interp("nb",0.08,feq.get_column_name(i));
      }
    }
    feq.line_of_data(feq.get_ncolumns(),line_trans);

  }

  if (true) {
    
    // Extend crust table beyond that provided originally

    table_units<> ccnew;
    bool apr_mode=false;
    if (sv[1]=="prc_feq/SLy4_cs01_feq.o2") {
      hf.open("prc_feq/SLy4_core.o2");
      hdf_input(hf,ccnew);
      hf.close();
    } else if (sv[1]=="prc_feq/APR_cs01_feq.o2") {
      hf.open("prc_feq/APR_core.o2");
      hdf_input(hf,ccnew);
      hf.close();
      apr_mode=true;
    } else if (sv[1]=="prc_feq/Gs_cs01_feq.o2") {
      hf.open("prc_feq/Gs_core.o2");
      hdf_input(hf,ccnew);
      hf.close();
    } else if (sv[1]=="prc_feq/Rs_cs01_feq.o2") {
      hf.open("prc_feq/Rs_core.o2");
      hdf_input(hf,ccnew);
      hf.close();
    } else {
      O2SCL_ERR("Crust core error.",gsl_einval);
    }

    // Set to linear interpolation temporarily
    size_t itype=feq.get_interp_type();
    feq.set_interp_type(itp_linear);

    // Last row
    size_t row=feq.get_nlines()-1;

    // Density grid
    double nb_last=feq.get("nb",row);
    double nb_start=(floor(nb_last*100.0)+1)/100.0;

    for(double nb=nb_start;nb<0.16001;nb+=0.01) {
      if (apr_mode) {
	// nb, rho, N
	double line[34]={nb,feq.interp("nb",nb,"rho"),feq.get("N",row),
			 // Z, nn
			 feq.get("Z",row),feq.interp("nb",nb,"nn"),
			 // bulk, surf, coul, pair
			 0.0,0.0,0.0,0.0,
			 // shell, exc, BEoA, ede
			 0.0,0.0,0.0,0.0,
			 // ne, fr
			 0.0,ccnew.interp("nb",nb,"ed"),
			 // fr_check
			 ccnew.interp("nb",nb,"ed"),
			 // nnuc, chi, mnuc, fr_x
			 feq.interp("nb",nb,"nnuc"),0.0,0.0,0.0,
			 // nnL, npL, kfn, nun, nuN, pr
			 0.0,0.0,0.0,0.0,0.0,ccnew.interp("nb",nb,"pr"),
			 // Rws, Rn, Rp, Ncell
			 0.0,0.0,0.0,0.0,
			 // pr_x, nb_x, dPdYe, dfdYe
			 0.0,0.0,0.0,0.0};
	feq.line_of_data(34,line);
      } else {
	// nb, rho, N
	double line[39]={nb,feq.interp("nb",nb,"rho"),feq.get("N",row),
			 // Z, nn
			 feq.get("Z",row),feq.interp("nb",nb,"nn"),
			 // bulk, surf, coul, pair
			 0.0,0.0,0.0,0.0,
			 // shell, exc, BEoA, ede
			 0.0,0.0,0.0,0.0,
			 // ne, fr
			 0.0,ccnew.interp("nb",nb,"ed"),
			 // fr_check
			 ccnew.interp("nb",nb,"ed"),
			 // nnuc, chi, mnuc, fr_x
			 feq.interp("nb",nb,"nnuc"),0.0,0.0,0.0,
			 // nnL, npL, kfn, nun, nuN, pr
			 0.0,0.0,0.0,0.0,0.0,ccnew.interp("nb",nb,"pr"),
			 // Rws, Rn, Rp, Ncell
			 0.0,0.0,0.0,0.0,
			 // pr_x, nb_x, dPdYe, dfdYe
			 0.0,0.0,0.0,0.0,
			 // EoA, dEdnb
			 0.0,0.0,
			 // pr_int, ed_nm
			 0.0,ccnew.interp("nb",nb,"ed"),
			 // pr_nm
			 ccnew.interp("nb",nb,"pr")};
	feq.line_of_data(39,line);
      }
    }

    // Return interpolation type
    feq.set_interp_type(itype);
  }
  
  // -------------------------------------------
  // Load core table

  table_units<> slb;
  hf.open("indata/slb11.o2");
  hdf_input(hf,slb);
  hf.close();
  if (verbose>0) cout << "Loaded core table 'indata/slb.o2'." << endl;

  // Create new table for use in TOV
  table_units<> te;
  te.line_of_names("ed pr");

  // -------------------------------------------
  // Add crust rows

  // Record last crust pressure and energy density
  double ed_last=0.0, pr_last=0.0;

  for(size_t i=0;i<feq.get_nlines();i++) {
    double ed=feq.get("fr",i);
    double pr_crust=feq.get("pr",i);

    if (feq.get("nb",i)<0.08) {
      double line[2]={ed,pr_crust};
      ed_last=ed;
      pr_last=pr_crust;
      te.line_of_data(2,line);
      if (verbose>1) cout << line[0] << " " << line[1] << endl;
    }
  }
  if (verbose>1) cout << "Switch." << endl;

  // -------------------------------------------
  // Add core rows

  string slb_col="pr";
  if (hd_flag==1) slb_col="pr_p1";
  else if (hd_flag==3) slb_col="pr_m1";
  else if (hd_flag==2) slb_col="pr_p2";
  else if (hd_flag==4) slb_col="pr_m2";
  if (verbose>0) {
    cout << "Using column: " << slb_col 
	 << " in high-density EOS file." << endl;
  }

  for(double ed=ed_last*1.2;ed<10.0;ed*=1.02) {
    double pr=slb.interp("ed",ed*hc_mev_fm,slb_col)/hc_mev_fm;
    if (hd_flag==0) {
      pr=(slb.interp("ed",ed*hc_mev_fm,"pr_m1")+
	  slb.interp("ed",ed*hc_mev_fm,"pr_p1"))/2.0/hc_mev_fm;
    }
    if (pr>pr_last) {
      double line[2]={ed,pr};
      te.line_of_data(2,line);
      if (verbose>1) cout << line[0] << " " << line[1] << endl;
    }
  }
  
  // -------------------------------------------
  // Output EOS table which is given to TOV solver

  string tov_fn=prefix+"_tov_in.o2";
  hf.open(tov_fn);
  hdf_output(hf,te,"tov");
  hf.close();
  if (verbose>0) cout << "Wrote EOS for TOV to '" << tov_fn << "'." << endl;
  
  // -------------------------------------------
  // Perform TOV solution

  tov_interp_eos teos;
  teos.set_units("1/fm^4","1/fm^4");
  teos.read_table(te,"ed","pr");

  tov_solve ts;
  ts.verbose=0;
  if (true) {
    ts.hstart/=5.0;
    ts.hmax/=5.0;
    ts.hmin/=5.0;
  }
  ts.set_units("1/fm^4","1/fm^4");
  ts.set_eos(teos);

  if (verbose>1) {
    cout << "Going to TOV." << endl;
  }

  // Compute maximum mass and check 
  ts.mvsr();
  shared_ptr<table_units<> > mr0=ts.get_results();
  if (freq_mass>mr0->max("gm")) {
    cout << mr0->max("gm") << endl;
    O2SCL_ERR_RET("Requested mass greater than maximum mass.",
		  o2scl::exc_efailed);
  }

  // Compute full M vs. R curve if target mass is 1.0
  if (fabs(freq_mass-1.0)<1.0e-4) {
    
    ts.mvsr();
    shared_ptr<table_units<> > mr=ts.get_results();
    
    string mvsr_fn=prefix+"_mvsr.o2";
    hf.open(mvsr_fn);
    hdf_output(hf,*mr,"tov");
    hf.close();
  }

  // Compute structure for requested mass
  ts.fixed(freq_mass);
  shared_ptr<table_units<> > mr=ts.get_results();

  // Output solution 
  string mfixed_fn=prefix+"_mfixed.o2";
  hf.open(mfixed_fn);
  hdf_output(hf,*mr,"tov");
  hf.close();
    
  if (verbose>0) cout << "Wrote TOV output in '" << mfixed_fn << "'." << endl;
    
  // -------------------------------------------
  // Construct table to compute frequencies

  table_units<> *shtab=(table_units<> *)(&feq);

  //shtab->set_convert(cng);

  shtab->add_constant("hc",o2scl_const::hc_mev_fm);
  shtab->add_constant("alpha",o2scl_const::fine_structure);
  shtab->add_constant("pi",acos(-1.0));
  shtab->add_constant("eB",mag_field*o2scl_const::ec_gauss_fm2);
  shtab->add_constant("Rschwarz",o2scl_mks::schwarzchild_radius/1.0e3);
  shtab->add_constant("c_km_s",o2scl_mks::speed_of_light/1.0e3);
  shtab->add_constant("ell",1.0);
  shtab->add_constant("TMeV",0.0258);

  bool freq_lint=false;
  if (freq_lint) {
    shtab->set_interp_type(itp_linear);
  }
  
  // Add radius and mass information from TOV solver
  shtab->new_column("r");
  shtab->set_unit("r","km");
  shtab->new_column("gm");
  shtab->set_unit("gm","Msun");
  for(size_t i=0;i<shtab->get_nlines();i++) {
    shtab->set("r",i,mr->interp("pr",shtab->get("pr",i),"r"));
    shtab->set("gm",i,mr->interp("pr",shtab->get("pr",i),"gm"));
  }

  shtab->add_constant("Rc",shtab->min("r"));
  shtab->add_constant("R",shtab->max("r"));
  shtab->add_constant("M",freq_mass);

  // Inter-ionic spacing
  shtab->function_column("(3/4/pi/nnuc)^(1/3)","a");
  shtab->set_unit("a","fm");

  // Coulomb Gamma parameter
  shtab->functions_columns("Gamma=Z^2/a/TMeV*hc*alpha");
  shtab->add_constant("R_solid",shtab->interp("Gamma",173.0,"r"));

  // Shear modulus, in fm^-4 and g/cm^3
  shtab->function_column("0.1194*nnuc*Z*Z*alpha/a/(1+0.595*(173/Gamma)^2)",
			 "mu");
  //  shtab->function_column("0.12*nnuc*Z*Z*alpha/a","mu");
  shtab->set_unit("mu","1/fm^4");
  shtab->function_column("mu","mu_cgs");
  for(size_t i=0;i<shtab->get_nlines();i++) {
    shtab->set("mu_cgs",i,shtab->get("mu_cgs",i)*3.5176723e14);
  }
  shtab->set_unit("mu_cgs","g/cm^3");

  // Mass density in 1/fm^4
  shtab->function_column("rho","rho_fm4");
  for(size_t i=0;i<shtab->get_nlines();i++) {
    shtab->set("rho_fm4",i,shtab->get("rho_fm4",i)*2.842789e-15);
  }
  shtab->set_unit("rho_fm4","1/fm^4");

  // Mass density of nuclei in g/cm^3
  shtab->function_column("nnuc*(939+BEoA)","rho_nuc");
  for(size_t i=0;i<shtab->get_nlines();i++) {
    shtab->set("rho_nuc",i,shtab->get("rho_nuc",i)*1.7826618e12);
  }
  shtab->set_unit("rho_nuc","g/cm^3");

  // Mass density of nuclei in 1/fm^4
  shtab->function_column("nnuc*(939+BEoA)/hc","rho_nuc_fm4");
  shtab->set_unit("rho_nuc_fm4","1/fm^4");
    
  for(size_t i=0;i<shtab->get_nlines();i++) {
    shtab->set("rho_nuc",i,shtab->get("rho",i));
    shtab->set("rho_nuc_fm4",i,shtab->get("rho_fm4",i));
  }

  // GR metric functions (unitless)
  shtab->function_column("1/(1-Rschwarz*M/R)","e2l");
  shtab->function_column("1/e2l","e2n");

  // Shear and Alfven speeds (unitless)
  shtab->function_column("sqrt(mu_cgs/rho_nuc)","vs");
  shtab->function_column("eB/sqrt(4.0*pi*rho_nuc_fm4)","vA");
  
  // If the Alfven speed is greater than c, set it equal to 'c'
  for(size_t i=0;i<shtab->get_nlines();i++) {
    if (shtab->get("vA",i)>1.0) shtab->set("vA",i,1.0);
  }

  if (verbose>1) {
    cout << "Computing derivative of shear modulus." << endl;
  }

  // Sort in increasing radius
  shtab->sort_table("r");

  // Derivative of shear modulus
  shtab->deriv("r","mu_cgs","mup_cgs");
  shtab->set_unit("mup_cgs","g/cm^3/km");

  // First derivative coefficient
  shtab->function_column("mup_cgs/rho/(vs^2+vA^2)","Asw");
  shtab->function_column("log(r^4*sqrt(e2n/e2l)*(fr+pr)*vs*vs/1.0e18^4)",
			 "Xdsb");
  shtab->deriv("r","Xdsb","Xdsbp");
  shtab->function_column("Xdsbp*vs^2/(vA^2+vs^2)","Adsb");
  // Units for first derivative coefficient
  shtab->set_unit("Asw","1/km");
  shtab->set_unit("Adsb","1/km");
    
  if (true) {
    cout << "R : " << shtab->get_constant("R") << endl;
    cout << "R_solid: " << shtab->get_constant("R_solid") << endl;
    cout << "Rc: " << shtab->get_constant("Rc") << endl;
  }

  string shear_fn=prefix+"_shear.o2";
  hf.open(shear_fn);
  hdf_output(hf,*shtab,"shear");
  hf.close();
  if (verbose>0) {
    cout << "Wrote shear properties to file '" << shear_fn << "'." << endl;
  }
    
  double freq1=150.0, freq0=1.0e1;
  shear_eigen se;
  se.freq_gr=freq_gr;

  if (false) {
    se.compute_freq(shtab,shtab->get_constant("Rc"),
		    shtab->get_constant("R_solid"),1.0,freq1,1.0e3,1);
    se.store_table(prefix+((string)"_xitab1.o2"));
    shtab->add_constant("freq1",freq1);

    string freq_fn=prefix+"_freq1.o2";
    hf.open(freq_fn);
    hdf_output(hf,*shtab,"freq");
    hf.close();
  }
  
  se.compute_freq(shtab,shtab->get_constant("Rc"),
		  shtab->get_constant("R_solid"),2.0,freq0,1.0e2,1);
  se.store_table(prefix+((string)"_xitab0.o2"));
  shtab->add_constant("freq0",freq0);
  
  string freq_fn=prefix+"_freq0.o2";
  hf.open(freq_fn);
  hdf_output(hf,*shtab,"freq");
  hf.close();
  
#endif
  
  return 0;
}

