/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2018, Andrew W. Steiner
  
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
#ifndef SNA_THERMO_H
#define SNA_THERMO_H

#include <o2scl/fermion_mag_zerot.h>

#include "matter.h"
#include "ldrop_crust.h"
#include "dist_thermo.h"

namespace crust {

  /** \brief Thermodynamics in the single-nucleus approximation
   */
  class sna_thermo {

  protected:

    bool inc_nuc_trans;
    
    /// Electron thermodynamics
    o2scl::fermion_rel relf;

    /// Classical thermodynamics
    o2scl::classical cla;

    /// Electron in a magnetic field
    o2scl::fermion_mag_zerot mfz;

    /// Electron 
    o2scl::fermion elec_B;

    /// Convert units
    o2scl::convert_units &cng;

    /// The mass model
    ldrop_crust *lda;
    
    /** \brief Compute the total fractional volume occupied by the
	the neutrons in all nuclei 

	This function computes
	\f[
	\phi \equiv \sum_i \phi_i \equiv \sum_i 
	\frac{4}{3} \pi R_{n,i}^3 n_i 
	\f]
	
	Used in rxns and in crust(), but can probably be called directly
	from gibbs_fixp_neutron() and/or gibbs_fixp().
    */
    double get_phi_sna(matter &m, double T);

    /** \brief Compute total baryon number density (with dripped neutrons)
     */
    int baryon_density_sna(matter &m, double T);

  public:

    /// \name Different parts of the free energy
    //@{
    double part1, part2, part3, part4;
    //@}
    
    /** \brief Create with the specified mass model and EOS
     */
  sna_thermo(ldrop_crust &lc, o2scl::eos_had_temp_base &he)
    : cng(o2scl::o2scl_settings.get_convert_units()) {
      
      het=&he;
      lda=&lc;

      elec_B.init(cng.convert("kg","1/fm",o2scl_mks::mass_electron),2.0);
      elec_B.inc_rest_mass=true;
      elec_B.non_interacting=true;

      inc_nuc_trans=false;
    }

    /// The magnetic field in Gauss
    double mag_field;

    /** \brief Compute the free energy density in fm^{-4} in the single
	nucleus approximation

	This is faster than free_energy_dist() for one species since
	no iteration is required.
    */
    int free_energy_sna(matter &m, double T);

    /** \brief For fixed nb and nn, determine nnuc and compute the
	free energy

	Used in \ref crust_driver::compute_sna().
    */
    int free_energy_sna_fix_nb_nn(double nb, matter &m, double T);

    /** \brief For fixed nb and nnuc, determine nn and compute the
	free energy

	Used by \ref sna_thermo::free_press_sna_fixed_ye .
    */
    int free_energy_sna_fix_nb_nnuc(double nb, matter &m, double T);

    /** \brief Check \ref free_energy_sna() by computing the
	free energy a different way
    */
    int check_free_energy_sna(dist_thermo &dt, o2scl::test_mgr &t);

    /** \brief The hadronic eos

	This needs to be public so it can be changed by \ref
	crust_driver::model() .
    */
    o2scl::eos_had_temp_base *het;

    /** \brief A base class for function objects
     */
    class snat_funct_base {

    public:

      /// Object of type \ref sna_thermo
      sna_thermo *stp;
      /// Matter object
      matter *mp;
      /// Temperature
      double T_;
      /// Create with specified
      snat_funct_base(sna_thermo &st, matter &m, double &T) {
	stp=&st;
	mp=&m;
	T_=T;
      }
    };

    /** \brief Class to minimize the single-nucleus free energy
	at fixed baryon density as a function of the neutron 
	number density  
	
	Used in \ref crust_driver::compute_sna().
    */
    class free_energy_sna_neut : public snat_funct_base {
      
    public:
      
    free_energy_sna_neut(sna_thermo &st, matter &m, double T,
			 double nb) : snat_funct_base(st,m,T) {
	nb_=nb;
      }
      
      /// Baryon number density
      double nb_;
      
      /** \brief Free energy as a function of neutron drip density
       */
      double operator()(double nn) const {
	mp->n->n=nn;
	stp->free_energy_sna_fix_nb_nn(nb_,*mp,T_);
	return mp->fr;
      }
      
    };

    /** \brief Return the free energy as a one-dimensional function
	of the neutron number density
      
	Used in \ref crust::crust_driver::compute_sna2().
    */
    class solve_nn_ni_fixed_nb_nnhat : public snat_funct_base {
    
    public:
      
      typedef boost::numeric::ublas::vector<double> ubvector;

    solve_nn_ni_fixed_nb_nnhat(sna_thermo &st, matter &m, double &T) : 
      snat_funct_base(st,m,T) {
      }
      
      /** \brief Specify the member function pointer
       */
      void set(sna_thermo &st, double nb, double nnhat, matter &m, 
	       double &T) {
	stp=&st;
	mp=&m;
	T_=T;
	nnhat_=nnhat;
	nb_=nb;
      }
      
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual int operator()(size_t nv, const ubvector &x, 
			     ubvector &y) {
	double y0, y1;
	mp->n->n=x[0];
	mp->dist[0].n=x[1];

	stp->free_energy_sna(*mp,T_);
	stp->baryon_density_sna(*mp,T_);
	
	double nnhat=mp->n->n*(1.0-stp->get_phi_sna(*mp,T_));
	
	y[0]=(nb_-mp->nb)/nb_;
	y[1]=(nnhat_-nnhat)/nnhat_;

	return 0;
      }
    
    protected:
    
      /// \name The parameters
      //@{
      double nb_, nnhat_;
      //@}
    
    };

    /** \brief To test minimization with respect to \f$ {\hat n}_n \f$
	instead of \f$ n_n \f$. 
	
	Used at the end of \ref crust::crust_driver::compute_sna().
    */
    class free_energy_sna_fix_nb_nnhat : snat_funct_base {
  
    public:
  
      typedef boost::numeric::ublas::vector<double> ubvector;

      /** \brief Specify the member function pointer
       */
      free_energy_sna_fix_nb_nnhat
	(sna_thermo &st, double nb, matter &m, double &T) : 
      snat_funct_base(st,m,T), mfna(st,m,T) {
      
	nb_=nb;
	x.resize(2);
	x[0]=m.n->n;
	x[1]=m.dist[0].n;

	gmh.tol_abs/=1.0e4;
	gmh.tol_rel/=1.0e4;
      }
  
      virtual ~free_energy_sna_fix_nb_nnhat() {};
  
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual double operator()(double nnhat) {
	mfna.set(*stp,nb_,nnhat,*mp,T_);
	o2scl::mm_funct fmp=std::bind
	  (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
	   (&solve_nn_ni_fixed_nb_nnhat::operator()),&mfna,
	   std::placeholders::_1,std::placeholders::_2,
	   std::placeholders::_3);
	gmh.msolve(2,x,fmp);
	stp->free_energy_sna(*mp,T_);
	return mp->fr;
      }

    protected:

      /// The solver
      o2scl::mroot_hybrids<> gmh;

      /// Hold the solution for the neutron and nuclear densities
      ubvector x;

      /// The object to solve
      solve_nn_ni_fixed_nb_nnhat mfna;
  
      /// Baryon density
      double nb_;

    };

    /** \brief Class to give the properties of matter as a function
	of baryon density at fixed electron fraction in the 
	single nucleus approximation
    */
    class free_press_sna_fixed_ye : public snat_funct_base {
      
    public:
      
    free_press_sna_fixed_ye(sna_thermo &st, matter &m, double T, 
			    dist_thermo &dt, bool return_fr) : 
      snat_funct_base(st,m,T) {
	dtp=&dt;
	// Evaluate current electron fraction
	st.baryon_density_sna(m,T);
	Ye_=m.dist[0].Z*m.dist[0].n/m.nb;
	return_fr_=return_fr;
      }
      
      /// Dist thermo object for gibbs energy
      dist_thermo *dtp;

      /// Electron fraction
      double Ye_;

      /// If true, return free energy. Otherwise, return Gibbs energy
      bool return_fr_;

      /** \brief Free energy as a function of neutron drip density
       */
      double operator()(double nb) const {

	// Vary nuclear density ensuring constant electron fraction
	mp->dist[0].n=nb*Ye_/mp->dist[0].Z;

	/*
	// Double check that the electron fraction is constant
	std::cout.precision(8);
	std::cout << "gesy1: " << mp->dist[0].n << " " 
	<< mp->dist[0].Z << " " << nb << " " 
	<< mp->dist[0].n*mp->dist[0].Z/nb << " " << Ye_ << std::endl;
	*/
	
	stp->free_energy_sna_fix_nb_nnuc(nb,*mp,T_);
	
	/*
	  std::cout << "gesy2: " << mp->dist[0].n << " " 
	  << mp->dist[0].Z << " " << nb << " " 
	  << mp->dist[0].n*mp->dist[0].Z/nb << " " << Ye_ << std::endl;
	*/

	if (return_fr_) return mp->fr;
	// Now compute the pressure
	dtp->gibbs_energy_dist(*mp,T_);
	return mp->pr;
      }
      
    };

  };

}

#endif
