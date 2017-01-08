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
#ifndef DIST_THERMO_H
#define DIST_THERMO_H

#include <o2scl/fermion_mag_zerot.h>
#include <o2scl/classical.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_nonrel.h>

#include "matter.h"
#include "ldrop_crust.h"

#ifndef DOXYGENP
namespace crust {
#endif

  /** \brief An object to compute the thermodynamic properties
      of a distribution

      The variable \f$ \chi \f$ is the fractional volume occupied by
      the neutrons relative to the rest of the Wigner-Seitz cell
      \f[
      \chi_i \equiv \left( \frac{R_{n,i}}{R_{\mathrm{WS},i}} \right)^3
      \f]

      \future It would be nice to simplify free_energy_dist(),
      free_energy_sna(), and other related functions so there
      isn't so much code repetition.
  */
  class dist_thermo {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:

    /// For the electron thermodynamics
    o2scl::fermion_rel relf;
    /// For neutron and protons
    o2scl::fermion_nonrel nrelf;
    /// For the thermodynamics of the nuclei
    o2scl::classical cla;
    /// Thermodynamics in the magnetic field
    o2scl::fermion_mag_zerot mfz;
    /// Electron in a magnetic field
    o2scl::fermion elec_B;

    /** \brief Function to solve for a fixed baryon density
     */
    double dist_sna_solve(double nb, matter &m, double T);

    /** \brief Function to solve to ensure the pressure is equal to \c
	pr0

	The variable \c n_fact is the factor that one should
	multiply all the number densities of the nuclei by
	in order to ensure that the new pressure is equal
	to \c pr0.

	This function returns <tt>(pr-pr0)/pr0</tt>

	Used in \ref gibbs_fixp().
    */
    double gibbs_pressure_func(double pr0, double n_fact, matter &m_new, 
			       double T);

    /** \brief Function to solve to ensure the pressure is equal to \c pr
	and that the "full" neutron density is \c nnf
    */
    int gibbs_press_neut_func(double pr, double nnf, double alpha,
			      double nn_new, matter &m_new, double T,
			      double &y1, double &y2);

    /** \brief A base class for function objects
     */
    class dt_funct_base {
    public:
      /// The base pointer
      dist_thermo *dtp;
      /// Matter data
      matter *mp;
      /// Temperature
      double Tv;
      /// Create a function object
      dt_funct_base(dist_thermo &dt, matter &m, double &T) {
	dtp=&dt;
	mp=&m;
	Tv=T;
      }
    };
    
    /// Convert units
    o2scl::convert_units *cng;

    /** \brief Solve for the pressure
	
	Used in crust_driver::gibbs_energy_pressure()
    */
    class funct_solve_pressure : public dt_funct_base {
      
    public:
  
      /** \brief Specify the member function pointer
       */
      funct_solve_pressure(dist_thermo &dt, matter &m, double &T,
			   double pr0) : dt_funct_base(dt,m,T) {
	ppr0=pr0;
      }
  
      virtual ~funct_solve_pressure() {};
  
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual double operator()(double alf) const {
	return dtp->gibbs_pressure_func(ppr0,alf,*mp,Tv);
      }

#ifndef DOXYGEN_INTERNAL
    
    protected:
  
      /// The pressure to fix
      double ppr0;

#endif

    };

    /** \brief Solve for the pressure
    
	Used in crust_driver::gibbs_energy_pressure_neutron()
    */
    class funct_solve_pressure2 : public o2scl::mm_funct11, dt_funct_base {
  
    public:
  
      /** \brief Specify the member function pointer
       */
    funct_solve_pressure2(dist_thermo &dt, matter &m, double T,
			  double pr0, double nn_full) : 
      dt_funct_base(dt,m,T) {
	
	ppr0=pr0;
	nnfull=nn_full;
      }
      
      virtual ~funct_solve_pressure2() {};
      
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual int operator()(size_t nv, const ubvector &x, 
			     ubvector &y) {
	// It's important not to return zero here because the 
	// function gibbs_pressure_func2 can return an error code 
	return dtp->gibbs_press_neut_func
	  (ppr0,nnfull,x[0],x[1],*mp,Tv,y[0],y[1]);
      }

#ifndef DOXYGEN_INTERNAL
    
    protected:
  
      /// \name The parameters
      //@{
      double ppr0, nnfull;
      //@}

#endif

    };

    /** \brief Solve for the pressure
    
	Used in crust_driver::gibbs_energy_pressure_neutron()
    */
    class funct_solve_pressure3 : public o2scl::mm_funct11, dt_funct_base {
  
    public:
  
      /** \brief Specify the member function pointer
       */
    funct_solve_pressure3(dist_thermo &dt, matter &m, double T,
			  double pr0, double nn_full) : 
      dt_funct_base(dt,m,T) {
	
	ppr0=pr0;
	nnfull=nn_full;
      }
      
      virtual ~funct_solve_pressure3() {};
      
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual int operator()(size_t nv, const ubvector &x, 
			     ubvector &y) {
	// It's important not to return zero here because the 
	// function gibbs_pressure_func2 can return an error code 
	return dtp->gibbs_press_neut_func
	  (ppr0,nnfull,x[0],exp(x[1]),*mp,Tv,y[0],y[1]);
      }

#ifndef DOXYGEN_INTERNAL
    
    protected:
  
      /// \name The parameters
      //@{
      double ppr0, nnfull;
      //@}

#endif

    };

  public:

    /** \brief Create with the specified unit conversion class
    */
    dist_thermo(o2scl::convert_units &conv) {
      elec_B.init(o2scl::o2scl_settings.get_convert_units().convert
		  ("kg","1/fm",o2scl_mks::mass_electron),2.0);
      elec_B.inc_rest_mass=true;
      elec_B.non_interacting=true;

      cng=&conv;

      check=0;
    }
    
    /// \name Different parts of the free energy
    //@{
    double part1, part2, part3, part4;
    //@}
    
    /// The mass model
    ldrop_crust *lda;
    
    /// The hadronic eos (default points to \c sk)
    o2scl::eos_had_temp_base *het;
    
    /// The magnetic field in Gauss
    double mag_field;

    /// Check the pressure computed in gibbs_energy_dist()
    int check_pressure(matter &m, double T);

    /// Given nb and nn, determine nnuc and compute the free energy
    int free_energy_dist_sna(double nb, matter &m, double T);

    /** \brief For fixed \f$ \alpha=n_n (1-\phi) \f$, solve for n_n
	and compute the free energy 

	This function uses iteration. 
    */
    int free_energy_dist_alpha(matter &m, double &T, double alpha);

    /** \brief Compute the free energy density in 
	\f$ \mathrm{fm}^{-4} \f$

	Also computes the baryon density and the rest mass
	energy density.

	If \c eval_chem_pots is true, this function evaluates the
	chemical potentials and stores them in \ref matter::mu and
	\ref matter::mun.

	\future Offer a direct method of computing derivatives
	from an o2scl deriv object.
    */
    int free_energy_dist(matter &m, double T, bool eval_chem_pots=false);

    /** \brief Compute the Gibbs energy density from the free energy density
	
	Also computes the free energy, pressure, rest mass energy
	density, and baryon density as a by-product.
    */
    int gibbs_energy_dist(matter &m, double T);

    /** \brief Compute the rest mass density in \f$
	\mathrm{g}/\mathrm{cm}^3 \f$ (including the nuclear binding
	energy) at temperature \c T

	Results are stored in \ref matter::rho .
    */
    int mass_density(matter &m, double T);

    /** \brief Compute the total free energy of one cell
	
	Compute the free energy of cell of species \c ix in 
	\f$ \mathrm{fm}^{-1} \f$, the volume in 
	\f$ \mathrm{fm}^{3} \f$, and the value of 
	\f$ N^{\prime} \equiv N - \phi n_{n,\mathrm{out}} / n_i \f$.
    */
    int free_energy_cell(matter &m, double T, int ix, double &fr, 
			 double &vol, double &Nprime);

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
    double get_phi(matter &m, double T);

    /** \brief Compute total baryon number density (with dripped neutrons)
     */
    int baryon_density(matter &m, double T);

    /** \brief Create a new configuration \c m_new with a fixed
	pressure of \c pr_target at temperature \c T 

	This function modifies the composition in \c m_new by
	multiplying all of the nuclear number densities and the
	external neutron density (in <tt>m.n.n</tt>) by the same
	factor until the pressure at temperature \c T is equal to that
	specified in <tt>pr_target</tt>.

	Because this function operates on <tt>m.n.n</tt> and not on
	the quantity <tt>m.n.n*(1-phi)</tt>, it is typically only
	used with <tt>m.n.n</tt> is zero. 
	
	Uses gibbs_pressure_func().
	
	Currently, this function uses a local \ref
	cern_mroot_root object to solve for pressure equality,
	but there is some legacy code for using \ref gsl_root_brent to
	do the solution.
    */
    int gibbs_fixp(double pr_target, matter &m_new, double T);

    /** \brief Create a new configuration \c m_new with pressure 
	of \c pr_target at temperature T and a fixed value
	of \f$ n_n (1-\phi) \f$ .

	Given the old neutron density and value of phi in \c nn_old
	and \c phi_old, this modifies the composition of \c m_new
	until its pressure at temperature \c T and at fixed \f$ \alpha
	= n_n (1-\phi) + n_{n,\mathrm{shift}} \f$ is the same as that
	given in \c pr_target.

	This uses the function gibbs_press_neut_func which
	has two unknowns, the new external neutron number density,
	<tt>m_new.n.n</tt> and the correction factor \c corr 
	to multiply the nuclear densities by. It solves for
	pressure equality, and for
	\f[
	(1-\phi^{[2]})n^{[2]}_{n} = \mathrm{corr} (1-\phi^{[1]})n^{[1]}_{n}
	\f]
	which ensures that the reduced number of external neutrons
	is properly scaled by the same factor.
    */
    int gibbs_fixp_neutron(double pr_old, double nn_full_old,
			   matter &m_new, double T, double nn_shift);

    /** \brief Compute the gibbs energy per baryon of cell with index
	\c ix at pressure \c P
	
	Uses free_energy_cell() to compute the free energy, etc.
    */
    int gibbs_energy_per_baryon_cell(matter &m, double T, double P, int ix,
				     double &gpb);

    /** \brief The total Gibbs energy of one cell
     */
    int gibbs_energy_cell(matter &m, double T, double P, int ix,
			  double &g_cell);

    /// \name Testing
    //@{
    int check;
    static const int check_none=0;
    static const int check_mass_density=1;
    static const int check_free_energy_sna=2;
    static const int check_free_energy_dist=3;
    static const int check_ldrop_derivs=5;
    static const int check_free_energy_cell=8;
    //@}

    /** \brief Return the free energy as a one-dimensional function
	of the neutron number density

	Used in \ref crust_driver::compute_sna_dist().
    */
    class free_energy_dist_neut : public dt_funct_base {

    protected:

      /// The baryon density
      double nbx;
      
    public:
      
      /** \brief Specify the member function pointer
       */
    free_energy_dist_neut(dist_thermo &dt, matter &m, double &T,
			  double nb) : dt_funct_base(dt,m,T) {
	nbx=nb;
      }
      
      virtual ~free_energy_dist_neut() {};
  
      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual double operator()(double nn) const {
	mp->n->n=nn;
	dtp->free_energy_dist_sna(nbx,*mp,Tv);
	return mp->fr;
      }

    };

    /** \brief Return the free energy as a one-dimensional function
	of the nuclear number density
    */
    class funct_dist_sna : public dt_funct_base {
	
    public:
	
      /** \brief Specify the member function pointer
       */
    funct_dist_sna(dist_thermo &dt, matter &m, double T,
		   double nb) : dt_funct_base(dt,m,T) {
	nb_=nb;
      }
	
      virtual ~funct_dist_sna() {};

      /** \brief Compute the function at point \c x, with result \c y
       */
      virtual double operator()(double ni) const {
	mp->dist[0].n=ni;
	return dtp->dist_sna_solve(nb_,*mp,Tv);
      }
	
#ifndef DOXYGEN_INTERNAL
	
    protected:
  
      /// Baryon density
      double nb_;

#endif

    };
    
    /// For computing the derivative of the free energy wrt alpha
    class free_dist_deriv_alpha {
    public:
      free_dist_deriv_alpha(dist_thermo &dt, matter &m, 
			    double &T, size_t pt=0) {
	dtp=&dt;
	m_=&m;
	T_=T;
	pt_=pt;
      }
      virtual ~free_dist_deriv_alpha() {}

      virtual double operator()(double x) {
	double ret;
	dtp->free_energy_dist_alpha(*m_,T_,x);
	ret=m_->fr;
	if (pt_==1) return dtp->part1;
	if (pt_==2) return dtp->part2;
	if (pt_==3) return dtp->part3;
	if (pt_==4) return dtp->part4;
	return ret;
      }
      dist_thermo *dtp;
      matter *m_;
      double T_;
      size_t pt_;
    };

    /** \brief For computing the derivative of the free energy wrt 
	the number density of one of the nuclei
    */
    class free_dist_deriv_nuc {
    public:
      free_dist_deriv_nuc(dist_thermo &dt, matter &m, double &T, double alpha, 
			  size_t i, size_t pt=0) {
	dtp=&dt;
	m_=&m;
	T_=T;
	i_=i;
	alpha_=alpha;
	pt_=pt;
      }
      virtual ~free_dist_deriv_nuc() {}
      virtual double operator()(double x) {
	double ret;
	m_->dist[i_].n=x;
	dtp->free_energy_dist_alpha(*m_,T_,alpha_);
	ret=m_->fr;
	if (pt_==1) return dtp->part1;
	if (pt_==2) return dtp->part2;
	if (pt_==3) return dtp->part3;
	if (pt_==4) return dtp->part4;
	return ret;
      }
      size_t i_;
      dist_thermo *dtp;
      matter *m_;
      double T_, alpha_;
      size_t pt_;
    };

  };

#ifndef DOXYGENP
}
#endif

#endif
