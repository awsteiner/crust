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
#include <iostream>
#include <vector>

#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/table.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/cli_readline.h>
#include <o2scl/convert_units.h>
#include <o2scl/min_cern.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/eos_tov.h>
#include <o2scl/tov_solve.h>
#include <o2scl/deriv_cern.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/nucmass_ldrop_shell.h>
#include <o2scl/hdf_io.h>

#include "dist_thermo.h"

namespace crust {

  /** \brief Compute pycnonuclear reaction rates

      The S-factor here is a function of \f$ Z_1, Z_2, A_1, A_2 \f$
      and the energy in the center of mass, \f$ E \f$.

      Ensure Z's and A's are even by increasing Z and decreasing A.
      Rates tend to decrease with Z and increase with A, so this
      choice minimizes fusion.
      
      \ref crust_Beard10 "Beard10" gives computes the S-factors for several
      reactions and fits each to the same functional form
      \f[
      S(E) = \exp \left\{ B_1 + B_2 E + B_3 E^2 + 
      \frac{C_1 + C_2 E + C_3 E^2 + C_4 E^3}
      {1 + \exp\left[\left(E_C - E\right)/D\right]}\right\}
      \f]
      these rates are tabulated for reactions among even-even isotopes
      of C, O, Ne, and Mg up to a maximum value of A (\f$
      A_{\mathrm{max}} \f$) which depends on which nuclei are being
      fused. S-factors for \f$ A > A_{\mathrm{max}} \f$ are assumed to
      be equal to those for \f$ A > A_{\mathrm{max}} \f$ (which
      increases fusion slightly for heavier isotopes).

      \ref crust_Yakovlev10 "Yakovlev10" 
      gives a fit in Eqs. 5 and 6 in terms of four
      quantities, \f$ S_0, E_C, \delta \f$ and \f$ \xi \f$:
      \f[
      S(E) = 
      S_0 \exp \left[ 2 \pi \eta + \Phi(E) \right]
      \f]
      for \f$ E \leq E_C \f$ and 
      \f[
      S(E) = S_0 \exp (2 \pi \eta) \sqrt{E/E_C} \left[ 
      1+\xi (E-E_C)/E\right]
      \f] 
      for \f$ E > E_C \f$
      where \f$ \eta \equiv \sqrt(E_R/E) \f$, \f$ E_R \equiv 
      \alpha^2 \mu /( 2 \hbar^2) \f$, \f$ \mu \f$ is the reduced 
      mass, \f$ \alpha \equiv Z_1 Z_2 e^2 \f$, and 
      \f$ \Phi(E) \equiv \f$ ...

      The fit and Eqs. 19 and 20 given in Table II of gives quantities
      \f{eqnarray*}
      E_C &=& \alpha/R_C^{(0)} \nonumber \\
      R_C^{(0)} &=& R + \Delta R_1 |A_1 - A_{10}| + \Delta R_2 
      | A_2 - A_{20}| \nonumber \\
      \xi &=& \xi_0 + \xi_1 (A_1 + A_2) \nonumber
      \f}
      where the table gives \f$ R, \Delta R_{1 a}, \Delta R_{2 a}, 
      \Delta R_{1 b}, \Delta R_{2 b}, \delta, S_0, \xi_0 \f$, and
      \f$ \xi_1 \f$  where \f$ A_{10} \equiv 2 Z_1 \f$ and 
      \f$ A_{20} \equiv 2 Z_2 \f$ and the quantities \f$ \Delta R_1 \f$
      and \f$ \Delta R_2 \f$ are given by
      \f{eqnarray*}
      \Delta R_1 &=& \theta(A_1-A_{10}) \Delta R_{1a} +
      \theta(A_{10}-A_1) \Delta R_{1b} \nonumber \\
      \Delta R_2 &=& \theta(A_2-A_{20}) \Delta R_{2a} +
      \theta(A_{20}-A_2) \Delta R_{2b} \nonumber
      \f}

      The default values of \ref Cpyc, \ref Cexp, and \ref Cpl
      are from the "optimal" model in Table II in \ref crust_Yakovlev06
      "Yakovlev06" .
  */
  class pyc_rates {

  public:

    bool loaded;

    /// Overall coefficient (default 3.9)
    double Cpyc;
    /// Coefficient in exponent (default 2.638)
    double Cexp;
    /// Coefficient controlling power of \f$ \lambda \f$ (default 1.25)
    double Cpl;
    /// Table of S-factor coefficients from \ref crust_Beard10 "Beard10"
    o2scl::table<> b10;
  
    /// Convert units
    o2scl::convert_units &cng;
    
    /// If true, use debug mode (default false)
    bool debug;

    /// Allow fusions involving high Z nuclei(default true)
    bool allow_highZ;
  
    /// If true, add fusions from Z=14 nuclei (default false)
    bool add14;

    /** \brief If true, use fit from \ref crust_Yakovlev10 "Yakovlev10" 
	(default false)
    */
    bool use_fit;

    /** \brief Fit coefficients from Table II in 
	\ref crust_Yakovlev10 "Yakovlev10"
    */
    double fit[10][13];
    
    /// Initialize with "Optimal" rate
    pyc_rates();

    /** \brief Load the S-factor data

	This is called by the constructor to set \ref b10.
    */
    int load_data(std::string dir);

    /** \brief Set fit coefficients from Table in Yakovlev et al. 2010

	This is called by the constructor to set \ref fit.
     */
    int set_fit();

    /// Test the S-factor data
    int test_Sfactor();

    /** \brief Compute the fusion times for C, O, Ne, and Mg assuming
	one-component composition
     */
    int fusion_times(dist_thermo &dt, matter &m, double T, 
		     double &t_acc, double &t_c_fusion,
		     double &t_o_fusion, double &t_ne_fusion,
		     double &t_mg_fusion);
    
    /** \brief Compute the S-factor in units of 
	\f$ \mathrm{MeV} \cdot \mathrm{barn} \f$
    */
    double Sfactor(size_t Z1, size_t Z2, size_t A1, size_t A2,
		   double E);

    /** \brief Compute the plasma temperature \f$ T_{12}^{(p)} \f$ in 1/fm

	Density \c n12 should be in \f$ \mathrm{fm}^{-3} \f$ and reduced
	mass \c mu12 should be in 1/fm.

	From Eq. 10 in \ref crust_Yakovlev06 "Yakovlev06"
	this gives the plasma temperature as
	\f[
	T_p = \sqrt{ \frac{4 \pi Z_1 Z_2 e^2 n_{12}}{2 \mu_{12}} }
	\f]
	where \f$ e^2 \f$ is replaced with \ref gsl_num::fine_structure.
    */
    double Tp(size_t Z1, size_t Z2, double n12, double mu12);
  
    /** \brief Compute \f$ r_{B12} \f$ in fm
      
	Reduced mass \c mu12 should be in 1/fm.
	
	From Eq. 11 in \ref crust_Yakovlev06 "Yakovlev06"
	\f[
	r_{B12} = \frac{1}{2 \mu_{12} Z_1 Z_2 e^2}
	\f]
	where \f$ e^2 \f$ is replaced with \ref gsl_num::fine_structure
	and \c m1 and \c m2 are ignored.
    */
    double rB(size_t Z1, size_t Z2, double m1, double m2,
	      double mu12);
  
    /** \brief Compute \f$ \lambda_{12} \f$ and \f$ \omega_{12} \f$
      
	Densities should be in \f$ \mathrm{fm}^{-3} \f$, masses should
	be in 1/fm, lambda is unitless, and omega is returned in 1/fm.

	The quantity \c omega is computed from \ref Tp() \c lambda is
	computed from Eq. 12 in \ref crust_Yakovlev06 "Yakovlev06"
	\f[
	\lambda_{12} = r_{B12} \left( \frac{n12}{2} \right)^{1/3}
	\f]
	where \f$ n_12 \equiv \frac{3}{4 \pi a_{12}^3} \f$, \f$ a_{12}
	\equiv (a_1 + a_2)/2 \f$ , \f$ a_i \equiv [3/(4 \pi
	n_i)]^{1/3} \f$ and \f$ r_{B12} \f$ is given from \ref rB().
    */
    int lambda_omega(size_t Z1, size_t Z2, 
		     double m1, double m2, double n1, double n2,
		     double &lambda, double &omega);

    /** \brief Compute \f$ \lambda_{12} \f$ and \f$ \omega_{12} \f$
      
	Densities should be in \f$ \mathrm{fm}^{-3} \f$, masses should
	be in 1/fm, lambda is unitless, and omega is returned in 1/fm.

	The quantity \c omega is computed from \ref Tp() \c lambda is
	computed from Eq. 12 in \ref crust_Yakovlev06 "Yakovlev06"
	\f[
	\lambda_{12} = \frac{A_1 + A_2}{A_1 A_2 Z_1 Z_2
	\left(Z_1^{1/3} + Z_2^{1/3}\right)}
	\left( \frac{ \rho X_N \left<Z\right> }
	{\left<A\right> \rho_0 } \right)^{1/3}
	\f]
	where \f$ \rho_0 \equiv 1.3574 \times 10^{11} \mathrm{g/cm}^{3} 
	\f$ and \f$ X_N \f$ is the total mass fraction contained in 
	nuclei.

	This function doesn't work currently.
    */
    int lambda_omega2(size_t Z1, size_t Z2, size_t A1, size_t A2,
		      double m1, double m2, double n1, double n2,
		      double rho, double XN, double aveZ, double aveA,
		      double &lambda, double &omega);

    /// Test the reaction rate
    int test_rate();

    /// Test the reaction rate
    int test_rate2();

    /** \brief Return the rate in fm^{-3}*s^{-1}

	Densities n1, n2, and ntot should be in 
	\f$ \mathrm{fm}^{-3} \f$, masses should be in 1/fm, 
	rho should be in g/cm3, 

	As in Eq. 33 in \ref crust_Yakovlev06 "Yakolev06" , we define 
	\f[
	R^{\mathrm{pyc}}_{ij} \equiv 
	8 \times 10^{7} C_{\mathrm{pyc}} 
	\frac{\rho X_N x_i x_j A_i A_j \left<A\right>
	Z_i^2 Z_j^2}{(1+\delta_{ij}) (A_1+A_2)^2}
	S(E) {\tilde{\lambda}}^{3-C_{\mathrm{pl}}}
	\exp \left( -\frac{C_{\mathrm{exp}}}
	{\sqrt{\tilde{\lambda}}}\right)~\mathrm{fm}^{-3}~\mathrm{s}
	\f]
    */
    double rate(size_t Z1, size_t Z2, size_t A1, size_t A2,
		double m1, double m2, double n1, double n2,
		double ntot, double avgZ, double avgA, double rho, 
		double XN);

    /// Test a mixture to see that the reaction rates are reasonable
    int test_mixture();

    /** \brief Return \c true if the rate is allowed and \c false 
	otherwise

	Densities n1, n2, and ntot should be in \f$ \mathrm{fm}^{-3}
	\f$, rho should be in g/cm3, Mdot should be in Msun/yr, mue
	should be in 1/fm

	We always replace with a smaller rate, that is we
	increase Z by one if Z is odd and decrease A by 
	one if A is odd. 
	
	If either of \c Z1 or \c Z2 are less than or equal 
	to 4, this function always returns \c true.

	Using the rate given in \ref rate(), the associated
	fusion timescale is computed as
	\f[
	t_{\mathrm{fusion}} = \frac{n_{\mathrm{tot}}}{3 R}
	\f]
	where \f$ n_{\mathrm{tot}} \f$ is the total number density of
	nuclei.
	
	The accretion timescale is 
	\f[
	t_{\mathrm{acc}} = \frac{4 \pi R^2 \dot{M}}{\Sigma} 
	\f]
	where \f$ \Sigma = P/g\f$ is the local column depth, \f$ P \f$
	is the pressure and \f$ g \f$ is the surface gravity.

	\comment
	sigma = P/(surf. grav.) = (g/cm/s^2) / (cm/s^2) is in units of g/cm^2
	sigma dot = g/cm^2/s 
	\endcomment
	
    */
    bool is_allowed(size_t Z1, size_t Z2, size_t A1, size_t A2,
		    double m1, double m2, double n1, double n2,
		    double ntot, double avgZ, double avgA, double rho,
		    double XN, double Mdot, double P);

    /** \brief Get fusion and accretion times 

	Densities n1, n2, and ntot should be in \f$ \mathrm{fm}^{-3}
	\f$, rho should be in g/cm3, Mdot should be in Msun/yr, mue
	should be in 1/fm

	We always replace with a smaller rate, that is we
	increase Z by one if Z is odd and decrease A by 
	one if A is odd. 
	
	If either of \c Z1 or \c Z2 are less than or equal 
	to 4, this function always returns \c true.

	Using the rate given in \ref rate(), the associated
	fusion timescale is computed as
	\f[
	t_{\mathrm{fusion}} = \frac{n_{\mathrm{tot}}}{3 R}
	\f]
	where \f$ n_{\mathrm{tot}} \f$ is the total number density of
	nuclei.
	
	The accretion timescale is 
	\f[
	t_{\mathrm{acc}} = \frac{4 \pi R^2 \dot{M}}{\Sigma} 
	\f]
	where \f$ \Sigma = P/g\f$ is the local column depth, \f$ P \f$
	is the pressure and \f$ g \f$ is the surface gravity.

     */
    void get_times(size_t Z1, size_t Z2, size_t A1, size_t A2,
		   double m1, double m2, double n1, double n2,
		   double ntot, double avgZ, double avgA, double rho,
		   double Xn, double Mdot, double P, double &t_fusion,
		   double &t_acc);

  };

}
