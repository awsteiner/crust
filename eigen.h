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

namespace crust {

  /** \brief Compute eigenfrequencies from a crust table

      <b>Steiner & Watts</b>

      \f[
      \xi^{\prime \prime} + A_{\mathrm{sw}} \xi^{\prime} + 
      B_{\mathrm{sw}} = 0
      \f]
      \f[
      A_{\mathrm{sw}} \equiv \frac{\mu^{\prime}}{\rho} \frac{1}{v_s^2+v_A^2}
      \f]
      \f[
      B_{1,\mathrm{sw}} \equiv \omega^2 \left( 1+v_A^2\right)
      \f]
      \f[
      B_{2,\mathrm{sw}} \equiv \frac{-(\ell-1)(\ell+2) v_s^2}{R^2}
      \f]
      \f[
      B_{\mathrm{sw}} \equiv \frac{B_{1,\mathrm{sw}}+B_{2,\mathrm{sw}}}
      {v_s^2+v_A^2}
      \f]

      <b>Term from Samuelsson and Andersson</b>
      
      \f[
      A_{\mathrm{sa}} \equiv \frac{d}{dr} \left\{ 
      \log \left[ r^4 e^{\nu-\lambda} 
      (\varepsilon + p_t) v_s^2 \right]
      \right\}
      \f]
      \f[
      B_{1,\mathrm{sa}} \equiv e^{2 \lambda-2\nu} \omega^2
      \f]
      \f[
      B_{2,\mathrm{sa}} \equiv -e^{2 \lambda} v_s^2
      \frac{\left(\ell-1\right)\left(\ell+2\right)}{r^2}
      \f]
      \f[
      B_{\mathrm{sa}} \equiv \frac{B_{1,\mathrm{sa}}+B_{2,\mathrm{sa}}}{v_s^2} 
      \f]

      <b>New expressions for Diebel, Steiner, and Brown</b>

      \f[
      \xi^{\prime \prime} + A_{\mathrm{dsb}} \xi^{\prime} + 
      B_{\mathrm{dsb}} = 0
      \f]
      \f[
      A_{\mathrm{dsb}} \equiv 
      \frac{v_s^2}{v_s^2+v_A^2} \left \{ 
      \log \left[ r^4 e^{\nu-\lambda} 
      \left(\varepsilon+p_t\right) v_s^2 \right]
      \right\}^{\prime}
      \f]
      \f[
      B_{1,\mathrm{dsb}} \equiv \omega^2 e^{2 \lambda-2 \nu} 
      (1+v_A^2)
      \f]
      \f[
      B_{2,\mathrm{dsb}} \equiv \frac{-(\ell-1)(\ell+2) v_s^2 e^{2 \lambda}}
      {r^2}
      \f]
      \f[
      B_{\mathrm{dsb}} \equiv \frac{B_{1,\mathrm{dsb}}+B_{2,\mathrm{dsb}}}
      {v_s^2+v_A^2}
      \f]

   */
  class shear_eigen {
  
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    shear_eigen();

    /// If true, include GR effects in frequency calculation (default false)
    bool freq_gr;

    /// Pointer to shear table
    o2scl::table_units<> *shtab;
  
    /// ODE solver
    o2scl::ode_iv_solve<> ois;
    
    /// Grid of radii
    ubvector r_grid; 

    /// Dispacement 
    ubmatrix xi;

    /// Dispacement at r=R_c
    ubvector xi_start;

    /// Derivatives
    ubmatrix dxi;

    /// Displacement uncertainty
    ubmatrix xi_err;

    /// Radius
    double Rad;

    /// Radius of bottom edge of crust
    double Radc;

    /// Table size 
    size_t ntab;

    /// Store output table
    void store_table(std::string fname);

    /// The number of nodes
    int nodes;

    /// Shear frequency differential equation
    int derivs(double r, size_t nv, const ubvector &y, ubvector &dydx);

    /// Solve for the eigenvalue problem
    int eigen_solve(size_t nv, const ubvector &x, ubvector &y); 

    /// Compute frequency
    int compute_freq(o2scl::table_units<> *shear, double Rc, double R, 
		     double ell, double &freq_lo, double freq_hi, 
		     int verbose=0);

    /// Interpolation
    o2scl::interp<ubvector> si;
  };

  /** \brief Compute neutron star and shear modulus
   */
  class tov_shear {

  public:

    /// Verbosity parameter
    int verbose;

    /// The magnetic field in Gauss
    double mag_field;

    /// Desc
    std::string in_dir;

    tov_shear() {
      verbose=1;
      hd_flag=0;
      freq_gr=false;
      freq_mass=1.4;
      mag_field=0.0;

      cng.units_cmd_string="units -f indata/units_hck.dat";
      in_dir="indata";
    }

    /// Convert units
    o2scl::convert_units<double> cng;

    /// If true, include GR effects in frequency calculation (default false)
    bool freq_gr;

    /// Mass for frequency calculation (default 1.4)
    double freq_mass;

    /** \brief High-density EOS flag (default 0)

        Used in \ref tov().

        0 is most probable, 1 is plus 1 sigma, 2 is plus 2 sigma,
	3 is minus 1 sigma, and 4 in minus 2 sigma
    */
    int hd_flag;
    
    /// Compute crust in full_equilibrium with a distribution
    int tov(std::vector<std::string> &sv, bool itive_com);
    
    /// Compute shear properties and magnetar oscillation frequencies
    int shear(std::vector<std::string> &sv, bool itive_com);


  };

}
