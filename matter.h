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
#ifndef MATTER_H
#define MATTER_H

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

#include "ldrop_crust.h"

#ifndef DOXYGENP
namespace crust {
#endif

  /** \brief An object which describes the composition and 
      thermodynamics of matter at some specified density or pressure
  */
  class matter {
  
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Quasi-free neutrons
    o2scl::fermion *n;
    /// Quasi-free protons
    o2scl::fermion *p;
    /// Distribution of nuclei
    std::vector<o2scl::nucleus> dist;
    /// Electrons
    o2scl::fermion e;
    /// Temperature
    double T;
    /// Energy density
    double ed;
    /// Pressure
    double pr;
    /// Free energy
    double fr;
    /// Gibbs energy
    double gb;
    /// Mass density
    double rho;
    /// Baryon number density
    double nb;
    /// Neutron chemical potential
    double mun;
    /// Proton chemical potential
    double mup;
    /// Nuclear chemical potentials
    ubvector mu;
    /// Ratio for \f$ d \chi_i / d n_i \f$ 
    ubvector zeta;
    /// Thermodynamic information for dripped particles
    o2scl::thermo drip_th;

    /// \name For speed of sound at constant Ye
    //@{
    double dPdnb_Ye;
    double dfdnb_Ye;
    //@}

    /** \brief Create a matter object (if \c rel is true, then ensure
	particles are relatvistic
    */
    matter(bool rel);
  
    /// Compute an average inter-ionic spacing
    double average_a();
  
    /// Compute the number-averaged mass number
    double average_A();
  
    /// Compute the number-averaged charge number
    double average_Z();
  
    /// Compute the impurity parameter 
    double impurity();

  private:

    // Make default constructor private
    matter();

    // Make default copy constructor private
    matter(const matter &m);

    // Make default copy constructor via assignment private
    matter &operator=(const matter &m);

  };
  
  /// Output a matter object
  std::ostream &operator<<(std::ostream &os, const matter &m);

#ifndef DOXYGENP
}
#endif

#endif
