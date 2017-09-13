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
#ifndef NM_THERMO_H
#define NM_THERMO_H

#include <o2scl/fermion_mag_zerot.h>
#include <o2scl/fermion_rel.h>

#include "matter.h"
#include "ldrop_crust.h"

namespace crust {

  /** \brief Compute the properties of nuclear matter
   */
  class nm_thermo {

  protected:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Electron thermodynamics
    o2scl::fermion_rel relf;

    /// Convert units
    o2scl::convert_units &cng;

    /// Solver for free_energy()
    o2scl::mroot_hybrids<> cmr2;

    /** \brief The function to solve for nuclear matter
     */
    int solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		  matter &nm);

    /// The mass model
    ldrop_crust *lda;

  public:

    /** \brief Create with the specified mass model, EOS, and 
	unit conversion class
    */
    nm_thermo(ldrop_crust &lc, o2scl::eos_had_temp_base &he, 
	      o2scl::convert_units &conv) :
    cng(o2scl::o2scl_settings.get_convert_units()) {
      het=&he;
      lda=&lc;
    }

    /** \brief Compute the free energy of beta-equilibrated matter 
	at fixed baryon density and temperature

	Also computes the pressure, gibbs energy, and rest mass
	density and internal energy density.
    */
    void calc(matter &nm);

    /** \brief Compute the free energy density of nuclear matter in
	\f$ \mathrm{fm}^{-4} \f$ with the same proton fraction and
	baryon density as a distribution
    */
    void calc_nm_from_dist(matter &m, matter &nm);

    /** \brief Test \ref calc_nm_from_dist() 
     */
    int check_free_energy_x(matter &m, double T, matter &nm);

    /** \brief The hadronic eos

	This needs to be public so it can be changed by \ref
	crust_driver::model() .
     */
    o2scl::eos_had_temp_base *het;

  };

}

#endif
