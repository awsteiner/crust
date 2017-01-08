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
#ifndef CRUST_FIT_H
#define CRUST_FIT_H

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

namespace crust {

/** \brief Fit a \ref crust::ldrop_crust object to nuclear mass data

    If \ref fit_moller is true, then the model is fit to
    the data from \ref moller, otherwise it is fit to the
    data from \ref ame.

 */
class crust_fit : public o2scl::nucmass_fit {

 public:

  crust_fit();

  /// Pointer to mass formula
  crust::ldrop_crust *lda;
  
  /// If true, fit FRDM instead of AME (default true)
  bool fit_moller;

  /// The experimental mass model
  o2scl::nucmass_ame_exp ame;
  
  /// Moller et al. mass model
  o2scl::nucmass_mnmsk moller;

  /** \brief Read fit from a file
   */
  int read_fit(std::vector<std::string> &sv, bool itive_com);

  /** \brief Fit the nuclear mass formula
   */
  int perform_fit(std::vector<std::string> &sv, bool itive_com);

};

}
 
#endif
