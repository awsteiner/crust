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
#ifndef RXNS_H
#define RXNS_H

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

#include "matter.h"

namespace crust {

  // Forward definition
  class crust_driver;
  
  /// Apply nuclear reactions to a matter distribution
  class rxns {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    /** \brief Compute electron captures by considering only one cell
	(default false)
     */
    bool ec_one_cell;
    
    /// True if carbon has been fused
    bool fused_6;

    /// True if oxygen has been fused
    bool fused_8;

    /// True if neon has been fused
    bool fused_10;

    /// True if magnesium has been fused
    bool fused_12;
  
    /** \brief The fractional shift in the total number of nuclei
	(default 0.01)
    */
    double delta_n;

    /** \brief Mass accretion rate in solar masses per year 
	(default \f$ 10^{-10} \f$)
    */
    double mdot;
    
    rxns();
    
    /** \brief Electron capture

        Tries to electron capture for every nucleus in the 
	distribution. 

	\comment

	1/9/17: I think it's always fixed pressure now
	
	By default, the new configuration is computed at fixed
	particle number. If \ref crust_driver::gibbs_press is true,
	then the new configuration is computed at fixed pressure.

	\endcomment
    */
    int elec_capture(crust_driver *a, matter &m, matter &m_new, double T, 
		     int &cnt, double &heat);

    /// Summarize the electron capture thresholds
    int ec_summary(crust_driver *a, matter &m, matter &m_new, 
		   double T, ubmatrix &gpb_store,
		   ubmatrix &m_store);
    
    /** \brief Add a nucleus to the distribution if it's missing, and 
	return its index in \c ix
    */
    int add_missing(std::vector<o2scl::nucleus> &trial, int newZ, int newN,
		    size_t &ix);

    /// Beta decay
    int beta_decay(crust_driver *a, matter &m, matter &m_new, double T, 
		   int &cnt, double &heat);
		   

    /// Neutron emission
    int emit_neutron(crust_driver *a, matter &m, matter &m_new, double T, 
		     int &cnt, double &heat);

    /// Summarize the neutron emission thresholds
    int en_summary(crust_driver *a, matter &m, matter &m_new, 
		   double T, ubmatrix &gpb_store,
		   ubmatrix &m_store);

    /// Neutron capture
    int neut_capture(crust_driver *a, matter &m, matter &m_new, double T, 
		     int &cnt, double &heat);

    /// Generalized reaction
    int gen_reaction(crust_driver *a, matter &m, matter &m_new, double T, 
		     int &cnt);

    /// Pycnonuclear fusion
    int pyc_fusion(crust_driver *a, matter &m, matter &m_new, double T, 
		   int &cnt, double &heat);

    /// Emit proton or alpha particle
    int emit_fragment(crust_driver *a, matter &m, matter &m_new, double T, 
		      int &cnt, double &heat);
  };
  
}

#endif
