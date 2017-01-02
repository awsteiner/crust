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
#ifndef CRUST_H
#define CRUST_H

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
#include "eigen.h"
#include "ldrop_crust.h"
#include "nm_thermo.h"
#include "sna_thermo.h"
#include "dist_thermo.h"
#include "pyc_rates.h"
#include "crust_fit.h"
#include "rxns.h"
#include "table3d_eos.h"
#include "multi_zone.h"

#ifndef DOXYGENP
namespace o2scl {
#endif

  /// To sort the nuclear distribution with density (currently unused)
  bool compare_density(const nucleus &n1, const nucleus &n2);

  /** \brief To sort the nuclear distribution with proton number and then 
      neutron number
  */
  bool compare_Z(const nucleus &n1, const nucleus &n2);

  /** \brief Main crust driver class
   */
  class crust_driver {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
  protected:

    /// \name Initial values for equilibrium crust
    //@{
    /// Proton number
    int feq_Z;
    /// Neutron number
    int feq_N;
    /// Neutron drip density
    double feq_nn;
    //@}

    /// \name Baryon density grid for equilibrium crust
    //@{
    /// Initial
    double feq_nb;
    /// Final
    double feq_end;
    /// Stepsize
    double feq_dnb;
    //@}
    
    /// Directory for input files
    std::string input_dir;
    
    /// True if the skyrme model has been set
    bool skyrme_set;

    /// The object to perform multizone calculations
    multi_zone mz;

    /// Full heating without neutrino losses
    double heat_full;

    /** \brief Density at which to compute the summary objects 
	(default \f$ 10^{12} \f$)
    */
    double rho_summary;

    /// Output filename for flow
    std::string flow_fn;

    /// Ensure all reactions happen e.g. capturing after fusing (default true)
    bool more_reactions;

    /// Object which fits the associated mass formula to experiment
    crust_fit cf;

    /// If true, include pasta deformation (default false - does not yet work)
    bool use_pasta;

    /// Test for reactions
    rxns rn;

    /** \brief If true, include the translational free energy for nuclei
	(default false) 
    */
    bool inc_nuc_trans;

    /** \brief Initial distribution type (default "nickel")
      
	Either "nickel", "heavy", or "ashes"
    */
    std::string dist_type;

    /// If true, compute the heating rate of matter (default true)
    bool calc_heating;

    /// \name The neutron and proton objects used by the ldrop_mass object
    //@{
    fermion n_lda, p_lda;
    fermion n_rel_lda, p_rel_lda;
    //@}

    /// The hadronic eos (default points to \c sk)
    o2scl::eos_had_temp_base *het;
  
    /// The homogeneous EOS (the default is SLy4)
    o2scl::eos_had_skyrme sk;
  
    /// APR EOS
    o2scl::eos_had_apr apr;

    /// RMF EOS
    o2scl::eos_had_rmf rmf;

    /// \name Tabulated EOS
    //@{
    table3d_eos t3d;
    //@}

    /** \brief Temperature (in \f$ \mathrm{fm}^{-1} \f$) (the default
	is \f$ 10^{8} \f$ K 
    **/
    double Tptr;

    /// To convert elements to strings
    nucmass_info nmi;

    /** \brief Minimizer for full equilibrium crust (used in \ref
	compute_sna() and \ref compute_sna_dist()).
    */
    min_cern<> cm;
    
    /// Verbosity (default 1)
    int verbose;

    /// Allow emission of protons and alpha particles (default false)
    bool allow_palpha;

    /// (default false)
    bool feq_fix_mode;

    /// Proton number for \ref feq_fix_mode
    int fix_Z;

    /// Neutron number for \ref feq_fix_mode
    int fix_N;
  
    /// \name Check parameter (default \ref check_none)
    //@{
    int check;
    static const int check_none=0;
    static const int check_mass_density=1;
    static const int check_free_energy_sna=2;
    static const int check_free_energy_dist=3;
    static const int check_free_energy_nm_x=4;
    static const int check_ldrop_derivs=5;
    static const int check_mixture=6;
    static const int check_sfactor=7;
    static const int check_free_energy_cell=8;
    static const int check_pressure=9;
    static const int check_rate2=10;
    //@}

    /** \brief Factor to increase pressure or density by for 
	accreted crust (default 1.01)
    */
    double acc_inc_factor;

    /// Set the directory for the input files
    int set_input_dir(std::vector<std::string> &sv, bool itive_com);

    /// Initialize the distribution (called by \ref acc() )
    int init_dist(std::string mode, matter &m);
  
    /** \brief Compute the single nucleus approximation at 
	the baryon density \c nb
	
	Called by full_eq(). Uses the given values of \c Z and \c N as
	initial guesses, updates these values and also returns \c rho
	and \c fr. Uses \ref sna_thermo::free_energy_sna_neut and
	\ref dist_thermo::mass_density().

	Sets temperature to zero automatically.
    */
    int compute_sna(double nb, double T, matter &m, bool debug=false);

    /** \brief Compute the single nucleus approximation at 
	the baryon density \c nb

	Called by full_eq2(). Uses the given values of \c Z and \c N
	as initial guesses, updates these values and also returns \c
	rho and \c fr. Uses \ref dist_thermo::free_energy_dist_sna()
	and \ref dist_thermo::mass_density().
    */
    int compute_sna_dist(double nb, matter &m, double T, bool debug=false);

    /// Select a model for the EOS
    int model(std::vector<std::string> &sv, bool itive_com);

    /// Make an EOS table for interpolation
    int make_table(std::vector<std::string> &sv, bool itive_com);

    /** \brief Compute crust in full_equilibrium
    */
    int full_eq(std::vector<std::string> &sv, bool itive_com);
  
    /** \brief Compute crust in full_equilibrium (V2)

	This version uses free_energy_dist() instead of the "sna" functions
    */
    int full_eq2(std::vector<std::string> &sv, bool itive_com);

    /// Compute crust in full_equilibrium with a distribution
    int full_eq_dist(std::vector<std::string> &sv, bool itive_com);

    /// Compute accreted crust
    int acc(std::vector<std::string> &sv, bool itive_com);

    /// Test neutron drip
    int test_ndrip(std::vector<std::string> &sv, bool itive_com);

    /** \brief Set temperature (takes argument in Kelvin and sets value in 
	\f$ \mathrm{fm}^{-1} \f$.
     */
    int tptr(std::vector<std::string> &sv, bool itive_com);
    
  public:

    /** \brief Remove nuclei in the distribution which have vanishing
	contribution

	Used by \ref dist_switch_gb() and in rxns for summary functions.
     */
    int prune_distribution(std::vector<nucleus> &dist);

    /// Compute shear modulus in a neutron star crust
    tov_shear tsh;

    /// Compute the thermodyanmics of a matter distribution
    dist_thermo dt;

    /// Properties of nuclear matter
    nm_thermo nmt;

    /// Properties of crust matter with a single nucleus
    sna_thermo snat;

    /// The mass model
    ldrop_crust lda;

    /// Perform a check 
    int check_fun(std::vector<std::string> &sv, bool itive_com);

    /// \name Variables needed by the reactions object \ref rn
    //@{
    /// Electron capture heating fraction (default 0.25, typically 1/6 to 1/4)
    double ec_heating;

    /// Convert units
    convert_units cng;

    /// If true, use simplified pycnonuclear rates (default true)
    bool simple_pyc;

    /// Compute pycnonuclear rates
    pyc_rates pyc;
    //@}

    /// Compute loop sizes for compute_sna()
    int delta_ZN(int &Z, int &N, int &deltaZ_lo, int &deltaZ_hi, 
		 int &deltaN_lo, int &deltaN_hi);
    
    /** \brief The limit for pycnonuclear fusion (old)

	This function returns the upper limit in proton number
	for fusion given a fixed mass density
    */
    bool pycno_allowed(double Z1, double Z2, double rho);

    /** \brief Update reaction flow with a reaction between 
	\f$ (Z_1,N_1) \rightarrow (Z_2,N_2) \f$
    */
    int update_flow(int Z1, int N1, int Z2, int N2, std::string type,
		    double dheat, double rho);
    
    /** \brief Check for consistency of the functions free_energy_cell()
	and free_energy_dist()
    */
    int check_free_energy_cell_fun();
  
    /** \brief Switch distributions and compute the heating, 
	assuming the pressure was fixed
    */
    bool dist_switch_gb(matter &m, matter &m_new, 
			int &cnt, bool &restart, double &heat, 
			double T, std::string type, int debug=0);

    /// Constructor
    crust_driver();

    ~crust_driver();

    /// Basic wrapper handling command-line (called from main())
    int run(int argc, char *argv[]);

  };

#ifndef DOXYGENP
}
#endif

#endif
