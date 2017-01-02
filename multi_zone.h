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

#include "dist_thermo.h"
#include "matter.h"

#include <o2scl/nucmass.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/table_units.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief Multi-zone stability test
   */
  class multi_zone {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Schwarzchild radius
    double schwarz_km;

    /// Object to convert units
    convert_units cng;

    /// To interpret element labels in accreted crust table
    nucmass_info nmi;

    /// Number of nuclei
    ubvector nN;

    /// Number of neutrons
    ubvector nn;

    /// Desc
    int sdebug;

    /// Desc
    ubvector x_solve;

    /// Desc
    size_t n_min;

    /// Desc
    ubvector min_scale;

    /// Desc
    ubvector solve_scale;

    /// Desc
    double nn_base;

    /// Conversion between 
    double conv;

    /// Target number of nuclei
    ubvector nN_base;

    /// Gravitational mass grid
    ubvector gm;

    /// Gravitational mass at outer boundary
    double gmtot;

    /// Total number of neutrons over all zones
    double nntot;

    /// Pressure grid (in \f$ M_{\odot}/km^3 \f$)
    ubvector P;

    /// Radial grid (in km)
    ubvector r;

    /// Energy density grid (in \f$ M_{\odot}/km^3 \f$)
    ubvector eps;

    /// Number of variables to solve for
    size_t n_solve;

    /// Number of zones
    size_t n_zones;

    /// The index of the first zone in the table
    size_t first_zone;
  
    /// If true, set the target number of nuclei in solve_fun()
    bool first_time;

    /// Storage for the zones
    std::vector<matter *> zones;

    /// To compute thermodynamics of each zone
    dist_thermo *dtp;

    /// Temperature
    double Tptr;

    /// Desc
    double Pbottom;

    /// Desc
    double Ptop;

    /// Solver
    o2scl::mroot_hybrids<> gmr;

    /// Minimizer
    o2scl::mmin_simp2<> gms;

    multi_zone();

    /** \brief Read table into a vector of \ref matter objects
     */
    void read_zones(size_t nz, o2scl::table_units<> &t, size_t fz);

    /** \brief Desc
    */
    int solve_fun(size_t nv, const ubvector &x, 
		  ubvector &y);
  
    /** \brief Desc
    */
    double min_fun(size_t nv, const ubvector &x);

    /** \brief 
     */
    void test_stability(double T, dist_thermo &dt, size_t nz, 
			o2scl::table_units<> &t, size_t fz);
			

  };

#ifndef DOXYGENP
}
#endif
