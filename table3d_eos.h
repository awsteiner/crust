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
#ifndef O2SCL_TABLE3D_EOS_H
#define O2SCL_TABLE3D_EOS_H

#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/fermion.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/table3d.h>

namespace crust {

  /** \brief EOS from a table3d for the NS crust

      The table, \ref tab, is created in \ref
      crust_driver::make_table().
  */
  class table3d_eos : public o2scl::eos_had_temp_base {
    
  public:

    /// The table which stores the EOS
    o2scl::table3d tab;

    table3d_eos() {
    }

    virtual ~table3d_eos() {
    }

    /** \brief Equation of state as a function of density
    */
    virtual int calc_e(o2scl::fermion &ne, o2scl::fermion &pr,
		       o2scl::thermo &th) {
      size_t ix, iy, ix2, iy2;

      // Manually perform two-dimensional linear interpolation

      tab.lookup_x(ne.n,ix);
      tab.lookup_y(pr.n,iy);
      ix2=ix++;
      iy2=iy++;
      if (ix2+1>tab.get_nx()) {
	ix2--;
	ix--;
      }
      if (iy2+1>tab.get_ny()) {
	iy2--;
	iy--;
      }

      double nn0=tab.get_grid_x(ix);
      double nn1=tab.get_grid_x(ix2);
      double np0=tab.get_grid_y(iy);
      double np1=tab.get_grid_y(iy2);

      double ed00=tab.get(ix,iy,"ed");
      double ed01=tab.get(ix,iy2,"ed");
      double ed10=tab.get(ix2,iy,"ed");
      double ed11=tab.get(ix2,iy2,"ed");

      double pr00=tab.get(ix,iy,"pr");
      double pr01=tab.get(ix,iy2,"pr");
      double pr10=tab.get(ix2,iy,"pr");
      double pr11=tab.get(ix2,iy2,"pr");

      double mun00=tab.get(ix,iy,"mun");
      double mun01=tab.get(ix,iy2,"mun");
      double mun10=tab.get(ix2,iy,"mun");
      double mun11=tab.get(ix2,iy2,"mun");

      double mup00=tab.get(ix,iy,"mup");
      double mup01=tab.get(ix,iy2,"mup");
      double mup10=tab.get(ix2,iy,"mup");
      double mup11=tab.get(ix2,iy2,"mup");

      double ed0=ed00+(ed01-ed00)*(pr.n-np0)/(np1-np0);
      double ed1=ed10+(ed11-ed10)*(pr.n-np0)/(np1-np0);
      th.ed=ed0+(ed1-ed0)*(ne.n-nn0)/(nn1-nn0);

      double pr0=pr00+(pr01-pr00)*(pr.n-np0)/(np1-np0);
      double pr1=pr10+(pr11-pr10)*(pr.n-np0)/(np1-np0);
      th.pr=pr0+(pr1-pr0)*(ne.n-nn0)/(nn1-nn0);

      double mun0=mun00+(mun01-mun00)*(pr.n-np0)/(np1-np0);
      double mun1=mun10+(mun11-mun10)*(pr.n-np0)/(np1-np0);
      ne.mu=mun0+(mun1-mun0)*(ne.n-nn0)/(nn1-nn0);

      double mup0=mup00+(mup01-mup00)*(pr.n-np0)/(np1-np0);
      double mup1=mup10+(mup11-mup10)*(pr.n-np0)/(np1-np0);
      pr.mu=mup0+(mup1-mup0)*(ne.n-nn0)/(nn1-nn0);

      th.ed=tab.interp(ne.n,pr.n,"ed");
      th.pr=tab.interp(ne.n,pr.n,"pr");
      ne.mu=tab.interp(ne.n,pr.n,"mun");
      pr.mu=tab.interp(ne.n,pr.n,"mup");

      th.en=0.0;
      ne.nu=ne.mu;
      pr.nu=pr.mu;
      ne.ms=ne.m;
      pr.ms=pr.m;

      return 0;
    }

    /// A wrapper for calc_e() since this is only zero temperature
    virtual int calc_temp_e(o2scl::fermion &ne, o2scl::fermion &pr, double T,
			    o2scl::thermo &th) {
      return calc_e(ne,pr,th);
    }

    virtual int calc_p(o2scl::fermion &ne, o2scl::fermion &pr,
		       o2scl::thermo &th) {
      return 0;
    }
    
    virtual int calc_temp_p(o2scl::fermion &ne, o2scl::fermion &pr,
			    double T, o2scl::thermo &th) {
      return 0;
    }
    

  };

}

#endif
