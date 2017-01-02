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
#include "matter.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

matter::matter(bool rel) {
  
  n=new fermion(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  p=new fermion(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_proton),2.0);
  
  n->inc_rest_mass=true;
  p->inc_rest_mass=true;
  n->non_interacting=false;
  p->non_interacting=false;
  
  e.init(o2scl_settings.get_convert_units().convert
	 ("kg","1/fm",o2scl_mks::mass_electron),2.0);
  e.inc_rest_mass=true;
  e.non_interacting=true;
}

double matter::average_a() {
  double ntot=0.0;
  for(size_t j=0;j<dist.size();j++) {
    ntot+=dist[j].n;
  }
  return cbrt(3.0/4.0/o2scl_const::pi/ntot);
}
  
double matter::average_A() {
    
  double ntot=0.0;
  double Antot=0.0;
    
  for(size_t j=0;j<dist.size();j++) {
    
    ntot+=dist[j].n;
    Antot+=(dist[j].N+dist[j].Z)*dist[j].n;
      
  }
    
  return Antot/ntot;
}
  
double matter::average_Z() {
    
  double ntot=0.0;
  double Zntot=0.0;

  for(size_t j=0;j<dist.size();j++) {
      
    ntot+=dist[j].n;
    Zntot+=dist[j].Z*dist[j].n;
      
  }
    
  return Zntot/ntot;
}
  
double matter::impurity() {
    
  if (dist.size()==1) return 0.0;

  // First compute average proton number
  double ave_Z=average_Z();

  double sum=0.0, ntot=0.0;
    
  for(size_t i=0;i<dist.size();i++) {
    sum+=dist[i].n*pow(dist[i].Z-ave_Z,2.0);
    ntot+=dist[i].n;
  }

  return sum/ntot;
}

std::ostream &o2scl::operator<<(std::ostream &os, const matter &m) {
  os << "----------------------------------------------------------" 
     << std::endl;
  os << "nn,np,ne: " << m.n->n << " " << m.p->n << " "
     << m.e.n << std::endl;
  os << "Dist size: " << m.dist.size() << std::endl;
  os << "index,Z,N,n,mu,zeta:" << std::endl;
  for(size_t i=0;i<m.dist.size();i++) {
    os << i << " " << m.dist[i].Z << " " << m.dist[i].N << " "
       << m.dist[i].n << " " << m.mu[i] << " " 
       << m.zeta[i] << std::endl;
  }
  os << "ed,pr,nb,rho,fr,gb: ";
  os << m.ed << " " << m.pr << std::endl << m.nb << " " 
     << m.rho << " " << m.fr << " " << m.gb << std::endl;
  os << "mun,mup: " << m.mun << " " << m.mup << std::endl;
  os << "----------------------------------------------------------";
  return os;
}
      
