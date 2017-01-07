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
/** \mainpage Documentation

    \todo
    - Compare to Newton's recent work?
    - Comment on uptick in Douchin et al.'s models?
    - Rerun base and Gs models to get flow figure?
    - Old entry: Make sure makefile targets set temperature to zero 
    when necessary?

    Figures:
    - Free energy of eq. crust with different EOSs (feq_fr)
    - Composition of eq. crust with different EOSs and w/o shell effects
      (feq_fr)
    - Free energy of acc. crust with different EOSs (acc_fr?)
    - Speed of sound in frozen and eq. modes ()
    - Figure separating individual contributions to heat
    - flow figure?

    Equilibrium list:
    - 4 EOSs
    - SLy4 w/o shell 

    Accreted List:
    - baseline SLy4, ashes, with Moller/ame (base)
    - different compositions (heavy, nickel)
    - different EOSs with ashes (Gs, Rs, APR)
    - check delta_n variation (Gs_dn10)
    - with and without experimental/Moller masses (ldrop)
    - with and without nuclear shell effects (noshell)
    - different temperatures (highT)
    - different accretion rates (lowacc, highacc)
    - with simplified rates (simppyc)
    - without highZ fusion (highZ)
    - single nucleus (sna)
    Later:
    - one run with summary info (extra)

    \future 
    - The computation of the pressure using gibbs_energy_density()
    in compute_sna() may be overkill, and it's possible that this
    function could be simplified
    - Move some public functions to protected
    - It would be great if full_eq2() was as accurate and or
    nearly as fast as full_eq()
    - Compute one proton separation energy
    - Directly compute entropy
    - Rewrite full_eq_dist() to properly handle constant baryon density?
    - Always use solver in nucleus_be()?
    - Move the hadronic_eos pointer and temperature into the matter object?

    <b>Basic parts of call stack</b>

    Full equilibrium crust (T=0):
    <ul>
    <li>\ref crust::crust_driver::full_eq()</li>
    <ul>
    <li>\ref crust::crust_driver::compute_sna()</li>
    <ul>
    <li>\ref crust::crust_driver::delta_ZN()</li>
    <li>Minimize \ref crust::sna_thermo::free_energy_sna_neut::operator()()</li>
    <ul>
    <li>\ref crust::sna_thermo::free_energy_sna_fix_nb_nn()</li>
    <li>\ref crust::sna_thermo::baryon_density_sna()</li>
    <li>\ref crust::sna_thermo::free_energy_sna()</li>
    </ul>
    </ul>
    </ul>
    </ul>
    
    Accreted crust:
    <ul>
    <li>\ref crust::crust_driver::acc()</li>
    <ul>
    <li>\ref crust::rxns::emit_neutron()</li>
    <ul>
    <li>\ref crust::dist_thermo::gibbs_energy_per_baryon_cell()</li>
    </ul>
    <li>\ref crust::rxns::elec_capture()</li>
    <li>\ref crust::rxns::beta_decay()</li>
    <li>\ref crust::rxns::pyc_fusion()</li>
    <li>\ref crust::dist_thermo::free_energy_dist()</li>
    <li>\ref crust::dist_thermo::gibbs_energy_dist()</li>
    </ul>
    </ul>

    <b>Generic Documentation</b>
    
    From the TOV equations
    \f[
    \frac{dP}{dr} = - \frac{G \varepsilon m}{r^2} 
    \left( 1-\frac{2 G m}{r}\right)^{-1}
    \f]
    and also 
    \f[
    \frac{d n_B}{dr} = - 4 \pi r^2 n_B \left( 1-\frac{2 G m}{r}\right)^{-1/2}
    \f]
    then we have 
    \f[
    d n_B = \frac{4 \pi r^4 n_B}{G \varepsilon m} 
    \left( 1- \frac{2 G m}{r}\right)^{1/2} d P
    \f]
    and discretizing the grid, then for a shell with index \f$ i \f$ ,
    \f[
    d P = \frac{P_{i+1}-P_{i-1}}{2}
    \f]
    and \f$ d n_B \f$ is the total number of baryons in shell \f$ i \f$ .
    To compute the gravitational mass in the shell,
    we use
    \f[
    d m_G = - 4 \pi r^2 \varepsilon
    \f]
    and then 
    \f[
    d m_G = 8 \pi r^3 \left( \frac{r}{2 G m} -1\right) d P
    \f]
    If we keep \f$ m \f$ and \f$ r \f$ constant over all of the shells,
    then there are constants \f$ A_n \f$ and \f$ A_m \f$ such that 
    \f[
    d n_B = \frac{A_n n_B}{\varepsilon} d P
    \f]
    \f[
    d m_G = A_m d P
    \f]

    <b>References</b>

    \anchor Beard10 Beard10:
    M. Beard, A.V. Afansjev, L.C. Chamon, L.R. Gasques, M. Wiescher, 
    D.G. Yakovlev, At. Data Nucl. Data Tables 96 (2010) 541.

    \anchor Yakovlev06 Yakovlev06:
    D.G. Yakovlev, L.R. Gasques, A.V. Afanasjev, M. Beard, and M. Wiescher,
    Phys. Rev. C 74 (2006) 035803.
    
    \anchor Yakovlev06b Yakovlev06b:
    D.G. Yakovlev, L. Gasques, and M. Wiescher
    Mon. Not. R. Astron. Soc. 371 (2006) 1322.

    \anchor Yakovlev10 Yakovlev10:
    D.G. Yakovlev, M. Beard, L.R. Gasques, M. Wiescher, Phys. Rev. C
    82 (2010) 044609

    

*/

