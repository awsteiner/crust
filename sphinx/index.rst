Crust
=====

A code for computing the structure of a neutron star 
crust. 

This documentation was generated from git commit

.. include:: commit.rst

Based partially on work in [Steiner08]_, this
code was written for use in [Steiner12]_ and
then used again in [Deibel14]_ (which was
itself based on [Steiner09]_ ).

Basic call stack
----------------

Full equilibrium crust (T=0):

- :cpp:func:`crust::crust_driver::full_eq()`
  
  * :cpp:func:`crust::crust_driver::compute_sna()`
    
- :cpp:func:`crust::crust_driver::delta_ZN()`
  
- Minimize
  :cpp:func:`crust::sna_thermo::free_energy_sna_neut::operator()`
       
  * crust::sna_thermo::free_energy_sna_fix_nb_nn()
  * crust::sna_thermo::baryon_density_sna()
  * crust::sna_thermo::free_energy_sna()

Accreted crust:

- :cpp:func:`crust::crust_driver::acc()`
  
  * :cpp:func:`crust::rxns::emit_neutron()`
    
  * :cpp:func:`crust::dist_thermo::gibbs_energy_per_baryon_cell()`
    
  * :cpp:func:`crust::rxns::elec_capture()`
    
  * :cpp:func:`crust::rxns::beta_decay()`
    
  * :cpp:func:`crust::rxns::pyc_fusion()`
    
  * :cpp:func:`crust::dist_thermo::free_energy_dist()`
    
  * :cpp:func:`crust::dist_thermo::gibbs_energy_dist()`

.. toctree::
   :maxdepth: 2

   classdoc
   bib
   
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

   
  
