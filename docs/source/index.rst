FP-Tools
========

`FP-Tools <https://github.com/lanl/fp-tools/>`__ is a molecular dynamics 
trajectory analysis code written in C++ and developed at Los Alamos National 
Laboratory. It was originally written with a focus on **first-passage 
diagnostics**, but has since grown to include a larger set of standard 
molecular dynamics analysis routines. At present, the code is a self-contained
toolkit for calculating structural and dynamical properties or materials from 
molecular dynamics trajectories.

Features
--------

- **Efficient trajectory parsing** for LAMMPS, VASP (XDATCAR), and XYZ formats.
- **Standard analyses**:

    - **Mean squared displacement (MSD)**
    - **Radial distribution function (RDF)**
    - A **correlation function calculator** with built-in support for:
    
        - Velocity autocorrelation function (VACF)
        - Longitudinal and transverse current correlation functions
        - Dynamic structure factor, :math:`S(q,\omega)`
        - Static structure factor, :math:`S(q)`
        - Intermediate scattering function, :math:`F(q,t)`

    - A **custom correlation function calculator** for arbitrary user-specified 
      time series

.. note::
   For details on the algorithms in FP-Tools see the master document ___

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   section#installation,
   section#usage,
   section#theory,
   section#data-analysis,
   section#community {
       display:none;
   }
   </style>

.. toctree::
   :hidden:

   acknowledge_fptools
   contribute
   contact

Getting Started
---------------
.. toctree::
   :caption: GETTING STARTED
   :maxdepth: 1
   :hidden:

   installation/building.rst
   installation/quick-start.rst
   installation/troubleshooting.rst
   installation/examples.rst
   installation/release.rst

Usage
-----
.. toctree::
   :caption: USAGE
   :maxdepth: 1
   :hidden:

   usage/parameters.rst
   usage/get_started.rst

Data Analysis
-------------
.. toctree::
   :caption: DATA ANALYSIS
   :maxdepth: 1
   :hidden:

   analysis/visualization.rst


Theoretical Background
----------------------
.. toctree


