.. README for Gitlab
.. Keep text up to date with top-level docs/readme.rst (for sphinx)
.. Those are two separate files, because
.. a) otherwise the links do not work and there is no convincing workaround
.. b) having different representations on gitlab and readthedocs could be helpful


=====
fflip
=====


.. image:: https://gitlab.com/alanyu/fflop/badges/master/build.svg
        :target: https://gitlab.com/alanyu/fflop/pipelines
        :alt: Continuous Integration

.. image:: https://img.shields.io/pypi/v/fflip.svg
        :target: https://pypi.python.org/pypi/fflop#
        :alt: pypi

.. image https://img.shields.io/travis/Olllom/coffe.svg
        :target: https://travis-ci.org/Olllom/coffe

.. image:: https://readthedocs.org/projects/coffe/badge/?version=latest
        :target: https://coffe.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. see https://anaconda.org/conda-forge/plotly/badges for conda badges

Force Field of Lipids Parametrization

A Python package for automation of molecular simulations and optimization of lipid force field parameters.

* Free software: GNU General Public License v3

* Documentation: https://fflop.readthedocs.io.


Features
--------

* Automated composure of simulation systems in Gromacs
* Automated simulation workflows in Gromacs
* Interaction with slurm and torque submission systems
* Advanced stopping conditions for Gromacs simulations (TBD)
* ... more to come


Getting started
---------------


Instructions for downloading and installing the code can be found here: Installation_.

.. _Installation: docs/notebooks/01_getting_started.ipynb



Third-party Software
--------------------

The following third-party software is **not** required to install coffe.
However, a lot of coffe's functionality depends on molecular simulation packages.
In order to use program-specific functions, these programs have to be installed.

1) CHARMM (optional): Required for everything in fflip/tail
   Make sure the command *charmm* can be called from a terminal.
2) OpenMM (optional): Required for everything in fflip/tail.
3) openmmtools (optional): Required for everything in fflip/tail.
#It is absolutely not necessary to have all of these programs installed to benefit from fflip.
#(In most cases you will use only one simulation program anyway.) fflip will (currently not) notify you, if you hit
#a dependency that is not installed on your system.




Tutorials
---------

-  Using fflip with CHARMM (and OpenMM):

   1) Building alkane/alkene systems: Tutorial1_.
   2) Running CHARMM simulations: Tutorial2_.
   3) Running OpenMM simulations: Tutorial3_.
   4) Extract information from simulations: Tutorial4_.

.. _Tutorial1: examples/01_building_alkane_system/placeholder1.ipynb
.. _Tutorial2: examples/02_01_charmm_sim/placeholder2.ipynb
.. _Tutorial3: examples/02_02_openmm_sim/placeholder3.ipynb
.. _Tutorial4: examples/03_getting_properties/placeholder4.ipynb


Developers' Guide
-----------------

1) `Coding conventions`_
2) `Notes on how to add modules, classes, functions, etc.`_
3) `Contributing to the code`_
4)  `Notes to core developers`_

.. _Coding conventions: docs/notebooks/02_coding_conventions.ipynb
.. _Notes on how to add modules, classes, functions, etc.: docs/notebooks/03_adding_stuff.ipynb
.. _Contributing to the code: CONTRIBUTING.rst
.. _Notes to core developers: docs/notebooks/04_mergerequests.ipynb


Related Work
------------

Currently none.


Credits
---------

Find out about the authors here_.

.. _here: AUTHORS.rst
