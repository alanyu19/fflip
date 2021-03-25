.. README for Gitlab
.. Keep text up to date with top-level docs/readme.rst (for sphinx)
.. Those are two separate files, because
.. a) otherwise the links do not work and there is no convincing workaround
.. b) having different representations on gitlab and readthedocs could be helpful


=====
fflip
=====


.. image:: https://gitlab.com/alanyu/fflip/badges/master/build.svg
        :target: https://gitlab.com/alanyu/fflip/pipelines
        :alt: Continuous Integration

.. image:: https://img.shields.io/pypi/v/fflip.svg
        :target: https://pypi.python.org/pypi/fflip#
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

* Documentation: https://fflip.readthedocs.io.


Features
--------

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

1) CHARMM (optional): Required for fflip/uachain. Make sure the command *charmm* can be called from a terminal.
2) OpenMM: Required.
3) openmmtools: Required.
4) mdtraj: Required.
5) Rickflow: Required (https://gitlab.com/Olllom/rickflow/tree/master).

Other comments here.



Tutorials
---------

-  Using fflip with OpenMM/CHARMM:

   1) Parameterize alkane/alkene systems (pending): Tutorial1_.
   2) Optmization for the CHARMM additive lipid FF: Tutorial2_.

.. _Tutorial1: examples/01_alkane_system/placeholder1.ipynb
.. _Tutorial2: examples/02_additive_lipid/


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
