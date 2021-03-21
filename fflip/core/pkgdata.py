# -*- coding: utf-8 -*-

"""Easy access to data files that are embedded in the caller's python package.
The argument relative_path is always a path relative to the caller's module.
For examples, take a look at the test_relpath module."""

from __future__ import absolute_import, division, print_function

import inspect
import os
import pkg_resources


# ========================================================================================
# Important notice: The functions in this module (relpath) should by no
# means call other functions from this module. As wrapper functions to
# pkg_resources function, they explicitly look into their call stack to
# extract the caller's module.
# ========================================================================================

# ========================================================================================
# The functions in this module only work for files in the tree of the coffe directory.
# The pkg_resource module only works on packaged files,
# as long as include_package_data=True is set in setup.py
# ========================================================================================

# ========================================================================================
# Further reading:
# - http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
# - http://setuptools.readthedocs.io/en/latest/pkg_resources.html#resourcemanager-api
# - http://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
# - https://stackoverflow.com/questions/1395593/managing-resources-in-a-python-project
# ========================================================================================

def abspath(relative_path): # type (object) -> object
    """Get file from a path that is relative to caller's module.
    Returns:    absolute path as string"""
    caller = inspect.stack()[1]
    mod = inspect.getmodule(caller[0])
    return os.path.normpath(pkg_resources.resource_filename(mod.__name__, relative_path))
