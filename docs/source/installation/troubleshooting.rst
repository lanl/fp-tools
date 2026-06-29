Troubleshooting
===============

This page covers common issues you might encounter when building or using 
FP-Tools, along with suggested solutions.

CMake Cannot Find FFTW
----------------------

**Problem:**
CMake fails to detect an installed FFTW library.

**Solution:**
Ensure FFTW is installed and accessible to CMake. You may need to specify the 
path explicitly using `FFTW_ROOT` or `CMAKE_PREFIX_PATH`:

.. code-block::bash

   cmake -S . -B build -DFFTW_ROOT=/path/to/fftw

Alternatively:

