.. _Anaconda Documentation: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _BOOST: https://www.boost.org/doc/libs/1_53_0/
.. _CMake: https://cmake.org/cmake/help/v3.14/
.. _CMake add_custom_target: https://cmake.org/cmake/help/latest/command/add_custom_target.html
.. _Doxygen: https://www.doxygen.nl/manual/docblocks.html
.. _Eigen: https://eigen.tuxfamily.org/dox/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _Breathe: https://breathe.readthedocs.io/en/latest/
.. _PEP-8: https://www.python.org/dev/peps/pep-0008/
.. _pipreqs: https://github.com/bndr/pipreqs 
.. _LaTeX: https://www.latex-project.org/help/documentation/
.. _W-13 DevOps Manual: https://xcp-confluence.lanl.gov/display/COM/W-13+DevOps
.. _upstream repository: https://re-git.lanl.gov/aea/material-models/cpp_stub
.. _Material Models: https://re-git.lanl.gov/aea/material-models
.. _UNIX group: https://ddw-confluence.lanl.gov/pages/viewpage.action?pageId=150929410

###################
C++ Stub repository
###################

*******************
Project Description
*******************

A stub repository for C++ development projects in W-13.

This repository will contain the necessary setup files to integrate C++ doc
strings, `CMake`_, `Doxygen`_, `Sphinx`_, and `Breathe`_ for a complete build
system with integrated documentation. It will also include the necessary hooks
to commonly used C++ libraries for constitutive modeling. This stub repository
also includes template hooks for integrating C++ code as Abaqus subroutines.

    **NOTE**

    You can use this repo as a stub for fortran projects as well! A step-by-step
    may never happen because c++ is the future of W-13 subroutines.  For now, refer
    to the following references:
   
    * `CMake for Fortran example CMakeLists.txt <https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/languages/fortran/ForFortranExample>`_
    * `CMake documentation starting point <https://cmake.org/cmake/help/v3.14/module/CheckFortranSourceRuns.html>`_
    * `Stack Overflow thread <https://stackoverflow.com/questions/12705562/using-cmake-with-fortran>`_
    * `Doxygen comments for Fortran <https://www.doxygen.nl/manual/docblocks.html#fortranblocks>`_

Information
===========

* Documentation: https://aea.re-pages.lanl.gov/material-models/cpp_stub
* Wiki: https://re-git.lanl.gov/aea/material-models/cpp_stub/-/wikis/home

Developers
==========

* Kyle Brindley: kbrindley@lanl.gov
* Nathan Miller: nathanm@lanl.gov

********************************************
Setting up a new project from this stub repo
********************************************

    **NOTE**

    The repository setup has moved out of the README and into the HTML
    documentation. You can find the Gitlab project setup guide here:
    https://aea.re-pages.lanl.gov/material-models/cpp_stub/gitlab_setup.html

************
Gitlab CI/CD
************

    **NOTE**

    The repository setup has moved out of the README and into the HTML
    documentation. You can find the Gitlab project setup guide here:
    https://aea.re-pages.lanl.gov/material-models/cpp_stub/gitlab_setup.html

************
Dependencies
************

Compilers
=========

* c++11 compiler (listed version number has been tested at some point)

  * g++ >= GNU 4.8.5

Executables
===========

* `CMake`_ >= 3.14
* `Doxygen`_ >= 1.8.5
* `LaTeX`_ >= 2017

Conda Environment
=================

For convenience, the minimal Python environment requirements for the
documentation build are included in ``configuration_files/environment.yaml``.
This file was created from the `pipreqs`_ command line tool and Sphinx
configuration inspection, e.g. the extension packages.

.. code-block:: bash

   $ pwd
   path/to/cpp_stub/
   $ pipreqs --use-local --print --no-pin .

A minimal anaconda environment for building the documentation can be created
from an existing anaconda installation with the following commands.

.. code-block:: bash

   $ conda env create --file configuration_files/environment.yaml

You can learn more about Anaconda Python environment creation and management in
the `Anaconda Documentation`_.

C++ Libraries
=============

    **NOTE**

    Non-admin installations for Eigen and Boost are no longer required.** This
    project is built and deployed against C++ libraries managed in Conda. See the
    Conda environment file and README discussion for non-admin environment
    management.

* `Eigen`_ >= 3.3.7
* `BOOST`_ >= 1.59.0
* error\_tools: https://re-git.lanl.gov/aea/material-models/error_tools
* vector\_tools: https://re-git.lanl.gov/aea/material-models/vector_tools
* abaqus\_tools: https://re-git.lanl.gov/aea/material-models/abaqus_tools
* constitutive\_tools: https://re-git.lanl.gov/aea/material-models/constitutive_tools
* stress\_tools: https://re-git.lanl.gov/aea/material-models/stress_tools
* solver\_tools: https://re-git.lanl.gov/aea/material-models/solver_tools

If not found on the current system or active Conda environment, all of the
``*_tools`` libraries are pulled from their git repos by branch name and built
with their respective cmake files as part of the cmake build for this project.

**************
Build and Test
**************

This project is built with `CMake`_ and uses `Sphinx`_ to build the
documentation with `Doxygen`_ + `Breathe`_ for the c++ API.

Build on sstelmo
================

1) Activate the correct python environment

   .. code-block:: bash

      $ module load python/2020.07-python-3.8
      $ sv3r

2) Create a build directory

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/

      $ mkdir build
      $ cd build

3) Configure ``cmake3``

       This step only needs to be performed once unless you need to specify a
       new CMake configuration for a re-build. Most command line arguments and
       environment variables are stored in the CMake cache. Anything found in cache
       will not be re-configured unless you remove the cache file or clobber the build
       directory.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake3 ..

4) Build various portions of the project

       Most of the project will re-build only as necessary after source updates. Some portions of the documentation
       require a ``make clean`` after documentation source file updates to force a re-build.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Build everything
      $ cmake3 --build .

      # Build only the c++ primary libraries
      $ cmake3 --build src/cpp

5) Locate build files

       The build directory structure may change between version releases. Developers and users are encouraged to become
       familiar with the bash ``find``, ``grep``, and ``tree`` commands to locate build files.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # find c++ libraries and ignore intermediate files with similar extensions
      $ find . \( -name "*.o" -o -name "*.so" -o -name "*.a" \) | grep -vE "\.cpp\."

6) Clean build directory to force a re-build

       **HEALTH WARNING**
      
       The abaqus input files and bash scripts used for integration testing are
       built with the `CMake add_custom_target`_ feature. Consequently, the integration
       test target is *always considered out of date*. The integration test target
       copies all registered input files and the integration test bash script from
       source to build directory. This means the file copy operation is always
       performed when the integration test target is requested in the cmake build
       command, e.g. ``cmake --build .`` or ``cmake --build src/abaqus/tests``. This
       operation is computationally inexpensive with respect to building the
       ``cpp_stub`` source code.
      
       Input files are registered in the ``src/abaqus/tests/CMakeLists.txt`` file
       under the ``ABAQUS_INPUT_FILES`` CMake variable.

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      $ make clean

Test on sstelmo
===============

4) Build tests of the project

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Build c++ tests
      $ cmake3 --build src/cpp/tests

      # Build Abaqus integration tests
      $ cmake3 --build src/abaqus/tests

5) Run the tests

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build

      # Run ctest
      $ ctest

      # Results print to screen
      # View details of most recent test execution including failure messages
      $ less Testing/Temporary/LastTest.log

Convenience build wrappers
==========================

Two build scripts have been created for convenience, ``new_build.sh`` and
``build_docs.sh``. The first will build everything including the library binary,
the test binary, and the documentation. This is the same build script used by
``jenkins_build.sh`` for CI builds and testing. The ``build_docs.sh`` script
only builds the documentation. Both build scripts clobber existing build
directories, reset any bash environment variables, and run the cmake
configuration from scratch.

2) Build everything and run tests

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/

      # Just perform the build (pick one)
      $ ./new_build.sh <cmake build type>
      $ ./new_build.sh None
      $ ./new_build.sh Release

      # Perform tests from PWD
      $ ./build/src/cpp/tests/test_cpp_stub

      # Build and perform tests
      $ ./jenkins_build.sh

3) View test results

   .. code-block:: bash

      # As built directly to PWD
      $ cat results.tex

      # As built by jenkins_build.sh
      $ cat build/src/cpp/tests/*_results.tex
      $ cat *results.tex

4) Display docs

   .. code-block:: bash

      # Sphinx
      $ firefox build/docs/sphinx/html/index.html &

      # Doxygen
      $ firefox build/docs/doxygen/html/index.html &

Building the documentation
==========================

    **HEALTH WARNING**
   
    The sphinx API docs are a work-in-progress. The doxygen API is much more
    useful.

The documentation can be built with ``build_docs.sh``. The steps used in that
shell script are repeated here.

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/
      $ mkdir build/
      $ cd build/

3) Run cmake3 configuration

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ cmake3 ..

4) Build the docs

   .. code-block:: bash

      $ cmake3 --build docs

5) Documentation builds to:

   .. code-block:: bash

      cpp_stub/build/docs/sphinx/html/index.html

6) Display docs

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/sphinx/html/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build/
      $ firefox docs/doxygen/html/index.html &

*******************
Install the library
*******************

Build the entire before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake3 --build .

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/cpp_stub/build
      $ cmake --install . --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install . --prefix /home/$USER/.local

      # Example install to conda environment
      $ conda active my_env
      $ cmake --install . --prefix ${CONDA_DEFAULT_ENV}

      # Example install to W-13 CI/CD conda environment performed by CI/CD institutional account
      $ cmake --install . --prefix /projects/python/release

***********************
Contribution Guidelines
***********************

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

`Sphinx`_ reads in docstrings and other special portions of the code as
reStructured text. Developers should follow
styles in this `Sphinx style guide
<https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#>`_.

Style Guide
===========

This project does not yet have a full style guide. Generally, wherever a style
can't be inferred from surrounding code this project falls back to `PEP-8`_-like
styles. There are two notable exceptions to the notional PEP-8 fall back:

1. `Doxygen`_ style docstrings are required for automated, API from source documentation.
2. This project prefers expansive whitespace surrounding parentheses, braces, and
   brackets.

   * No leading space between a function and the argument list.
   * One space following an open paranthesis ``(``, brace ``{``, or bracket
     ``[``
   * One space leading a close paranthesis ``)``, brace ``}``, or bracket ``]``

An example of the whitespace style:

.. code-block:: bash

   my_function( arg1, { arg2, arg3 }, arg4 );

The following ``sed`` commands may be useful for updating white space, but must
be used with care. The developer is recommended to use a unique git commit
between each command with a corresponding review of the changes and a unit test
run.

* Trailing space for open paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([({[]\)\([^ ]\)/\1 \2/g' <list of files to update>

* Leading space for close paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([^ ]\)\([)}\]]\)/\1 \2/g' <list of files to update>

* White space between adjacent paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([)}\]]\)\([)}\]]\)/\1 \2/g' <list of files to update>
