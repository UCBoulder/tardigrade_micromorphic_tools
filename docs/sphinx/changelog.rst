.. _changelog:


#########
Changelog
#########

******************
1.4.0 (unreleased)
******************

Internal Changes
================
- Changed the fuzzyEquals in the tests to use BOOST_TEST (:pull:`7`). By `Nathan Miller`_.
- Completed first pass of improving the computational efficiency (:pull:`9`). By `Nathan Miller`_.

Breaking changes
================
- Added vector-based Jacobians (:pull:`1`, :pull:`4`, :pull:`5`). By `Nathan Miller`_.

******************
1.3.2 (2023-09-29)
******************

Internal Changes
================
- Add draft GitHub release action. By `Kyle Brindley`_.

******************
1.3.1 (2023-07-24)
******************

Breaking changes
================
- Change project, package, and namespace from 'micromporphic tools' to 'tardigrade micromorphic tools' (:issue:`14`,
  :merge:`18`). By `Kyle Brindley`_.

******************
1.2.1 (2023-07-11)
******************

Internal Changes
================
- Updated build to build the `dev` documentation before `main` (:merge:`5`). By `Nathan Miller`_.
- Updated build to enable documentation deployment (:merge:`4`). By `Nathan Miller`_.
- Updated build to latest cpp_stub (:merge:`3`). By `Nathan Miller`_.
- Updated tests to BOOST test (:merge:`2`). By `Nathan Miller`_.
- Removed the older build scripts in favor of documented, direct cmake commands (:issue:`4`, :merge:`12`). By `Kyle
  Brindley`_.
- Run CI tests with a project specific CI Conda environment (:issue:`3`, :merge:`13`). By `Kyle Brindley`_.
- Use setuptools_scm for Git tag version number (:issue:`5`, :merge:`14`). By `Kyle Brindley`_.
- Add the tardigrade license and meta data (:issue:`10`, :merge:`16`). By `Kyle Brindley`_.


******************
0.1.6 (2021-08-26)
******************

Internal Changes
================
-upgrade build system to CMake3

******************
0.1.5 (2021-07-19)
******************

Documentation
=============
- Update project setup instructions from Atlassian to Gitlab workflows (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.

Internal Changes
================
- Convert README from markdown to restructured text (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.
- Separate Abaqus integration test setup from Abaqus integration ctest declaration. Enables documentation build
  dependencies on Abaqus integration test input files without requiring Abaqus test execution on systems with no Abaqus
  installation (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.


******************
0.1.4 (2021-07-13)
******************

Internal Changes
================
- Upstream project settings update to set default merge-request branch. By `Kyle Brindley`_.

******************
0.1.3 (2021-07-13)
******************

- Migrate from ddw-bibucket.lanl.gov to re-git.lanl.gov and convert to Gitlab CI/CD (:issue:`1`, :merge:`1`). By `Kyle
  Brindley`_.

******************
0.1.2 (2021-07-01)
******************

Internal Changes
================
- Use Git SCM tags for semantic versioning (:jira:`702`, :pull:`50`). By `Kyle Brindley`_.
- Master branch production release logic for CD, including automated micro-version bumps (:jira:`702`, :pull:`50`). By `Kyle
  Brindley`_.


******************
0.1.1 (2021-06-15)
******************

Bug Fixes
=========
- Corrected bug in `cpp_stub.cpp` in the map of `ddsdde` to `DDSDDE` due to using `spatialDimensions` instead
  of `NTENS` (:jira:`685`, :pull:`47`). By `Nathan Miller`_.

Documentation
=============
- Add camelCase project name replacement instructions to project setup. By `Kyle Brindley`_.


******************
0.1.0 (2021-05-28)
******************

New Features
============
- Add CMake install configuration and CI/CD scripts for build, test, and installation to a Conda environment
  (:jira:`654`, :pull:`41`). By `Kyle Brindley`_.

Documentation
=============
- Update the Python package dependencies and add an example approach to future updates to the documentation
  (:jira:`636`, :pull:`37`). By `Kyle Brindley`_.
- Add file renaming commands to the project setup instructions (:jira:`634`, :pull:`38`). By `Kyle Brindley`_.
- Update the user manual to reflect required environment variable ``LD_LIBRARY_PATH`` (:jira:`662`, :pull:`43`). By
  `Kyle Brindley`_.

Internal Changes
================
- Update markdown syntax in README for wider compatibility (:jira:`604`, :pull:`36`). By `Kyle Brindley`_.
- Maintenance on ReST style guide updates (:jira:`604`, :pull:`36`). By `Kyle Brindley`_.
- Address BOOST output test stream deprecations and update minimum version
  (:jira:`654`, :pull:`41`). By `Kyle Brindley`_.
- Change project UMAT library name to avoid conflicts with external projects (:jira:`661`, :pull:`42`). By `Kyle
  Brindley`_.
- Remove the ``CXX`` compiler variable settings for build scripts (:jira:`671`, :pull:`44`). By `Kyle Brindley`_.

Enhancements
============
- Add multi-host and multi-environment CI/CD (:jira:`630`, :pull:`39`). By `Kyle Brindley`_.


******************
0.0.4 (2021-04-30)
******************

Documentation
=============
- Clarify behavior for custom target for the integration tests (:jira:`557`, :pull:`29`). By `Kyle Brindley`_.
- Add template documentation for the Abaqus material input definition (:jira:`575`, :pull:`31`). By `Kyle Brindley`_.
- Major overhaul of documentation organization to single source the Jenkins setup information from markdown files.  Adds
  the ``myst-parser`` Python package dependency and a pull request reviewer guide (:jira:`601`, :pull:`33`). By `Kyle
  Brindley`_.

Internal Changes
================
- Update Jenkins CI configuration to build and test for PRs to both ``master`` and ``dev`` branches (:jira:`544`,
  :pull:`26`). By `Kyle Brindley`_.
- Minor cleanup to root directory files. Move configuration and environment files to a subdirectory (:jira:`544`,
  :pull:`26`). By `Kyle Brindley`_.
- Add integration test CMake target for conditional rebuilds and file copy (:jira:`551`, :pull:`27`). By `Kyle
  Brindley`_.
- Add one ctest per abaqus input file (:jira:`551`, :pull:`27`). By `Kyle Brindley`_.
- Accept paths for input file in integration test shell script and check for errors in the abaqus stdout/stderr log
  (:jira:`551`, :pull:`27`). By `Kyle Brindley`_.
- Enable parallel CMake builds for continuous integration (CI) tests (:jira:`518`, :pull:`28`). By `Kyle Brindley`_.
- Add c++ source files ``*.cpp`` as dependencies for the Doxygen CMake target (:jira:`569`, :pull:`30`). By `Kyle
  Brindley`_.
- Add checks for ``STATEV`` and ``PROPS`` vector lengths to the abaqus interface. Throw exceptions with file and
  function name to interrupt Abaqus execution on input errors (:jira:`575`, :pull:`31`). By `Kyle Brindley`_.
- Add Abaqus interface unit tests for checking the ``STATEV`` and ``PROPS`` vector lengths (:jira:`575`, :pull:`31`). By
  `Kyle Brindley`_.
- Add unit tests for error codes in ``cpp_stub::sayHello`` (:jira:`334`, :pull:`32`). By `Kyle Brindley`_.

Enhancements
============
- Add error reporting to the Abaqus interface from the ``tardigrade_error_tools`` package (:jira:`334`, :pull:`32`). By `Kyle Brindley`_.


******************
0.0.3 (2021-04-13)
******************

Internal Changes
================
- Use ``tardigrade_abaqus_tools`` from a dedicated project (:jira:`535`, :pull:`23`). By `Kyle Brindley`_.
- Add ``bibtex_bibfiles`` variable to Sphinx configuration for newer version of ``sphinxcontrib.bibtex`` extension in
  Anaconda 2020 (:jira:`526`, :pull:`21`). By `Kyle Brindley`_.
- Add explicit list of documentation source files for better conditional CMake documentation re-builds (:jira:`526`,
  :pull:`21`). By `Kyle Brindley`_.


******************
0.0.2 (2021-02-11)
******************

Breaking changes
================
- Remove testing and support for intel ``icpc`` compiler (:jira:`516`, :pull:`9`). By `Kyle Brindley`_.

New Features
============
- Add do-nothing template c++ Abaqus UMAT interface and sample Abaqus input file (:jira:`502`, :pull:`6`). By `Kyle Brindley`_.
- Use example c++ library in Abaqus UMAT template (:jira:`505`, :pull:`8`). By `Kyle Brindley`_.
- Add c++ to fortran variable conversion and Abaqus variable return template (:jira:`521`, :pull:`15`, :pull:`16`). By
  `Kyle Brindley`_.
- Add common abaqus tensor handling tools and a c++ converted umat interface (:jira:`522`, :pull:`17`). By `Kyle
  Brindley`_.

Bug fixes
=========

Documentation
=============
- Add changelog to documentation (:jira:`450`, :pull:`11`). By `Kyle Brindley`_.
- Add direct CMake build instructions and minimal user manual (:jira:`519`, :pull:`12`). By `Kyle Brindley`_.
- Add release guidance and release branch instructions (:jira:`520`, :pull:`13`). By `Kyle Brindley`_.

Internal Changes
================
- Use BOOST and ctest for unit testing (:jira:`357`, :pull:`4`). By `Kyle Brindley`_.
- Update Jenkins CI configuration and store with version controlled repository (:jira:`442`, :pull:`5`). By `Kyle Brindley`_.
- Demonstrate c++ ``tardigrade_vector_tools`` library for unit testing (:jira:`506`, :pull:`7`). By `Kyle Brindley`_.
- Add integration tests for Abaqus UMAT interface (:jira:`504`, :pull:`10`). By `Kyle Brindley`_.
- Move project Abaqus interface into project files. Treat UMAT Fortran/c++ subroutine as a UMAT selection and pass
  through subroutine (:jira:`523`, :pull:`18`). By `Kyle Brindley`_.
- Bump micro version number for release (:jira:`524`). By `Kyle Brindley`_.

Enhancements
============


******************
0.0.1 (2020-10-26)
******************

Breaking changes
================

New Features
============
- Create c++ stub repository targeting constitutive modeling (:jira:`332`, :pull:`1`). By `Kyle Brindley`_.

Bug fixes
=========

Documentation
=============

Internal Changes
================
- Add continuous integration scripts (:jira:`333`, :pull:`2`). By `Kyle Brindley`_.

Enhancements
============
