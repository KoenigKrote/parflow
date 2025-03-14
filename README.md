# ParFlow

[![Build Status](https://travis-ci.org/parflow/parflow.svg?branch=master)](https://travis-ci.org/parflow/parflow)

ParFlow is an open-source, modular, parallel watershed flow model. It
includes fully-integrated overland flow, the ability to simulate
complex topography, geology and heterogeneity and coupled land-surface
processes including the land-energy budget, biogeochemistry and snow
(via CLM). It is multi-platform and runs with a common I/O structure
from laptop to supercomputer. ParFlow is the result of a long,
multi-institutional development history and is now a collaborative
effort between CSM, LLNL, UniBonn and UCB. ParFlow has been coupled to
the mesoscale, meteorological code ARPS and the NCAR code WRF.

The Parflow User Manual is available at [Parflow Users
Manual](https://github.com/parflow/parflow/blob/master/parflow-manual.pdf).
The manual contains additional documentation on how to use ParFlow and
setup input files.  A quick start is included below.

### Citing Parflow

To cite Parflow, please use the following reference.

If you use ParFlow in a publication, please cite the these papers that describe model physics:

* Ashby S.F. and R.D. Falgout, Nuclear Science and Engineering 124:145-159, 1996
* Jones, J.E. and C.S. Woodward, Advances in Water Resources 24:763-774, 2001
* Kollet, S.J. and R.M. Maxwell, Advances in Water Resources 29:945-958, 2006
* Maxwell, R.M. Advances in Water Resources 53:109-117, 2013

If you use ParFlow coupled to CLM in a publication, please also cite
two additional papers that describe the coupled model physics:

* Maxwell, R.M. and N.L. Miller, Journal of Hydrometeorology 6(3):233-247, 2005
* Kollet, S.J. and R.M. Maxwell, Water Resources Research 44:W02402, 2008

### Additional Parflow resources

The ParFlow website has additional information on the project:
- [Parflow Web Site](https://parflow.org/)

You can join the Parflow mailing list to communicate with the Parflow
developers and users.  Join our mailing list over via:
- [Parflow-Users](https://mailman.mines.edu/mailman/listinfo/parflow-users)

A Parflow blog is available with notes from users on how to compile and use Parflow:
- [Parflow Blog](http://parflow.blogspot.com/)

To report Parflow bugs, please use the GitHub issue tracker for Parflow:
- [Parflow Issue Tracker](https://github.com/parflow/parflow/issues)

## Quick Start on Unix/Linux

Important note for users that have built with Autoconf, the CMake
configure process is one step by default.  Most builds of of ParFlow
are on MPP architectures or workstations where the login node and
compute nodes are same architecture the default build process builds
both the ParFlow executable and tools with the same compilers and
libraries in one step.  This will hopefully make building easier for
the majority of users.  It is still possible to build the two
components separately; see instruction below for building pftools and
pfsimulator separately.

CMake supports builds for several operating systems and IDE tools
(like Visual Studio on Windows and XCode on MacOS).  The ParFlow team
has not tested building on platforms other than Linux; there will
likely be some issues on other platforms.  The ParFlow team welcomes
bug reports and patches if you attempt other builds.

### Step 1: Setup

Decide where to install ParFlow and associated libraries.

Set the environment variable `PARFLOW_DIR` to the chosen location:

For bash:

```shell
   export PARFLOW_DIR=/home/snoopy/parflow
```   

For csh and tcsh:

```shell
   setenv PARFLOW_DIR /home/snoopy/parflow
```

### Step 2: Extract the Source

Extract the source files from the compressed tar file.

Obtain the release from the ParFlow GitHub web site:

https://github.com/parflow/parflow/releases

and extract the release.  Here we assume you are building in new
subdirectory in your home directory:

```shell
   mkdir ~/parflow 
   cd ~/parflow 
   tar -xvf ../parflow.tar.gz
```

Note the ParFlow tar file will be have a different name based on the
version number.

If you are not using GNU tar or have a very old version GNU tar you
will need to uncompress the file first:

```shell
   mkdir ~/parflow 
   cd ~/parflow 
   gunzip ../parflow.tar.gz
   tar -xvf ../parflow.tar
```

### Step 3: Running CMake to configure ParFlow

CMake is a utility that sets up makefiles for building ParFlow.  CMake
allows setting of compiler to use and other options.  First create a
directory for the build.  It is generally recommend to build outside
of the source directory to make it keep things clean.  For example,
restarting a failed build with a separate build directory simply
involves removing the build directory.

#### Building with the ccmake GUI

You can control build options for ParFlow using the ccmake GUI.

```shell
   mkdir build
   cd build
   ccmake ../parflow 
```
At a minimum, you will want to set the CMAKE_INSTALL_PREFIX value to the same thing
as PARFLOW_DIR was set to above.  Other variables should be set as desired.

After setting a variable 'c' will configure `ParFlow.  When you are
completely done setting configuration options, use 'g' to generate the
configuration and exit ccmake.

If you are new to CMake, the creators of CMake provide some additional ccmake usage notes here:

https://cmake.org/runningcmake/

#### Building with the cmake command line

CMake may also be configured from the command line using the cmake
command.  The default will configure a sequential version of ParFlow
using MPI libraries.  CLM is being enabled.

```shell
   mkdir build
   cd build
   cmake ../parflow \
   	 -DCMAKE_INSTALL_PREFIX=$(PARFLOW_DIR) \
   	 -DPARFLOW_HAVE_CLM=ON
```

If TCL is not installed in the standard locations (/usr or /usr/local)
you need to specify the path to the tclsh location:

```shell
	-DTCL_TCLSH=${PARFLOW_TCL_DIR}/bin/tclsh8.6
```

Building a parallel version of ParFlow requires the communications
layer to use must be set.  The most common option will be MPI.  Here
is a minimal example of an MPI build with CLM:

```shell
   mkdir build
   cd build
   cmake ../parflow \
      	 -DCMAKE_INSTALL_PREFIX=$(PARFLOW_DIR) \
   	 -DPARFLOW_HAVE_CLM=ON \
	 -DPARFLOW_AMPS_LAYER=mpi1
```

Here is a more complex example where location of various external
packages are being specified and some features are being enabled:

```shell
   mkdir build
   cd build
   cmake ../parflow \
        -DPARFLOW_AMPS_LAYER=mpi1 \
	-DHYPRE_ROOT=$(PARFLOW_HYPRE_DIR) \
	-DHDF5_ROOT=$(PARFLOW_HDF5_DIR) \
	-DSILO_ROOT=$(PARFLOW_SILO_DIR) \
	-DSUNDIALS_ROOT=$(PARFLOW_SUNDIALS_DIR) \
	-DCMAKE_BUILD_TYPE=Debug \
	-DPARFLOW_ENABLE_TIMING=TRUE \
	-DPARFLOW_HAVE_CLM=ON \
	-DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR)
```

### Step 4: Building and installing

Once CMake has configured and created a set of Makefiles; building is
easy:

```shell
   cd build
   make 
   make install
```

### Step 5: Running a sample problem

If all went well a sample ParFlow problem can be run using:

```shell
cd parflow/test
tclsh default_single.tcl 1 1 1
```

Note that the environment variable `PAFLOW_DIR` must be set for this
to work and it assumes tclsh is in your path.  Make sure to use the
same TCL shell as was used in the cmake configure.

Some parallel machines do not allow launching a parallel executable
from the login node; you may need to run this command in a batch file
or by starting a parallel interactive session.

## Building documentation

### User Manual

A version of the user manual is available at github : [Parflow Users Manual](https://github.com/parflow/parflow/blob/master/parflow-manual.pdf)

The user manual for Parflow may be built as part of the build when
Latex is available on the system. Adding the
-DPARFLOW_ENABLE_LATEX=TRUE option to the CMake configure will enable
building of the documentation.

```shell
   mkdir build
   cd build
   cmake ../parflow \
        <other cmake options> \
	-DPARFLOW_ENABLE_LATEX=TRUE \
	-DCMAKE_INSTALL_PREFIX=$(INSTALL_DIR)
```

When make is run the documenation will be built and installed in
$(INSTALL_DIR)/docs/user_manual.pdf.

### Code documentation

Parflow is moving to using Doxygen for code documenation.  The documentation is currently very sparse.

Adding the -DPARFLOW_ENABLE_DOXYGEN=TRUE option to the CMake configure
will enable building of the code documentation.  After CMake has been
run the Doxygen code documenation is built with:

```shell
   cd build
   make doxygen
```

HTML pages are generated in the build/docs/doxygen/html directory.

## Configure options

A number of packages are optional for building ParFlow.  The optional
packages are enabled by PARFLOW_ENABLE_<package> value to be `TRUE` or
setting the <package>_ROOT=<directory> value.  If a package is enabled
with the using an ENABLE flag CMake will attempt to find the package
in standard locations.  Explicitly setting the location using the ROOT
variable for a package automatically enables it, you don't need to
specify both values.

### How to specify command to run MPI applications

There are multiple ways to run MPI applications such as mpiexec,
mpirun, srun, and aprun.  The command used is dependent on the job
submission system used.  By default CMake will attempt to determine an
appropriate tool; a process that does not always yield the correct result.

There are several ways to modify the CMake guess on how applications
should be run.  At configure time you may overwride the MPI launcher
using:

```shell 
   -DMPIEXEC="<launcher-name>"
   -DMPIEXEC_NUMPROC_FLAG="<flag used to set number of tasks>"
```

An example for mpiexec is -DMPIEXEC="mpiexec" -DMPIEXEC_NUMPROC_FLAG="-n".

The ParFlow script to run MPI applications will also include options
specified in the environment variable PARFLOW_MPIEXEC_EXTRA_FLAGS on
the MPI execution command line.  For example when running with OpenMPI
on a single workstation the following will enable running more MPI
tasks than cores and disable the busy loop waiting to improve
performance:

```shell
   export PARFLOW_MPIEXEC_EXTRA_FLAGS="--mca mpi_yield_when_idle 1 --oversubscribe"
```

Last the TCL script can explicity set the command to invoke for
running ParFlow.  This is done by setting the Process.Command key in
the input database.  For example to use the mpiexec command and
control the cpu set used the following command string can be used:

```shell
   pfset Process.Command "mpiexec -cpu-set 1 -n %d parflow %s"
```

The '%d' will be replaced with the number of processes (computed using
the Process.Topology values : P * Q * R) and the '%s' will be replaced
by the name supplied to the pfrun command for the input database name.
The following shows how the default_single.tcl script could be
modified to use the custom command string:

```shell
   pfset Process.Command "mpiexec -cpu-set 1 -n %d parflow %s"
   pfrun default_single
   pfundist default_single
```
## Building simulator and tools support separately

This section is for advanced users runing on heterogenous HPC architectures.

ParFlow is composed of two main components that maybe configured and
built separately.  Some HPC platforms are heterogeneous with the login
node being different than the compute nodes.  The ParFlow system has
an executable for the simulator which needs to run on the compute
nodes and a set of TCL libraries used for problem setup that can be
run on the login node for problem setup.

The CMake variables PARFLOW_ENABLE_SIMULATOR and PARFLOW_ENABLE_TOOLS
control which component is configured.  By default both are `TRUE`.  To
build separately use two build directories and run cmake in each to
build the simulator and tools components separately. By specifying
different compilers and options for each, one can target different
architectures for each component.

## Release

Copyright (c) 1995-2019, Lawrence Livermore National Security LLC. 

Produced at the Lawrence Livermore National Laboratory. 

Written by the Parflow Team (see the CONTRIBUTORS file)

CODE-OCEC-08-103. All rights reserved.

Parflow is released under the GNU General Public License version 2.1

For details and restrictions, please read the LICENSE.txt file.
- [LICENSE](./LICENSE.txt)
