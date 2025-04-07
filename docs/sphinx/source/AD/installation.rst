.. _ad_installation:

Installation
************

SICOPOLIS-AD v2 requires the installation of Tapenade as well as SICOPOLIS. It is mandatory to install the external libraries such as NetCDF, LIS to access the full functionality of the code, as well as git, to be able to clone and contribute to the repository.

.. _tapenade:

Open source AD Tool Tapenade
============================

`Tapenade <https://team.inria.fr/ecuador/tapenade/>`__ is an Automatic Differentiation Engine developed at `Inria at Sophia Antipolis <https://www.inria.fr/en/inria-centre-universite-cote-azur>`__ by the `Ecuador team <https://team.inria.fr/ecuador/>`__ (formerly `TROPICS team <https://www-sop.inria.fr/tropics/>`__). Tapenade takes as input a computer source program, plus a request for differentiation. Tapenade builds and returns the differentiated source program, that evaluates the required derivatives.

While the SICOPOLIS source files are prepared to generate adjoint sensitivities, they will not be able to do so without an operable installation of Tapenade. Fortunately the Tapenade installation procedure is straightforward.

We detail the instructions here, but the latest instructions can always be found `here. <http://www-sop.inria.fr/ecuador/tapenade/distrib/README.html>`__

Prerequisites for Linux or Mac OS X
-----------------------------------

Before installing Tapenade, you must check that an up-to-date Java Runtime Environment is installed. Tapenade will not run with older Java Runtime Environment.

Steps for Mac OS
----------------

Tapenade 3.16 distribution does not contain a fortranParser executable for MacOS. It uses a docker image from `here <https://gitlab.inria.fr/tapenade/tapenade>`__. You need docker on your Mac to run the Tapenade distribution with Fortran programs. Details on how to build fortranParser is `here <https://tapenade.gitlabpages.inria.fr/tapenade/docs/html/src/frontf/README.html?highlight=mac>`__. You may also build Tapenade on your Mac from the `gitlab repository <https://tapenade.gitlabpages.inria.fr/tapenade/docs/html/distrib/README.html>`__.

Steps for Linux
---------------

1. Read `the Tapenade license. <https://tapenade.gitlabpages.inria.fr/userdoc/build/html/LICENSE.html>`__

2. Download `tapenade_3.16.tar <https://tapenade.gitlabpages.inria.fr/tapenade/distrib/tapenade_3.16.tar>`__ into your chosen installation directory *install_dir*.

  **NOTE**: Alternatively, a Tapenade version that works correctly with SICOPOLIS-AD v2 is always available in the ``test_ad/tapenade_supported`` directory.

3. Go to your chosen installation directory *install_dir*, and extract Tapenade from the tar file::

    % tar xvfz tapenade_3.16.tar

4. On Linux, depending on your distribution, Tapenade may require you to set the shell variable ``JAVA_HOME`` to your java installation directory. It is often ``JAVA_HOME=/usr/java/default``. You might also need to modify the ``PATH`` by adding the bin directory from the Tapenade installation. An example can be found :ref:`here <tapenade_bashrc_snippet>`.

.. _tapenade_bashrc_snippet:

**NOTE**: Every time you wish to use the adjoint capability of SICOPOLIS-AD, you must re-source the environment. We recommend that this be done automatically in your bash or c-shell profile upon login. An example of an addition to a ``.bashrc`` file from a Linux server is given below. Luckily, shell variable ``JAVA_HOME`` was not required to be explicitly set for this particular Linux distribution, but might be necessary for some other distributions. ::

  ##set some env variables for SICOPOLIS tapenade

  export TAPENADE_HOME="/home/shreyas/tapenade_3.16"
  export PATH="$PATH:$TAPENADE_HOME/bin"

  ##Modules

  module use /share/modulefiles/
  module load java/jdk/16.0.1 # Java required by Tapenade

You should now have a working copy of Tapenade.

For more information on the tapenade command and its arguments, type::

  tapenade -?

Prerequisites for Windows
-------------------------

**NOTE**: Although Tapenade can be built on Windows, SICOPOLIS requires a Unix-like system (e.g., Linux), as mentioned in ":ref:`getting_started`".

Before installing Tapenade, you must check that an up-to-date Java Runtime Environment is installed. Tapenade will not run with older Java Runtime Environment. The Fortran parser of Tapenade uses `cygwin <https://www.cygwin.com/>`__.

Steps for Windows
-----------------

1. Read `the Tapenade license. <https://tapenade.gitlabpages.inria.fr/userdoc/build/html/LICENSE.html>`__

2. Download `tapenade_3.16.zip <https://tapenade.gitlabpages.inria.fr/tapenade/distrib/tapenade_3.16.zip>`__ into your chosen installation directory *install_dir*.

3. Go to your chosen installation directory *install_dir*, and extract Tapenade from the zip file.

4. Save a copy of the ``install_dir\tapenade_3.16\bin\tapenade.bat`` file and modify ``install_dir\tapenade_3.16\bin\tapenade.bat`` according to your installation parameters\:

  | replace ``TAPENADE_HOME=..`` by ``TAPENADE_HOME="install_dir"\tapenade_3.16``
  | replace ``JAVA_HOME="C:\Progra~1\Java\jdkXXXX"`` by your current java directory
  | replace ``BROWSER="C:\Program Files\Internet Explorer\iexplore.exe"`` by your current browser.

.. _sicopolis_ad_config:

SICOPOLIS-AD v2
===============

See the ":ref:`getting_started`" section. However, we will mention here some steps since using the Automatic Differentiation (AD) capabilities with Tapenade requires a slightly modified setup.

The satisfaction of the following prerequisites is highly recommended to access all the features of the code. Details can differ from the ":ref:`getting_started`" section, since there are multiple ways to do things. We detail one of them here.

GNU GCC Compiler (gfortran+gcc) or Intel Compiler (ifort+icc)
-------------------------------------------------------------

We have tested SICOPOLIS-AD v2 on gfortran/gcc v5.4.0, v7.2.0 and v8.5.0, any intermediate versions should work just as well. We have also tested it on ifort/icc v18.0.0 (however, it should be noted that we have not tested the external Lis solver with Intel compilers).

Lis (1.4.43 or newer)
---------------------

Lis installation is mandatory to use shallow-shelf/shelfy-stream dynamics in simulations. Install Lis as explained in :ref:`Dependencies/Lis <dependencies-lis>`.

The following commands might be helpful, they are written for the latest version at the time of writing::

  wget https://www.ssisc.org/lis/dl/lis-2.0.30.zip
  unzip lis-2.0.30.zip
  cd lis-2.0.30
  ./configure --enable-fortran --prefix=/lis-2.0.30/installation/ --enable-shared && make && make check && make install

For AD purposes, we compile the code using the ``src/MakefileTapenade`` makefile. This makefile requires the following environment variables set\:

1. ``LISDIR`` - The installation directory for the Lis version to be used. You can change this variable, either in the Makefile directly, or automatically in your bash or c-shell profile upon login (for example, ``.bashrc``). Examples for both are shown here.

  ``src/MakefileTapenade``::

    ifndef LISDIR
    LISDIR=/home/shreyas/lis-2.0.30/installation
    endif

  ``.bashrc``::

    export LISDIR="/home/shreyas/lis-2.0.30/installation"

2. ``LIBLIS`` - Absolute path to ``liblis.so``. By default in ``src/MakefileTapenade``, it is ``${LISDIR}/lib/liblis.so``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``. 

3. ``LIBLISFLAG`` - Add directory using ``-L`` to be searched for ``-llis``. By default in ``src/MakefileTapenade``, it is ``-L${LISDIR}/lib -llis``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``.

4. ``LISFLAG`` - This flag declares directories to be searched for Lis ``#include`` header file ``lisf.h``, as well as defines the ``BUILD_LIS`` as a macro with value 1. By default in ``src/MakefileTapenade``, it is ``-DBUILD_LIS -I${LISDIR}/include/``. If you follow the original instructions to install Lis, this should work, else one can set it manually within ``src/MakefileTapenade``.

**NOTE**: Some users have reported needing to extend their ``LD_LIBRARY_PATH`` with the location of ``${LISDIR}/lib`` in order to find ``liblis.so.0``.

NetCDF (3.6.x or newer)
-----------------------

NetCDF installation is mandatory since it is a powerful library with widespread use for I/O with a machine-independent data format. Install NetCDF as explained in :ref:`Dependencies/NetCDF <dependencies-netcdf>`.

In some cases, for example while working on a shared server which uses a module manager or Docker container, thing have to be set up differently. ``src/MakefileTapenade`` needs either the ``NETCDF_FORTRAN_DIR`` macro set or both ``NETCDF_F90_FLAG`` and ``LIB_NETCDF_F90_FLAG`` set (see code snippet from ``src/MakefileTapenade`` here). ::

  ifndef NETCDF_F90_FLAG
  ifndef LIB_NETCDF_F90_FLAG
  ifdef NETCDF_FORTRAN_DIR
  LIB_NETCDF_F90=${NETCDF_FORTRAN_DIR}/lib/libnetcdff.so
  LIB_NETCDF_F90_FLAG=-L${NETCDF_FORTRAN_DIR}/lib -lnetcdff
  NETCDF_F90_FLAG=-I${NETCDF_FORTRAN_DIR}/include/
  endif
  endif
  endif

1. ``NETCDF_FORTRAN_DIR`` - The installation directory for netcdf-fortran. You can change this variable, either in the Makefile directly, or automatically in your bash or c-shell profile upon login (for example, ``.bashrc``). Examples for both are shown here.

  ``src/MakefileTapenade`` ::

    ifndef NETCDF_FORTRAN_DIR
    NETCDF_FORTRAN_DIR=/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4
    endif

  ``.bashrc`` ::

    export NETCDF_FORTRAN_DIR="/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4"

2. ``LIB_NETCDF_F90`` - Absolute path to ``libnetcdff.so``. By default in ``src/MakefileTapenade``, it is ``{NETCDF_FORTRAN_DIR}/lib/libnetcdff.so``.

3. ``LIB_NETCDF_F90_FLAG`` - Add directory using ``-L`` to be searched for ``-lnetcdff``. By default in ``src/MakefileTapenade``, it is ``-L${NETCDF_FORTRAN_DIR}/lib -lnetcdff``. See some examples below where this has to be set explicitly in case of a Docker container.

4. ``NETCDF_F90_FLAG`` - This flag declares directories to be searched for netcdf-fortran ``#include`` header files. By default in ``src/MakefileTapenade``, it is ``-I${NETCDF_FORTRAN_DIR}/include/``. See some examples below where this has to be set explicitly in case of a Docker container. 

For a server that uses modules, you can load the relevant modules using commands like these (can also make permanent by adding to login script like ``.bashrc``::

  % module use /share/modulefiles/
  % module load openmpi
  % module load netcdf-fortran
  % module load netcdf

You then have to give the ``NETCDF_FORTRAN_DIR`` macro to ``src/MakefileTapenade``, either by adding to a login script or directly inside the makefile. If your system uses a module manager, you can query to find the exact directory location. ::

  % module show netcdf-fortran
  ----------------------------------------------------------
  /opt/ohpc/pub/moduledeps/gnu-openmpi/netcdf-fortran/4.4.4:
  ----------------------------------------------------------

  whatis("Name: NETCDF_FORTRAN built with gnu toolchain ")
  whatis("Version: 4.4.4 ")
  whatis("Category: runtime library ")
  whatis("Description: Fortran Libraries for the Unidata network Common Data Form ")
  whatis("http://www.unidata.ucar.edu/software/netcdf/ ")
  prepend_path("PATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/bin")
  prepend_path("MANPATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/share/man")
  prepend_path("INCLUDE","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/include")
  prepend_path("LD_LIBRARY_PATH","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/lib")
  setenv("NETCDF_FORTRAN_DIR","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4")
  setenv("NETCDF_FORTRAN_BIN","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/bin")
  setenv("NETCDF_FORTRAN_LIB","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/lib")
  setenv("NETCDF_FORTRAN_INC","/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4/include")
  help([[ 
  This module loads the NetCDF Fortran API built with the gnu compiler toolchain.
   
  Note that this build of NetCDF leverages the HDF I/O library and requires linkage
  against hdf5 and the native C NetCDF library. Consequently, phdf5 and the standard C
  version of NetCDF are loaded automatically via this module. A typical compilation
  example for Fortran applications requiring NetCDF is as follows:
   
  ]])

In this case ``NETCDF_FORTRAN_DIR=/opt/ohpc/pub/libs/gnu/openmpi/netcdf-fortran/4.4.4``.

For a Docker container, for example one with a ``centos:8`` distribution, and the ``dnf`` package manager, NetCDF is typically installed as follows::

  RUN dnf install -y https://github.com/openhpc/ohpc/releases/download/v1.3.GA/ohpc-release-1.3-1.el7.x86_64.rpm
  
  # Add some packages
  RUN dnf -y install epel-release
  RUN dnf -y install dnf-plugins-core
  RUN dnf config-manager --set-enabled powertools
  RUN dnf -y install make which git
  RUN dnf -y install diffutils
  RUN dnf -y install vim
  RUN dnf -y install autoconf automake
  RUN dnf -y install valgrind-ohpc
  RUN dnf -y install gnu8-compilers-ohpc
  RUN dnf -y install gsl-gnu8-ohpc hdf5-gnu8-ohpc
  RUN dnf -y install openmpi-devel
  RUN dnf -y install bc wget zlib-devel perl-Digest-MD5
  RUN dnf -y --enablerepo=powertools install netcdf-fortran netcdf-devel # NetCDF installation
  RUN dnf -y install netcdf-fortran-devel # NetCDF installation

In this case, you will find that the ``./usr/lib64/gfortran/modules/netcdf.mod`` exists in your docker environment. In this case, you can just directly set  ``NETCDF_F90_FLAG=-I/usr/lib64/gfortran/modules`` either the makefile or the login script (no need to set ``NETCDF_FORTRAN_DIR`` macro). 

You can also confirm that the files ``/usr/lib64/libnetcdff.so*`` and ``/usr/lib64/libnetcdf.so*`` exist, which means you have to set ``LIB_NETCDF_F90_FLAG=-L/usr/lib64 -lnetcdff``.

The instructions given in :ref:`Dependencies/NetCDF <dependencies-netcdf>`, and these two cases should help cover most of the issues with the installation of NetCDF.

Downloading SICOPOLIS-AD v2
---------------------------

As described in the ":ref:`getting_started`" section. However, when using Git, the ``ad`` branch should be cloned::

  git clone --branch ad \
      https://github.com/sicopolis/sicopolis.git

Tagged versions of SICOPOLIS-AD are also available from `Zenodo <https://doi.org/10.5281/zenodo.3686392>`__.

Initial configuration
---------------------

As described in the ":ref:`getting_started`" section.

Now you should be ready to use SICOPOLIS-AD v2, as described in :ref:`Running SICOPOLIS-AD v2 <ad_running>`.
