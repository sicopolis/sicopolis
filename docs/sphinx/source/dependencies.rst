.. _dependencies:

Dependencies
************

.. _dependencies-netcdf:

NetCDF
======

NetCDF (Network Common Data Form) is a common format for scientific data (Unidata :cite:`unidata_2023`) that is also used by SICOPOLIS. The NetCDF C and Fortran libraries are required.

For **GCC**, installation from a package manager is recommended. Under openSUSE Leap 15.3, install ``netcdf``, ``netcdf-devel``, ``netcdf-devel-static``, ``netcdf-fortran``, ``netcdf-fortran-devel``, ``netcdf-fortran-static``, ``ncview``. This requires the repositories "Software for Scientists and Engineers" and "sebschub's Home Project". Details (especially the required repositories) will differ for other systems.

For the **Intel compiler**, building from the source code is required. The C and Fortran sources are available for download on the `NetCDF website <https://doi.org/10.5065/D6H70CW6>`__ as zip or tar archives. Unzip the archives into temporary directories.

* Prior to version 4.2, a single archive contained both the C and Fortran libraries. A minimal installation for version 4.1.3 (without NetCDF-4 support) can be done by changing to the source directory, then::

    export NCDIR=/opt/netcdf
    export CC=icc
    export FC=ifort
    export CFLAGS="-O2"
    export CPPFLAGS=
    export FCFLAGS="-O2"
    export FFLAGS=${FCFLAGS}
    ./configure --prefix=${NCDIR} --libdir=${NCDIR}/lib \
                --disable-netcdf-4
    make install

* Since version 4.2, the C and Fortran libraries must be installed separately. If NetCDF-4 support is dispensable, the following installation should work (tested under openSUSE Leap 15.2 and icc/ifort 19.1 with versions netcdf-c-4.8.0 and netcdf-fortran-4.5.3).

  Change to the source directory of the C library, then::

    export NCDIR=/opt/netcdf
    export CC=icc
    export FC=ifort
    export CFLAGS="-O2"
    export CPPFLAGS=
    export FCFLAGS="-O2"
    export FFLAGS=${FCFLAGS}
    ./configure --prefix=${NCDIR} --libdir=${NCDIR}/lib \
                --disable-netcdf-4 --enable-logging
    make install

  Change to the source directory of the Fortran library, then::

    export NFDIR=/opt/netcdf
    export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
    export CPPFLAGS=-I${NCDIR}/include
    export LDFLAGS=-L${NCDIR}/lib
    ./configure --prefix=${NFDIR} --libdir=${NFDIR}/lib \
                --disable-netcdf-4 --enable-logging
    make install

* The recommended, complete build with NetCDF-4 support requires additional libraries. Download the sources of the data-compression library **zlib** (https://www.zlib.net/), the Hierarchical Data Format **HDF5** (https://www.hdfgroup.org/) and the **NetCDF-4** C and Fortran libraries (https://doi.org/10.5065/D6H70CW6). Unpack all archives.

  Step 1: Installation of **zlib**. Change to the source directory, then::

    export ZDIR=/opt/zlib
    export CC=icc
    export FC=ifort
    export CFLAGS="-O2"
    export CPPFLAGS=
    export FCFLAGS="-O2"
    export FFLAGS=${FCFLAGS}
    ./configure --prefix=${ZDIR} --libdir=${ZDIR}/lib
    make check
    make install

  Step 2: Installation of **HDF5**. Change to the source directory, then::

    export ZDIR=/opt/zlib
    export H5DIR=/opt/hdf5
    export LD_LIBRARY_PATH=${ZDIR}/lib:${LD_LIBRARY_PATH}
    export CC=icc
    export FC=ifort
    export CFLAGS="-O2"
    export CPPFLAGS=
    export FCFLAGS="-O2"
    export FFLAGS=${FCFLAGS}
    ./configure --with-zlib=${ZDIR} \
                --prefix=${H5DIR} --libdir=${H5DIR}/lib \
                --enable-hl \
                --with-default-api-version=v18
    make check
    make install 

  Step 3: Installation of **NetCDF-4**. Change to the source directory of the C library, then::

    export ZDIR=/opt/zlib
    export H5DIR=/opt/hdf5
    export NCDIR=/opt/netcdf
    export LD_LIBRARY_PATH=${H5DIR}/lib:${LD_LIBRARY_PATH}
    export CC=icc
    export FC=ifort
    export CFLAGS="-O2"
    export CPPFLAGS="-I${H5DIR}/include -I${ZDIR}/include"
    export LDFLAGS="-L${H5DIR}/lib -L${ZDIR}/lib"
    export FCFLAGS="-O2"
    export FFLAGS=${FCFLAGS}
    ./configure --prefix=${NCDIR} --libdir=${NCDIR}/lib \
                --enable-logging --disable-dap-remote-tests
    make install
    make check

  Change to the source directory of the Fortran library, then::

    export NFDIR=/opt/netcdf
    export LD_LIBRARY_PATH=${NCDIR}/lib:${LD_LIBRARY_PATH}
    export CPPFLAGS=-I${NCDIR}/include
    export LDFLAGS=-L${NCDIR}/lib
    ./configure --prefix=${NFDIR} --libdir=${NFDIR}/lib \
                --enable-logging --disable-dap-remote-tests
    make install
    make check

  This was tested under openSUSE Leap 15.3 and icc/ifort 19.1 with versions zlib-1.2.11, hdf5-1.12.1, netcdf-c-4.8.1 and netcdf-fortran-4.5.4.

  If installation from a package manager does not work out for GCC, try the above procedure, replacing everywhere ``icc`` and ``ifort`` by ``gcc`` and ``gfortran``, respectively.

If needed, further instructions can be found on the `NetCDF website <nc>`__.

Installation under ``/opt`` usually requires admin rights. The same holds for the common alternative ``/usr/local``. For a local installation, replace it by ``/home/<my_user_name>/local``.

.. _dependencies-lis:

Lis
===

Lis (Library of Iterative Solvers for linear systems) is a software library for solving discretized linear equations (Nishida :cite:`nishida_2010`).

Download the source of Lis as a zip archive from https://www.ssisc.org/lis/. Unzip the archive and change to the source directory.

For **GCC**, install lis by executing::

  export LISDIR=/opt/lis
  ./configure --prefix=${LISDIR} --libdir=${LISDIR}/lib \
              --enable-fortran --enable-f90 \
              --enable-omp --enable-saamg --enable-fma \
              CC=gcc FC=gfortran F77=gfortran \
              CFLAGS="-mcmodel=medium" CPPFLAGS="-mcmodel=medium" \
              FCFLAGS="-mcmodel=medium" FFLAGS="-mcmodel=medium"
  make install

This was tested under openSUSE Leap 15.2 and Linux Mint 20.1 with lis-2.0.30. Some modifications might be needed under different systems.

For the **Intel compiler**, replace ``gcc`` and ``gfortran`` by ``icc`` and ``ifort``, respectively.

Installation under ``/opt`` usually requires admin rights. The same holds for the common alternative ``/usr/local``. For a local installation, replace it by ``/home/<my_user_name>/local``.
