#-------------------------------------------------------------------------------
# Configuration file for SICOPOLIS
#-------------------------------------------------------------------------------

#-------- Compiler --------

export FC=gfortran
###    Can be set here if needed.
###    So far, gfortran and ifort are supported.

if [ "$FC" != "" ] ; then
   echo "Fortran compiler >$FC< found."
else
   echo "No Fortran compiler found (environment variable FC not set)."
   echo "Trying gfortran..."
   FC=gfortran
fi

#-------- Flags --------

LIS_FLAG="true"
#   "true"  -> linking of SICOPOLIS with Lis library
#              (requires also OPENMP_FLAG="true")
#   "false" -> linking of SICOPOLIS without Lis library

OPENMP_FLAG="true"
#   "true"  -> OpenMP will be used
#              (required if LIS_FLAG="true")
#   "false" -> OpenMP will not be used

LARGE_DATA_FLAG="false"
#   "true"  -> no memory restriction on data (more than 2GB possible);
#              required for very-high-resolution simulations
#   "false" -> memory restriction on data (max. 2GB possible);
#              OK for medium- and low-resolution simulations

#-------- NetCDF settings --------

# module load ...
### A 'module load' command is needed on some systems.

# module load netcdf/4.4.1_gnu5.4.0
### needed for GEO-UiO Linux environment (as of Aug 2020)

export NETCDFHOME=/opt/netcdf

# export NETCDFHOME=/usr
#        ### often works if NetCDF was installed from a repository
#        ### rather than manually

# export NETCDFHOME=/opt/uio/modules/packages/netcdf/4.4.1_gnu5.4.0
#        ### setting for GEO-UiO Linux environment (as of Aug 2020)

if [ -z "${LD_LIBRARY_PATH}" ]; then
   export LD_LIBRARY_PATH=${NETCDFHOME}/lib
else
   export LD_LIBRARY_PATH=${NETCDFHOME}/lib\:${LD_LIBRARY_PATH}
fi
   ### !!! On some systems, lib may have to be changed to lib64 !!!

if [ -f $NETCDFHOME/bin/nf-config ] ; then
   # Automatic settings with nf-config
   if [ "$($NETCDFHOME/bin/nf-config --has-f90)" = "yes" ] || [ "$($NETCDFHOME/bin/nf-config --has-f03)" = "yes" ]; then
      NETCDF_FLAGS="$($NETCDFHOME/bin/nf-config --fflags)"
      NETCDF_LIBS="$($NETCDFHOME/bin/nf-config --flibs)"
   else
      echo "NetCDF error: compiled without either f90 or f03 interface."
      exit 1
   fi
elif [ -f $NETCDFHOME/bin/nc-config ] ; then
   # Automatic settings with nc-config
   if [ "$($NETCDFHOME/bin/nc-config --has-f90)" = "yes" ] ; then
      NETCDF_FLAGS="$($NETCDFHOME/bin/nc-config --cflags)"
      NETCDF_LIBS="$($NETCDFHOME/bin/nc-config --flibs)"
   else
      echo "NetCDF error: compiled without f90 interface."
      exit 1
   fi
else
   # Manual settings
   NETCDFINCLUDE=${NETCDFHOME}'/include'
   NETCDFLIB=${NETCDFHOME}'/lib'
fi

#-------- Lis settings --------

# module load ...
### A 'module load' command is needed on some systems.

if [ "$LIS_FLAG" = "true" ] ; then
   LISHOME=/opt/lis
   LISINCLUDE=${LISHOME}'/include'
   LISLIB=${LISHOME}'/lib'
fi

#-------- Compiler flags --------

if [ "$FC" = "ifort" ] ; then

   case $PROGNAME in
        "sicopolis")
           FCFLAGS='-xHOST -O3 -no-prec-div'
           # This is '-fast' without '-static' and '-ipo'
           ;;
        "sicograph")
           FCFLAGS='-O2 -Vaxlib'
           ;;
        *)
           FCFLAGS='-O2'
           ;;
   esac            

   if [ "$LARGE_DATA_FLAG" = "true" ] ; then
      FCFLAGS=${FCFLAGS}' -mcmodel=medium -shared-intel'
      # No memory restriction on data (more than 2GB possible)
   fi

   if [ "$OPENMP_FLAG" = "true" ] ; then
      case $PROGNAME in
           "sicopolis")
              FCFLAGS=${FCFLAGS}' -qopenmp'
              # Change to ' -openmp' for older versions of the compiler
              ;;
           *) ;;
      esac            
   fi

elif [ "$FC" = "gfortran" ] ; then

   case $PROGNAME in
        "sicopolis")
           FCFLAGS='-O2 -ffree-line-length-none'
           ;;
        "sicograph")
           FCFLAGS='-O2 -ffree-line-length-none'
           ;;
        *)
           FCFLAGS='-O2 -ffree-line-length-none'
           ;;
   esac            

   if [ "$LARGE_DATA_FLAG" = "true" ] ; then
      FCFLAGS=${FCFLAGS}' -mcmodel=medium'
      # No memory restriction on data (more than 2GB possible)
   fi

   if [ "$OPENMP_FLAG" = "true" ] ; then
      case $PROGNAME in
           "sicopolis")
              FCFLAGS=${FCFLAGS}' -fopenmp'
              ;;
           *) ;;
      esac            
   fi

else
   echo "Unknown compiler flags for >$FC<, must exit."
   echo "Add flags to >sico_configs.sh< and try again."
   exit 1
fi

if [ "$LIS_FLAG" = "true" ] ; then
   case $PROGNAME in
        "sicopolis")
           FCFLAGS=${FCFLAGS}' -I'${LISINCLUDE}' -L'${LISLIB}' -llis'
           ;;
        *) ;;
   esac            
fi

if [ -f $NETCDFHOME/bin/nf-config ] ; then
   # Automatic settings with nf-config
   FCFLAGS=${FCFLAGS}' '${NETCDF_FLAGS}' '${NETCDF_LIBS}
elif [ -f $NETCDFHOME/bin/nc-config ] ; then
   # Automatic settings with nc-config
   FCFLAGS=${FCFLAGS}' '${NETCDF_FLAGS}' '${NETCDF_LIBS}
else
   # Manual settings
   FCFLAGS=${FCFLAGS}' -I'${NETCDFINCLUDE}' -L'${NETCDFLIB}' -lnetcdf'
fi

echo "Flags: $FCFLAGS"

#-------------------------------------------------------------------------------
