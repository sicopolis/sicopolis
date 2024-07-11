#!/bin/bash

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# make_climatology.sh   Date: 2024-07-11
#
# -----------
# Description
# -----------
#
#   Generate climatology (SMB, ST, reference elevation) for a given
#   time interval from a SICOPOLIS output file.
#
#   SICOPOLIS must have been run with the settings
#      #define OUTPUT_FLUX_VARS 2
#      #define CLIMATOLOGY_EXTRACTION_HACK
#   in the run-specs header, and the output file must be valid for the
#   desired time interval (this is not checked!).
#
# --------
# Examples
# --------
#
#   ./make_climatology.sh \
#      -d /work/sicopolis/ant08_spinup_clim \
#      -m ant08_spinup_clim -n 0003 -y 1960 -z 1989
#
#     -> create climatology (ant08_spinup_clim_1960-1989.nc)
#           from output file ant08_spinup_clim0003.nc
#              in directory /work/sicopolis/ant08_spinup_clim
#                 for the period 1960-1989.
#
#   ---------
#
#   ./make_climatology.sh \
#      -d /work/sicopolis/ant08_hist_clim \
#      -m ant08_hist_clim -n 0002 -t -y 1995 -z 2014
#
#     -> create climatology (ant08_hist_clim_1995-2014.nc)
#           from output file ant08_hist_clim_2d_0002.nc
#              in directory /work/sicopolis/ant08_hist_clim
#                 for the period 1995-2014.
#
#   ---------
#
#   ./make_climatology.sh \
#      -d /work/sicopolis/grl05_paleo_clim \
#      -m grl05_paleo_clim -n 0004 -y 1960 -z 1989
#
#     -> create climatology (grl05_paleo_clim_1960-1989.nc)
#           from output file grl05_paleo_clim0004.nc
#              in directory /work/sicopolis/grl05_paleo_clim
#                 for the period 1960-1989.
#
# ----
# Note
# ----
#
#   The resulting climatology file must be moved to
#   sico_in/ant/ (for Antarctica) or sico_in/grl/ (for Greenland).
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# set -x

#-------- Flags --------

while getopts d:m:n:ty:z: flag
do
   case ${flag} in
      d) outdir=${OPTARG};;
      m) runname=${OPTARG};;
      n) ergnum=${OPTARG};;
      t) two_d_file="TRUE";;
      y) year1=${OPTARG};;
      z) year2=${OPTARG};;
   esac            
done

if [ $two_d_file ]; then
   two_d="_2d_"
else
   two_d=""
fi

#-------- Settings --------

dir1="${outdir}"
dir2="."

filename1="${runname}${two_d}${ergnum}.nc"
filename2="${runname}_${year1}-${year2}.nc"
filename3="${dir1}/${filename1}"
filename4="${dir2}/${filename2}"

#-------- Extracting climatology from SICOPOLIS output file --------

echo "Creating ${filename2}..." ;

# ncks -O -F -v crs,x,y,lon,lat,temp_maat,as_perp,zs ${filename3} ${filename4}
ncks -O -F -v mapping,x,y,lon,lat,temp_maat,as_perp,zs ${filename3} ${filename4}
ncrename -v temp_maat,ST_clim ${filename4}
ncrename -v as_perp,SMB_clim ${filename4}
ncrename -v zs,surf_elev_ref ${filename4}

#-------- End of script --------

echo "Done." ;

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
