#!/bin/bash

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# make_smb_corr.sh   Date: 2024-07-12
#
# -----------
# Description
# -----------
#
#   Extract surface mass balance (SMB) correction for a given
#   time interval from a SICOPOLIS output file.
#
#   SICOPOLIS must have been run with the setting
#
#      #define OUTPUT_FLUX_VARS 2
#
#   in the run-specs header, and the output file must be valid for the
#   desired time interval (this is not checked!).
#
# -------
# Example
# -------
#
#   ./make_smb_corr.sh \
#      -d /work/sicopolis/grl05_spinup \
#      -m grl05_spinup -n 0001 -y 1960 -z 1989
#
#     -> extract SMB correction (grl05_spinup0001_smb_corr_1960-1989.nc)
#           from output file grl05_spinup0001.nc
#              in directory /work/sicopolis/grl05_spinup
#                 for the period 1960-1989.
#
# ----
# Note
# ----
#
#   The resulting SMB correction file must be moved to
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
filename2="${runname}${two_d}${ergnum}_smb_corr_${year1}-${year2}.nc"
filename3="${dir1}/${filename1}"
filename4="${dir2}/${filename2}"

#-------- Extracting SMB correction from SICOPOLIS output file --------

echo "Creating ${filename2}..."

# ncks -O -F -v crs,x,y,lon,lat,smb_corr ${filename3} ${filename4}
ncks -O -F -v mapping,x,y,lon,lat,smb_corr ${filename3} ${filename4}
ncrename -v smb_corr,DSMB ${filename4}

#-------- End of script --------

echo "Done."

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
