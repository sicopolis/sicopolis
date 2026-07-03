#!/bin/bash

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# make_c_slide_dimless.sh
#
# Description
# -----------
# Generate c_slide_dimless (2D array of dimensionless sliding coefficient)
# from a SICOPOLIS output file.
#
# Example
# -------
#
# ./make_c_slide_dimless.sh \
#    -d /work/sicopolis/grl04_spinup \
#    -m grl04_spinup -n 0002 -s 10
#
#   -> create c_slide_dimless file ('grl04_spinup0002_c_slide_dimless.nc')
#         from output file 'grl04_spinup0002.nc'
#            in directory '/work/sicopolis/grl04_spinup'
#               using a scale of 10 [m a-1 Pa-(p-q)]
#               for the non-dimensionalization.
#
# Note
# ----
# The resulting c_slide_dimless file must be moved to
# sico_in/ant/ (for Antarctica) or sico_in/grl/ (for Greenland).

# Author: Ralf Greve
# Date:   2026-07-03
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------- Flags --------

while getopts d:m:n:s: flag
do
   case ${flag} in
      d) outdir=${OPTARG};;
      m) runname=${OPTARG};;
      n) ergnum=${OPTARG};;
      s) c_slide_scale=${OPTARG};;
   esac            
done

#-------- Settings --------

dir1="${outdir}"
dir2="."

filename1="${runname}${ergnum}.nc"
filename2="${runname}${ergnum}_c_slide_dimless.nc"
filename3="${dir1}/${filename1}"
filename4="${dir2}/${filename2}"

#-------- Computing c_slide_dimless from SICOPOLIS output file --------

ncks -O -F -v mapping,x,y,c_slide ${filename3} c_slide_tmp1.nc
ncap2 -O -F -s "c_slide_dimless=c_slide/${c_slide_scale}f" c_slide_tmp1.nc c_slide_tmp2.nc
ncks -O -F -x -v c_slide c_slide_tmp2.nc ${filename4}
ncatted -O -a units,c_slide_dimless,d,, ${filename4}
ncatted -O -a source,global,d,, ${filename4}
ncatted -O -a title,global,d,, ${filename4}
ncatted -O -a references,global,d,, ${filename4}
ncatted -O -a institution,global,d,, ${filename4}

RM=/bin/rm

$RM -f c_slide_tmp*.nc

#-------- End of script --------

echo "Done." ;

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
