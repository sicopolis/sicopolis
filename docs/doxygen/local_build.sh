#!/bin/bash
# (Selection of shell)

#--------------------------------------------------------------------------
# Building the SICOPOLIS developer manual with Doxygen locally.
#--------------------------------------------------------------------------

# Requirements:
# -------------
# Ensure that you have Doxygen version 1.10.0 or later installed.

# Doxygen config file (in docs/doxygen/doxygen-config):
# -----------------------------------------------------
# Copy the config template:
# cp Doxyfile_Template.txt my_Doxyfile.txt
#
# Edit my_Doxyfile.txt:
#  * Search for "Revision ...", and replace "..." with the current
#    revision number
#    (execute `./rev_id.sh` in the main directory to find out).
#  * Search for "/home/username/Documents/sicopolis/src", and replace
#    this by the actual path of the src directory.

# Building the HTML manual:
# -------------------------

echo " "

RM=/bin/rm
MV=/bin/mv

if [ -d html ]
then
   $RM -rf html/*
fi

if [ -d latex ]
then
   $RM -rf latex/*
fi

cd ../../

xyz_dir=src/subroutines/xyz
if ls "$xyz_dir"/*90 > /dev/null 2>&1
then
   tmp_dir="tmp_`date --iso-8601='seconds'`"
   mkdir "$tmp_dir"
   $MV "$xyz_dir"/*90 "$tmp_dir"/
fi

cd docs/doxygen/doxygen-config/

doxygen my_Doxyfile.txt

cd ../../../

if [ -d "$tmp_dir" ]
then
   $MV "$tmp_dir"/*90 "$xyz_dir"/
   $RM -rf "$tmp_dir"
fi

cd docs/doxygen/

# Open html/index.html to see the locally built manual.

# Building the PDF manual (optional, requires LaTeX):
# ---------------------------------------------------

# cd latex
# make
# cd ../

# Open latex/refman.pdf to see the locally built manual.

echo " "
echo "==> Done."
echo " "

#--------------------------------------------------------------------------
#
