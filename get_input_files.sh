#!/bin/bash
# (Selection of shell)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get_input_files.sh
#
# Description:
#   Downloading the input files for SICOPOLIS,
#   copying them to the corresponding directories.
#
# Author: Ralf Greve
# Date:   2024-04-15
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#-------- Settings (to be customized) --------

ANT_FLAG=1       # Antarctica:
                 #    1 - get files, 0 - don't get files
GRL_FLAG=1       # Greenland:
                 #    1 - get files, 0 - don't get files
NHEM_FLAG=1      # Entire northern hemisphere:
                 #    1 - get files, 0 - don't get files
SCAND_FLAG=1     # Fennoscandian and Eurasian ice sheets:
                 #    1 - get files, 0 - don't get files
TIBET_FLAG=1     # Tibetan ice sheet:
                 #    1 - get files, 0 - don't get files
ASF_FLAG=1       # Austfonna:
                 #    1 - get files, 0 - don't get files
MOCHO_FLAG=1     # Mocho-Choshuenco ice cap:
                 #    1 - get files, 0 - don't get files
EISMINT_FLAG=1   # EISMINT:
                 #    1 - get files, 0 - don't get files
HEINO_FLAG=1     # ISMIP HEINO:
                 #    1 - get files, 0 - don't get files
NMARS_FLAG=1     # North polar cap of Mars:
                 #    1 - get files, 0 - don't get files
SMARS_FLAG=1     # South polar cap of Mars:
                 #    1 - get files, 0 - don't get files

#-------- Initialization --------

# Zenodo repo:
REPO_URL=https://zenodo.org/record/10976314/files

# Backup repo (no version control, latest version only):
# REPO_URL=https://www2.lowtem.hokudai.ac.jp/gisg/repo/sicopolis/sico_in

SICOPOLIS_HOME=$PWD

domains=

if [[ $ANT_FLAG -eq 1 ]]; then
   domains="`echo $domains` ant"; fi
if [[ $GRL_FLAG -eq 1 ]]; then
   domains="`echo $domains` grl"; fi
if [[ $NHEM_FLAG -eq 1 ]]; then
   domains="`echo $domains` nhem"; fi
if [[ $SCAND_FLAG -eq 1 ]]; then
   domains="`echo $domains` scand"; fi
if [[ $TIBET_FLAG -eq 1 ]]; then
   domains="`echo $domains` tibet"; fi
if [[ $ASF_FLAG -eq 1 ]]; then
   domains="`echo $domains` asf"; fi
if [[ $MOCHO_FLAG -eq 1 ]]; then
   domains="`echo $domains` mocho"; fi
if [[ $EISMINT_FLAG -eq 1 ]]; then
   domains="`echo $domains` eismint"; fi
if [[ $HEINO_FLAG -eq 1 ]]; then
   domains="`echo $domains` heino"; fi
if [[ $NMARS_FLAG -eq 1 ]]; then
   domains="`echo $domains` nmars"; fi
if [[ $SMARS_FLAG -eq 1 ]]; then
   domains="`echo $domains` smars"; fi

TMP_DIR="tmp_`date --iso-8601='seconds'`"
mkdir $TMP_DIR

#-------- Downloading and unpacking files --------

cd $TMP_DIR

echo; echo "Downloading and unpacking files:"

for domain in ${domains}; do
   echo; echo "  ${domain} ..."
   wget ${REPO_URL}/${domain}.tgz
   tar -x -v -z -f ${domain}.tgz
done

#-------- Copying files --------

MV=/bin/mv

SICOPOLIS_INPUT=sico_in

echo; echo "Copying files:"

for domain in ${domains}; do
   echo; echo "  ${domain} ..."
   $MV -f ${domain}/* "${SICOPOLIS_HOME}/${SICOPOLIS_INPUT}/${domain}/"
done

#-------- Clean-up --------

RM=/bin/rm

cd $OLDPWD

$RM -rf $TMP_DIR

echo; echo "... done!"; echo

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
