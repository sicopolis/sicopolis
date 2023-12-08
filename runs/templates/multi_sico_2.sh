#!/bin/bash
LANG=C

################################################################################

function error()
{
   echo -e "\n******************************************************************************">&2
   echo -e   "** ERROR:                                                                   **">&2
   echo -e   "$1" | awk '{printf("** %-73s**\n",$0)}'>&2
   echo -e "******************************************************************************\n">&2
   exit 1
}

################################################################################

function info()
{
   echo -e "$1" >&2
}

################################################################################

function usage()
{
   echo -e "\n  Usage: `basename $0` [options...]\n"\
   "     [-i <dir>] => individual input directory, default is sico_in\n"\
   "     [-d <dir>] => individual output directory, default is sico_out/<run_name>\n"\
   "     [-f] => force overwriting the output directory\n"\
   "     [-n] => skip make clean\n"\
   "     [-b] => skip execution, build only\n"\
   "     [-c <FILE>] => configuration file FILE instead of sico_configs.sh\n"\
   "     [-u] => (redundant option)\n"\
   "     [-z] => (redundant option)\n"\

      if [ "$1" ]; then error "$1"; fi
}

################################################################################

function check_args()
{
   while getopts bc:d:fhi:nuz? OPT ; do
     case $OPT in
       b) local BUILD_ONLY="TRUE";;
       c) local CONFIG=$OPTARG ;;
       d) local MULTI_OUTDIR_ARG=$OPTARG ;;
       f) local FORCE="TRUE";;
       i) local INDIR=$OPTARG ;;
       n) local SKIP_MAKECLEAN="TRUE";;
       u) local REDUNDANT_OPTION_U="TRUE";;
       z) local REDUNDANT_OPTION_Z="TRUE";;
     h|?) usage; exit 1;;
     esac            
   done

   MULTI_OPTIONS_1=' '
   MULTI_OPTIONS_2=' '

   if [ ! "$MULTI_OUTDIR_ARG" ]; then
      MULTI_OUTDIR=${PWD}"/../sico_out"
   else
      lastch=`echo $MULTI_OUTDIR_ARG | sed -e 's/\(^.*\)\(.$\)/\2/'`
      if [ ${lastch} == "/" ]; then
         MULTI_OUTDIR_ARG=`echo $MULTI_OUTDIR_ARG  | sed '$s/.$//'`
      fi
      MULTI_OUTDIR=${MULTI_OUTDIR_ARG}
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -d $MULTI_OUTDIR"
      MULTI_OPTIONS_2="$MULTI_OPTIONS_2 -d $MULTI_OUTDIR"
   fi

   if [ "$FORCE" ]; then
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -f"
   fi

   if [ "$INDIR" ]; then
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -i $INDIR"
   fi

   if [ "$SKIP_MAKECLEAN" ]; then
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -n"
      MULTI_OPTIONS_2="$MULTI_OPTIONS_2 -n"
   fi

   if [ "$BUILD_ONLY" ]; then
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -b"
      MULTI_OPTIONS_2="$MULTI_OPTIONS_2 -b"
   fi

   if [ "$CONFIG" ]; then
      MULTI_OPTIONS_1="$MULTI_OPTIONS_1 -c $CONFIG"
      MULTI_OPTIONS_2="$MULTI_OPTIONS_2 -c $CONFIG"
   fi

   info "Options for  sico.sh:$MULTI_OPTIONS_1"
   info "Options for tools.sh:$MULTI_OPTIONS_2"

   # Redundant options
   if [ $REDUNDANT_OPTION_U ]; then info "\nOption -u is redundant."; fi
   if [ $REDUNDANT_OPTION_Z ]; then info "\nOption -z is redundant."; fi
}

################################################################################

function run()
{
   (./sico.sh ${MULTI_OPTIONS_1} -m repo_asf2_steady) \
              >out_multi_201.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_asf2_surge) \
              >out_multi_202.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_scand_test) \
              >out_multi_203.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_tibet_test) \
              >out_multi_204.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_nmars10_steady) \
              >out_multi_205.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_smars10_steady) \
              >out_multi_206.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_nhem80_nt012_new) \
              >out_multi_207.dat 2>&1

   #--------

   ## OMP_NUM_THREADS=1; export OMP_NUM_THREADS
   ## #              (number of threads for the SSA solver using OpenMP)
   ## 
   ## (./sico.sh ${MULTI_OPTIONS_1} -m repo_emtshelf25_expH) \
   ##            >out_multi_211.dat 2>&1

   #--------

   # !!! WARNING: Uncommenting the following will overwrite any self-written
   #              routines in sicopolis/src/subroutines/xyz !!!

   ## cd $PWD/../src/subroutines/xyz ; $CP -f ./heino/*90 ./ ; cd $OLDPWD
   ## 
   ## (./sico.sh ${MULTI_OPTIONS_1} -m repo_heino50_st) \
   ##            >out_multi_212.dat 2>&1
}

################################################################################

RM=/bin/rm
CP=/bin/cp
MV=/bin/mv

check_args $*
run

######################## End of multi_sico_2.sh ################################
