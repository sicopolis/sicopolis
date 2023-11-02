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
   (./sico.sh ${MULTI_OPTIONS_1} -m repo_vialov3d25) \
              >out_multi_101.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_emtp2sge25_expA) \
              >out_multi_102.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl16_bm5_ss25ka) \
              >out_multi_103.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl16_bm5_init100a) \
              >out_multi_104.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl16_bm5_ss25ka_nudged \
              -t ${MULTI_OUTDIR}/repo_grl16_bm5_init100a) \
              >out_multi_105.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_ant40_b2_ss25ka) \
              >out_multi_106.dat 2>&1

   #--------

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl20_b2_paleo21) \
              >out_multi_111.dat 2>&1

   cd $PWD/../tools ; echo 0004 | \
   (./tools.sh -p resolution_doubler ${MULTI_OPTIONS_2} \
               -m repo_grl20_b2_paleo21) \
               >$OLDPWD/out_multi_112.dat 2>&1 ; cd $OLDPWD

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl10_b2_paleo21 \
              -a ${MULTI_OUTDIR}/repo_grl20_b2_paleo21) \
              >out_multi_113.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl10_b2_future21_ctrl \
              -a ${MULTI_OUTDIR}/repo_grl10_b2_paleo21) \
              >out_multi_114.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_grl10_b2_future21_asmb \
              -a ${MULTI_OUTDIR}/repo_grl10_b2_paleo21) \
              >out_multi_115.dat 2>&1

   #--------

   # NOTE: The following simulations require Lis and OpenMP.
   #       In order to run them, make sure that both LIS_FLAG and OPENMP_FLAG
   #       are set to "true" in sico_configs.sh.

   OMP_NUM_THREADS=1; export OMP_NUM_THREADS
   #              (number of threads for the SSA solver using OpenMP)

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_ant64_b2_spinup09_init100a) \
              >out_multi_121.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_ant64_b2_spinup09_fixtopo \
              -a ${MULTI_OUTDIR}/repo_ant64_b2_spinup09_init100a \
              -t ${MULTI_OUTDIR}/repo_ant64_b2_spinup09_init100a) \
              >out_multi_122.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_ant64_b2_spinup09 \
              -a ${MULTI_OUTDIR}/repo_ant64_b2_spinup09_fixtopo) \
              >out_multi_123.dat 2>&1

   (./sico.sh ${MULTI_OPTIONS_1} -m repo_ant64_b2_future09_ctrl \
              -a ${MULTI_OUTDIR}/repo_ant64_b2_spinup09) \
              >out_multi_124.dat 2>&1
}

################################################################################

check_args $*
run

######################## End of multi_sico_1.sh ################################
