#!/bin/bash
LANG=C

################################################################################
#
#  t o o l s . s h
#
#  bash script for
#  compilation, linking and execution of the tools.
#
#  Authors: Malte Thoma, Thomas Goelles, Ralf Greve, Fuyuki Saito
#
#  Date: 2020-04-02
#
#    Execute script 
#       ./tools.sh -p <program> -m <run_name> [further options...]
#    where <program> is the name of the program to be executed,
#    and <run_name> is the name of the simulation.
#
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
   echo -e "\n  Usage: `basename $0` -p <program> -m <run_name> [further options...]\n"\
   "     use -p? to get a list of available programs\n"\
   "     use -m? to get a list of available simulations\n"\
   "     [-p <program>]  => name of the program\n"\
   "     [-m <run_name>] => name of the simulation\n"\
   "     [-d <dir>] => individual output directory used for the simulation,\n"\
   "                   default is sico_out/<run_name>\n"\
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
   while getopts bc:d:hm:np:uz? OPT ; do
     case $OPT in
       b) BUILD_ONLY="TRUE";;
       c) CONFIG=$OPTARG ;;
       d) local OUTDIR=$OPTARG ;;
       m) RUN=$OPTARG ;;
       n) SKIP_MAKECLEAN="TRUE";;
       p) PROGNAME=$OPTARG ;;
       u) local REDUNDANT_OPTION_U="TRUE";;
       z) local REDUNDANT_OPTION_Z="TRUE";;
     h|?) usage; exit 1;;
     
     esac            
   done

   # Check if PROGNAME is set correctly
   if [ ! "$PROGNAME" ]; then error "No program chosen. Try option -h."; fi
   if [ "$PROGNAME" == "?" ]; then 
      info "-----------------------"
      info "Available programs are:"
      info "-----------------------"
      # ls *.F90 | sed 's/.F90//g' 
      info "make_ismip_output"
      info "resolution_doubler"
      # info "sicograph"
      exit 1
   fi
   if [ ! -e $PROGNAME/$PROGNAME.F90 ]; then 
      error "Program $PROGNAME.F90 not found. Try options -p? or -h."
   fi 

   # Output directory, absolute paths
   if [ ! "$OUTDIR" ]; then      
      OUTDIR=`pwd`"/../sico_out/"
   else
      lastch=`echo $OUTDIR | sed -e 's/\(^.*\)\(.$\)/\2/'`
      if [ ${lastch} == "/" ]; then OUTDIR=`echo $OUTDIR | sed '$s/.$//'`; fi
      if [ ! -e $OUTDIR ]; then error "$OUTDIR does not exist."; fi       
   fi

   # Check if RUN is set correctly
   if [ ! "$RUN" ]; then error "No simulation set. Try option -h."; fi
   if [ $RUN == "?" ]; then 
      info "Available results in $OUTDIR are:"
      ls $OUTDIR  | sed 's/F_//g' | sed 's/R_.*//g' 
      exit 1
   fi

   RUN_SPECS_HEADER="sico_specs_$RUN.h"

   DATAPATH="$OUTDIR/$RUN" 

   if [ ! -e $DATAPATH ]; then 
      error "$DATAPATH does not exist."
   fi

   # Configuration file
   if [ ! "$CONFIG" ]; then
      CONFIG="../sico_configs.sh"
   else
      CONFIG="../`basename $CONFIG`"
      if [ ! -e $CONFIG ]; then error "$CONFIG does not exist."; fi
   fi

   # Redundant options
   if [ $REDUNDANT_OPTION_U ]; then info "\nOption -u is redundant."; fi
   if [ $REDUNDANT_OPTION_Z ]; then info "\nOption -z is redundant."; fi
}

################################################################################

function compile()
{
   source ../sico_environment.sh

   info "\nConfiguration file ${CONFIG}."
   source $CONFIG

   cd ./${PROGNAME}

   EXE_FILE=${PROGNAME}.x

   $RM -f ${EXE_FILE}
   $RM -f *.mod

   $RM -f $RUN_SPECS_HEADER
   $CP $DATAPATH/$RUN_SPECS_HEADER .

   FCFLAGS="${FCFLAGS} -DRUN_SPECS_HEADER=\"${RUN_SPECS_HEADER}\""

   FCFLAGS="${FCFLAGS} -DOUT_PATH=\"${DATAPATH}\""

   FCFLAGS="${FCFLAGS} -o ${EXE_FILE}"

   ${FC} ./${PROGNAME}.F90 ${FCFLAGS}
}

################################################################################

function run()
{
   if [ $BUILD_ONLY ]; then
      info "Skip execution."
      return 0
   fi

   info "Starting ./${EXE_FILE} ..."
   ./${EXE_FILE}

   if [ ! $SKIP_MAKECLEAN ]; then
      $RM -f ${EXE_FILE}
      $RM -f ${RUN_SPECS_HEADER}
      $RM -f *.mod
   fi

   info "\n... finished."

   cd $OLDPWD
}

################################################################################

RM=/bin/rm
CP=/bin/cp
MV=/bin/mv

check_args $*
compile
run

############################# End of tools.sh ##################################
