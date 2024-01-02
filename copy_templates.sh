#!/bin/bash
# (Selection of shell)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# copy_templates.sh
#
# Description:
#   Making local, not version-controlled copies of the scripts for executing
#   SICOPOLIS and the run-specification header files.
#
# Author: Ralf Greve
# Date:   2024-01-02
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CP=/bin/cp

RUN_DIR=./runs
TOOLS_DIR=./tools
HEADER_DIR=${RUN_DIR}/headers
TEMPLATE_DIR_NAME=templates

echo "Copying scripts *.sh"
echo "        from ${RUN_DIR}/${TEMPLATE_DIR_NAME}"
echo "        to ${RUN_DIR} ..."
echo " "
$CP -f ${RUN_DIR}/${TEMPLATE_DIR_NAME}/*.sh ${RUN_DIR}

echo "Copying scripts *.sh"
echo "        from ${TOOLS_DIR}/${TEMPLATE_DIR_NAME}"
echo "        to ${TOOLS_DIR} ..."
echo " "
$CP -f ${TOOLS_DIR}/${TEMPLATE_DIR_NAME}/*.sh ${TOOLS_DIR}

echo "Copying run-specification headers *.h"
echo "        from ${HEADER_DIR}/${TEMPLATE_DIR_NAME}"
echo "        to ${HEADER_DIR} ..."
echo " "
$CP -f ${HEADER_DIR}/${TEMPLATE_DIR_NAME}/*.h ${HEADER_DIR}

echo "... done!"
echo " "

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
