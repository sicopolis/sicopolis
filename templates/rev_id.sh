#!/bin/bash

################################################################################
#
#  r e v _ i d . s h
#
#  bash script for determining a unique revision identifier from Git.
#
#  Author: Ralf Greve
#
#  Date: 2023-03-25
#
################################################################################

git status >/dev/null 2>&1

if [[ $(echo $?) -eq 0 ]]; then

   BUILD_BRANCH=$(git rev-parse --abbrev-ref HEAD)
   BUILD_REV_ID=$(git rev-parse HEAD)
   BUILD_REV_ID_SHORT=$(git rev-parse --short=9 HEAD)
   COMM_DATE=$(git show -s --date=format:'%Y%m%dT%H%M%z' --format=%cd)

   REPO_REVISION="${BUILD_BRANCH}_${BUILD_REV_ID_SHORT}_${COMM_DATE}"

   if [[ -n $(git status --untracked-files=no --porcelain) ]]; then
      REPO_REVISION="${REPO_REVISION}_dirty"
   fi

else

   REPO_REVISION="no_Git_repository"

fi

export REPO_REVISION

echo "Git revision identifier: ${REPO_REVISION}"

########################### End of rev_id.sh ###################################
