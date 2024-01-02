#!/bin/bash
# (Selection of shell)

#--------------------------------------------------------------------------
# Building the SICOPOLIS manual with Sphinx locally.
#--------------------------------------------------------------------------

# Requirements:
# -------------
# Ensure that you have pip (package installer for Python)
# and Pandoc (document converter) installed.

# cd docs/sphinx
# pip install -r requirements.txt   # \!/ not as root

# Actual build:
# -------------

make clean
make html
# make latexpdf PAPER=a4

# Open docs/sphinx/build/html/index.html to see the locally built manual.

#--------------------------------------------------------------------------
#
