#!/bin/bash
#
# SUBMODULE_ADD
# A git routine for adding mat-tfer-pma as a submodule into a program.
# Adds submodule to tfer-pma folder.
# Author:  Timothy Sipkens, 2021-01-31
#===========================================================#

git submodule add -b master https://github.com/tsipkens/mat-tfer-pma tfer-pma
git submodule init
