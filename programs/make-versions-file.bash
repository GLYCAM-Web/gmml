#!/usr/bin/env bash
#
# This can be called from anywhere, but is intended to be called from 
#    the top-level GMML directory.  That is, the intended method
#    of use looks like:
#
#    programs/make-versions-file.bash
# 
echo """
GMML_GIT_BRANCH=\"$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')\"
GMML_GIT_COMMIT_HASH=\"$(git rev-parse HEAD)\"
""" >> VERSIONS.sh

