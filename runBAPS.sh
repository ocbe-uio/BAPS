#!/bin/bash
# Runs BAPS on terminal mode
MATLAB_SCRIPTS="run('general/add_BAPS_to_path.m'); run('baps.m');"
matlab -nosplash -nodesktop -r "$MATLAB_SCRIPTS"
