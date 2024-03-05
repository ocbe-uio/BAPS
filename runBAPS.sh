#!/bin/bash
# Runs BAPS on terminal mode
MATLAB_SCRIPTS="run('add_BAPS_to_path.m'); run('general/baps6_terminal.m');"
matlab -nosplash -nodesktop -r "$MATLAB_SCRIPTS"
