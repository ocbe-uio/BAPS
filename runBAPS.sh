#!/bin/bash
# Runs BAPS on terminal mode
MATLAB_STARTUP_SCRIPT="run('general/add_BAPS_to_path.m'); run('general/welcome.m');"
MATLAB_PATH="/usr/local/MATLAB/R2023b/bin/matlab"
"$MATLAB_PATH" -nosplash -nodesktop -r "$MATLAB_STARTUP_SCRIPT"
