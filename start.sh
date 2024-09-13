# Script to run both versions of BAPS (GUI and terminal)

MATLABDIR=/usr/local/MATLAB/R2023b

# Ask the user if they want to run the GUI version of BAPS
read -r -p "Select which version of BAPS you want to run ([t]erminal [g]ui): " choice

# Run the appropriate script
case "$choice" in
  [tT])
    bash runBAPS.sh
    ;;
  [gG])
    echo "Using MATLAB on $MATLABDIR."
    echo " If this fails, please change the MATLABDIR variable in this script."
    bash BAPS_package/run_baps6.sh $MATLABDIR
    ;;
  *)
    echo "Invalid choice. Please select either 't' or 'g'."
    ;;
esac
