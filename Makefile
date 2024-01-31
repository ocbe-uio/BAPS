MATLABDIR = /usr/local/MATLAB/R2023b
MATLABEXEC = $(MATLABDIR)/bin/matlab
BAPSFOLDERS = ./admixture ./general ./graph ./independent ./linkage ./parallel ./spatial

BAPS_package/run_baps6.sh: $(BAPSFOLDERS)
	@if [ ! -d $(MATLABDIR) ]; then \
		echo "$(MATLABDIR) does not exist"; \
		echo "Please edit the Makefile and set MATLABDIR to a valid path"; \
		exit 1; \
	fi

	@echo -n "Adding BAPS to MATLAB path... "
	@$(MATLABEXEC) -nodisplay -nosplash -nodesktop -batch "run('add_BAPS_to_path.m'); exit;"
	@echo "done"

	@echo -n "Compiling BAPS... "
	@$(MATLABEXEC) -nodisplay -nosplash -nodesktop -batch "run('compileBaps6.m'); exit;"
	@echo "done"

run:
	@echo "Running BAPS"
	bash BAPS_package/run_baps6.sh $(MATLABDIR)
