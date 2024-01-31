MATLABDIR = /usr/local/MATLAB/R2023b
MATLABEXEC = $(MATLABDIR)/bin/matlab
BAPSFOLDERS = ./admixture ./general ./graph ./independent ./linkage ./parallel ./spatial
MATLABOPTS = -nodisplay -nosplash -nodesktop -batch  # -batch comes last

BAPS_package/run_baps6.sh: $(BAPSFOLDERS)
	@if [ ! -d $(MATLABDIR) ]; then \
		echo "$(MATLABDIR) does not exist"; \
		echo "Please edit the Makefile and set MATLABDIR to a valid path"; \
		exit 1; \
	fi

	@echo -n "Adding BAPS to MATLAB path... "
	@$(MATLABEXEC) $(MATLABOPTS) "run('add_BAPS_to_path.m'); exit;"
	@echo "done"

	@echo "Compiling BAPS"
	@$(MATLABEXEC) $(MATLABOPTS) "run('compileBaps6.m'); exit;"
	@echo "done"

run:
	@echo "Running BAPS"
	bash BAPS_package/run_baps6.sh $(MATLABDIR)
