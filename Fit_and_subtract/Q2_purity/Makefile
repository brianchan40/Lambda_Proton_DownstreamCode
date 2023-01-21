CENTRALITIES = 0 1 2 3 4 5 6 7 8
LAMBDA = lam antilam
#CENTRALITIES = 0

submit_analysis:
	for lam in $(LAMBDA) ; do for cen in $(CENTRALITIES) ; do ./run.sh $$cen $$lam; done ; done ;