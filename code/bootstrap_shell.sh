#!/bin/tcsh -f


foreach a (`seq 0 1 374`)

./spear_xycorr -bt -z ./charts/log_chart$a.txt | grep -v "#" >> Data/log_sampled_corr_SRC.txt

end
