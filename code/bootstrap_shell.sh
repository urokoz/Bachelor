#!/bin/tcsh -f


foreach a (`seq 0 1 374`)

./spear_xycorr -bt -z ./charts/chart$a.txt | grep -v "#" >> sampled_corr.txt

end
