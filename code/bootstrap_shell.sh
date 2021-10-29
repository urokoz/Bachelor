#!/bin/tcsh -f


foreach a (`seq 0 1 98`)

./spear_xycorr -bt -z ./Data/birch/charts/log_chart$a.txt | grep -v "#" >> Data/birch/significances/birch_log_SCC_sig.txt
./xycorr -bt -z ./Data/birch/charts/log_chart$a.txt | grep -v "#" >> Data/birch/significances/birch_log_PCC_sig.txt
./spear_xycorr -bt -z ./Data/birch/charts/chart$a.txt | grep -v "#" >> Data/birch/significances/birch_SCC_sig.txt
./xycorr -bt -z ./Data/birch/charts/chart$a.txt | grep -v "#" >> Data/birch/significances/birch_PCC_sig.txt

end
