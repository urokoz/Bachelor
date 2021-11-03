#!/bin/tcsh -f


foreach a (`seq 0 1 375`)

./spear_xycorr -bt -z ./Data/ragweed/charts/log_chart$a.txt | grep -v "#" >> Data/ragweed/significances/ragweed_log_SCC_sig.txt
./xycorr -bt -z ./Data/ragweed/charts/log_chart$a.txt | grep -v "#" >> Data/ragweed/significances/ragweed_log_PCC_sig.txt
./spear_xycorr -bt -z ./Data/ragweed/charts/chart$a.txt | grep -v "#" >> Data/ragweed/significances/ragweed_SCC_sig.txt
./xycorr -bt -z ./Data/ragweed/charts/chart$a.txt | grep -v "#" >> Data/ragweed/significances/ragweed_PCC_sig.txt

end
