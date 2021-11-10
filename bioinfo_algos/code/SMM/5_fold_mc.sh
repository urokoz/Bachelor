#! /bin/tcsh -f

## Define path to your code directory
set RDIR = /home/mathias/bioinfo_algos/code/

## Define path you where you have placed the HLA data sets
set DDIR = /home/mathias/bioinfo_algos/data/SMM

# Here you can type your allele names A0301 A6802 B4403
foreach a ( B4403 )

mkdir -p $a.res

cd $a.res

mkdir -p MC

cd MC

# Here you can type the lambdas to test
foreach l ( 0 0.02 )

mkdir -p l.$l

cd l.$l

# Loop over the 5 cross validation configurations
foreach n ( 0 1 2 3 4 )

# Do training
if ( ! -e mat.$n ) then
	python $RDIR/SMM/smm_monte_carlo.py -l $l -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n | grep -v "#" > mat.$n
endif

# Do evaluation
if ( ! -e c00$n.pred ) then
	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
endif

end

# Do concatinated evaluation
echo $a $l `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | /home/mathias/bioinfo_algos/code/Gibbs/xycorr` \
	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' `

cd ..

end

cd ..

cd ..

end
