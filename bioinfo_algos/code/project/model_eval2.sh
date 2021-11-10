#! /bin/tcsh -f

## Define path to your code directory
set RDIR = /home/mathias/bioinfo_algos/code/

## Define path you where you have placed the HLA data sets
set DDIR = /home/mathias/bioinfo_algos/data/project

set prog = True

# Here you can type your allele names A0301 A6802 B4403
foreach a ( $argv )

# create allele folder
mkdir -p $a.res

cd $a.res

# # create PSSM folder
# mkdir -p PSSM
#
# cd PSSM
#
# # Sequence weighting
# foreach sw ( 0 1 )
#
# mkdir -p sw_$sw
#
# cd sw_$sw
#
# # Vary beta values
# foreach beta ( 0 50 100 200)
#
# mkdir -p beta_$beta
#
# cd beta_$beta
#
# foreach n ( 0 1 2 3 4 )
#
# if ( $prog == True ) echo $a PSSM beta: $beta SW: $sw cv: $n
#
# # Do training
# if ( $sw == 0 ) then
# 	if ( ! -e mat.$n ) then
# 		python $RDIR/PSSM/pep2mat.py -b $beta -f $DDIR/$a/pt00$n | grep -v "#" > mat.$n
# 	endif
# else
# 	if ( ! -e mat.$n ) then
# 		python $RDIR/PSSM/pep2mat.py -b $beta -w -f $DDIR/$a/pt00$n | grep -v "#" > mat.$n
# 	endif
# endif
# # Do evaluation
# if ( ! -e c00$n.pred ) then
# 	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
# endif
#
# end
#
# echo $a PSSM beta: $beta SW: $sw `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/Gibbs/xycorr` \
# 	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' ` >> $RDIR/project/$a\_models.data
#
# cd .. 	# beta -> sw
#
# end  		# beta
#
# cd .. 	# sw -> PSSM
#
# end  		# sw
#
#
# cd ..		# PSSM -> a

# create SMM folder
mkdir -p SMM

cd SMM

# gadient decent
mkdir -p GD

cd GD

# Here you can type the lambdas to test
foreach l ( 5 10 50 )

mkdir -p l.$l

cd l.$l

# Loop over the 5 cross validation configurations
foreach n ( 0 1 2 3 4 )

if ( $prog == True ) echo $a SMM GD lambda: $l cv: $n

# Do training
if ( ! -e mat.$n ) then
	python $RDIR/SMM/smm_gradient_descent.py -l $l -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n | grep -v "#" > mat.$n
endif

# Do evaluation
if ( ! -e c00$n.pred ) then
	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
endif

end		# cv

# Do concatinated evaluation
echo $a SMM GD lambda: $l `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/Gibbs/xycorr` \
	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' ` >> $RDIR/project/$a\_models.data

cd ..		# lambda -> GD

end			# lambda loop

cd ..		# GD -> SMM

mkdir -p MC

cd MC

# Here you can type the lambdas to test
foreach l (5 10 50)

mkdir -p l.$l

cd l.$l

# Loop over the 5 cross validation configurations
foreach n ( 0 1 2 3 4 )

if ( $prog == True ) echo $a SMM MC lambda: $l cv: $n

# Do training
if ( ! -e mat.$n ) then
	python $RDIR/SMM/smm_monte_carlo.py -l $l -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n | grep -v "#" > mat.$n
endif

# Do evaluation
if ( ! -e c00$n.pred ) then
	python $RDIR/PSSM/pep2score.py -mat mat.$n -f  $DDIR/$a/c00$n | grep -v "PCC:" > c00$n.pred
endif

end		# cv

# Do concatinated evaluation
echo $a SMM MC lambda: $l `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/Gibbs/xycorr` \
	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' ` >> $RDIR/project/$a\_models.data

cd ..		# lambda -> MC

end			# lambda

cd ..		# MC -> SMM

cd ..		# SMM -> a
#
# # make ANN folder
# mkdir -p ANN
#
# cd ANN
#
# foreach eta (0.1 0.5 1)
#
# # make dir for each eta
# mkdir -p eta.$eta
#
# cd eta.$eta
#
# foreach nh (2 3 5)
#
# # make dir for each hidden neuron configuration
# mkdir -p nh.$nh
#
# cd nh.$nh
#
# # make dir for blosum encoding
# mkdir -p bl
#
# cd bl
#
# # Loop over the 5 cross validation configurations
# foreach n ( 0 1 2 3 4 )
#
# if ( $prog == True ) echo $a ANN eta: $eta nh: $nh bl cv: $n
#
# if ( ! -e c00$n.out) then
# 	python $RDIR/ANN/ANN_train.py -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n -epi $eta  -syn syn.$n -bl -nh $nh -stop | grep -v "#" > c00$n.out
# endif
#
# if ( ! -e c00$n.pred ) then
# 	python $RDIR/ANN/ANN_forward.py -e $DDIR/$a/c00$n -syn syn.$n -bl | grep -v "#" > c00$n.pred
# endif
#
# end			# cv
#
# echo $a ANN eta: $eta nh: $nh bl `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/Gibbs/xycorr` \
# 	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' ` >> $RDIR/project/$a\_models.data
#
#
# cd ..		# bl -> nh
#
# # make dir for sparse encoding
# mkdir -p sp
#
# cd sp
#
# # Loop over the 5 cross validation configurations
# foreach n ( 0 1 2 3 4 )
#
# if ( $prog == True ) echo $a ANN eta: $eta nh: $nh sp cv: $n
#
# if ( ! -e c00$n.out) then
# 	python $RDIR/ANN/ANN_train.py -t $DDIR/$a/f00$n -e $DDIR/$a/c00$n -epi $eta  -syn syn.$n -nh $nh -stop | grep -v "#" > c00$n.out
# endif
#
# if ( ! -e c00$n.pred ) then
# 	python $RDIR/ANN/ANN_forward.py -e $DDIR/$a/c00$n -syn syn.$n | grep -v "#" > c00$n.pred
# endif
#
# end			# cv
#
# echo $a ANN eta: $eta nh: $nh sp `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | $RDIR/Gibbs/xycorr` \
# 	   			 `cat c00?.pred | grep -v "#" | gawk '{print $2,$3}' | gawk 'BEGIN{n+0; e=0.0}{n++; e += ($1-$2)*($1-$2)}END{print e/n}' ` >> $RDIR/project/$a\_models.data
#
# cd .. 	# sp -> nh
#
# cd ..		# nh -> eta
#
# end			# nh
#
# cd ..		# eta -> ANN
#
# end			# eta
#
# cd ..		# ANN -> a

cd ..		# a -> /

end			# a
