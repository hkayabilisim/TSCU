#!/bin/sh

ls ../../UCR/*/*_TRAIN | grep -v AVIRIS | grep TwoLeadECG | while read trnfile
do
	echo "********* $trnfile"
	cp $trnfile TRNFILE
	cat createfeatures.R | R --no-save	
	mv TRNFILE.features $trnfile.features
done
