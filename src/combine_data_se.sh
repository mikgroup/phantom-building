#!/bin/bash

# without .dat extension
data=( meas_MID36_se_FID281892  meas_MID41_se_FID281897 meas_MID37_se_FID281893 meas_MID39_se_FID281895 meas_MID40_se_FID281896 )

for dat in "${data[@]}" ; do
	bart twixread -A $dat.dat $dat-ksp.ra
	bart transpose 2 13 $dat-ksp.ra $dat-ksp-trp.ra
	rm $dat-ksp.ra
	#for i in `seq 0 24` ; do
		#bart slice 2 $i $dat-ksp-trp $dat-ksp-trp-$i
		#bart ecalib -m 1 -k 5 -c .1 $dat-ksp-trp-$i $dat-sens-$i.ra
	#done
	#bart join 3 $dat-sens-* $dat-sens
	bart fft -iu 3 $dat-ksp-trp.ra $dat
	rm $dat-ksp-trp.ra
	#bart fmac -C -s 8 $dat-cimg $dat-sens $dat
done

bart join 5 ${data[*]} se_cimg

for dat in "${data[@]}" ; do
	rm $dat.cfl $dat.hdr
done
