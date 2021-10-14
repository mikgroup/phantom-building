#!/bin/bash

# put data files in order of increasing TI, without .dat extension
data=( meas_MID43_t1_se_FID281899 meas_MID44_t1_se_FID281900 meas_MID45_t1_se_FID281901 meas_MID46_t1_se_FID281902 meas_MID47_t1_se_FID281903 )

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

bart join 5 ${data[*]} ir_cimg

for dat in "${data[@]}" ; do
	rm $dat.cfl $dat.hdr
done
