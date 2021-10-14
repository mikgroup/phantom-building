#!/bin/bash

# without .dat extension
data=( meas_MID48_t1_ir_se_TI100_FID287556 meas_MID49_t1_ir_se_TI300_FID287557 meas_MID50_t1_ir_se_TI500_FID287558 meas_MID51_t1_ir_se_TI700_FID287559 meas_MID52_t1_ir_se_TI900_FID287560  )

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
