
OUTNAME=fullPWOut

cat ${OUTNAME}_* > ${OUTNAME}_tmp.dat
sed '/^\#/d' ${OUTNAME}_tmp.dat > ${OUTNAME}.dat 
# (un)comment to keep (remove) extra files, no going back
rm ${OUTNAME}_*

sort -t $'\t' -k 6,6 -g ${OUTNAME}.dat > ${OUTNAME}_sorted.dat
# (un)comment to keep (remove) extra files, no going back
rm ${OUTNAME}.dat


