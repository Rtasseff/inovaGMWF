
OUTNAME=fullPWOut

# cat small files together in order
cat $(find ./ -name "${OUTNAME}_*" | sort -V) > ${OUTNAME}_tmp.dat
# remove comment lines
sed '/^\#/d' ${OUTNAME}_tmp.dat > ${OUTNAME}.dat 
# (un)comment to keep (remove) extra files, no going back
rm ${OUTNAME}_*

sort -t $'\t' -k 6,6 -g ${OUTNAME}.dat > ${OUTNAME}_sorted.dat
# (un)comment to keep (remove) extra files, no going back
rm ${OUTNAME}.dat


