# while read line
# do
#     echo $line
# done < inactiveBarcodes.txt

# grep -v '^#' inactiveBarcodes.txt | while read line ; do
#     echo $line
# done

barcodes=$1

# for ln in "$barcodes":
#     if  ln.startswith('#'):
#         lnInactive = ln[1:]
#         print (lnInactive)

while read old new; do
  for fastq in ./*/"$old"*.fastq; do
    new_name=$new${fastq##*/"$old"}
    echo mv "$fastq" "${fastq%/*}/$new_name"
  done
done <name-change.txt