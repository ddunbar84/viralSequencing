file=$1; 

for line in "$file";

do 
	sed 's/^[^ ]*[|]\([^|]*\)[|] .*$/>\1/'

done

sed 's/^[^ ]*[|]\([^|]*\)[|] .*$/>\1/' fhvGenbank_alignStrip > fhvGenbank_alignStripped 