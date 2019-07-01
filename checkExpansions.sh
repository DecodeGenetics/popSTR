#!/bin/bash
BAM=$1
MaxRepLength=$2
for i in {2,3,4,5,6,7,9,12,13,14,16,18,19,22}
do
    ./computeReadAttributes $BAM . ./panelMarkerInfo/chr${i}markerInfo 4 $MaxRepLength chr${i}
done

while read -r fileName limit
do
    lineAboveLimit=`awk -v limit="$limit" 'NF==8 && substr($1,1,3) != "chr" && $1>=limit' ./attributes/${fileName} | wc -l`
    motifLength=`echo $fileName | cut -d '_' -f 2 | awk '{print length($1)}'`
    lineAtLengthLimit=`awk -v motifLength="$motifLength" -v maxRepLength="$MaxRepLength" 'NF==8 && substr($1,1,3) != "chr" && $1>=maxRepLength/motifLength' ./attributes/${fileName} | wc -l`
    echo "Marker $fileName has $lineAboveLimit reads supporting a pathogenic expansion and $lineAtLengthLimit reads with repeat longer than $MaxRepLength bases."
done < ./markersAndLimits
