#!/bin/bash

# Input:
# $1: file containing a list of taxids
# e.g. "/home/kcw2/data/testing/kraken2_test/dataset/ncbi_dataset/data/taxids_pseudomonas.txt"

# Output: json patterns file for grepq (same name as input file but with .json extension)

name="${1%.*}"
outname="${name}.json"

patternsPlain=$(mktemp)
for taxid in $(cat $1); do
  echo "taxid\|${taxid}\\b" >> $patternsPlain # enforce word boundary so that taxid|123 doesn't match taxid|12345
done

# make yaml (easier than directly making json)
yamlTemp=$(mktemp)
echo "regexSet:" > $yamlTemp
echo "  regexSetName: taxids in headers" >> $yamlTemp
echo "  regex:" >> $yamlTemp
echo "  - regexName: dummy" >> $yamlTemp
echo "    regexString: '[A-Za-z]+'" >> $yamlTemp

patterns=$(paste -sd "|" $patternsPlain)
echo "  headerRegex: '$patterns'" >> $yamlTemp

echo "  minimumSequenceLength: 0" >> $yamlTemp
echo "  minimumAverageQuality: 0" >> $yamlTemp
echo "  qualityEncoding: Phred+33" >> $yamlTemp
rm -f $patternsPlain

# now convert yaml to json
yq . $yamlTemp > $outname
rm -f $yamlTemp