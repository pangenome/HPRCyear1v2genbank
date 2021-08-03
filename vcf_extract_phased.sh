#!/bin/bash

set -eo pipefail

input=$1
ref=$2
threads=$3

#echo "rendering GFA"
gfa=$(basename $input .gz)
# fixme: hack; add "haplotype identifiers" to the reference paths
# so that all paths can be loaded into the GBWT with sample annotations
zcat $input | sed 's/chm13#/chm13#1#/g' | sed 's/grch38#/grch38#1#/g' >$gfa

#echo "building the XG index for $gfa"
#xg=$(basename $input .og).xg
#TEMPDIR=$(pwd) vg convert -t $threads -g -x $gfa >$xg
timer=/home/erikg/.guix-profile/bin/time
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"

echo "building the GBWT index for $gfa"
gbwt=$(basename $input .gfa.gz).gbwt
TEMPDIR=$(pwd) $timer vg gbwt --num-threads $threads --num-jobs $threads -G -o $gbwt --path-regex '(.+?)#(.+?)#(.+?)' --path-fields 'CSH' $gfa

echo "building VCF from $gfa and $gbwt"
vcf=$(basename $input .gfa.gz).vcf
TEMPDIR=$(pwd) $timer vg deconstruct -a -P $ref -g $gbwt -t $threads $gfa >$vcf

rm $gfa
pigz $vcf
pigz $gbwt

