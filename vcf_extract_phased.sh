#!/bin/bash

set -eo pipefail

input=$1
ref_spec=$2
threads=$3

prefix=$(basename $input .gfa.gz)

gfa=$(basename $input .gz)
# fixme: hack; add "haplotype identifiers" to the reference paths
# so that all paths can be loaded into the GBWT with sample annotations
zcat $input | sed 's/chm13#/chm13#1#/g' | sed 's/grch38#/grch38#1#/g' >$gfa

timer=/home/erikg/.guix-profile/bin/time
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"

echo "building the GBWTGraph index for $gfa"
gbwt=$prefix.gbwt
gg=$prefix.gg
TEMPDIR=$(pwd) $timer vg gbwt --num-threads $threads --num-jobs $threads -G -g $gg -o $gbwt --path-regex '(.+?)#(.+?)#(.+?)' --path-fields 'CSH' $gfa

echo "building the XG index for $gfa"
xg=$prefix.xg
TEMPDIR=$(pwd) $timer vg convert $gg -b $gbwt -t $threads -x >$xg

for ref in $( echo "$ref_spec" | tr ',' ' ' );
do
    vcf=$prefix.$ref.vcf
    echo "building VCF $vcf from $gfa and $gbwt"
    TEMPDIR=$(pwd) $timer vg deconstruct -a -P $ref -g $gbwt -t $threads $xg >$vcf
    pigz -p $threads $vcf
done

rm $gfa
pigz -p $threads $gbwt
pigz -p $threads $gg
pigz -p $threads $xg
