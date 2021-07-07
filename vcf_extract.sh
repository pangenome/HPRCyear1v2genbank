#!/bin/bash

input=$1
ref=$2
threads=$3

#echo "rendering GFA"
gfa=$(basename $input .gz)
zcat $input >$gfa
#odgi view -i $input -g >$gfa

#echo "building the XG index for $gfa"
#xg=$(basename $input .og).xg
#TEMPDIR=$(pwd) vg convert -t $threads -g -x $gfa >$xg
timer=/home/erikg/.guix-profile/bin/time
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"

echo "building VCF from $gfa"
vcf=$(basename $input .gfa.gz).vcf
TEMPDIR=$(pwd) $timer -f "$fmt" vg deconstruct -P $ref \
   $(for i in $(awk '$1 == "P" { print $2 }' $gfa | grep -v Cons | cut -f 1 -d '#' | sort | uniq ); do echo -n ' -A '$i; done) \
   -e -a -t $threads $gfa >$vcf

pigz $vcf
rm $gfa

