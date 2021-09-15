#!/bin/bash

input=$1
vcf_spec=$2
threads=$3

gfa=$(basename $input .gz)
zcat $input >$gfa

timer=/home/erikg/.guix-profile/bin/time
fmt="%C\n%Us user %Ss system %P cpu %es total %MKb max memory"

for s in $( echo "$vcf_spec" | tr ',' ' ' );
do
    ref=$(echo "$s" | cut -f 1 -d: )
    samples=$(echo "$s" | cut -f 2 -d: )
    echo "[vg::deconstruct] making VCF with reference=$ref and samples=$samples"
    vcf="$gfa".smooth.$(echo $ref | tr '/|' '_').$(echo $samples | tr '/|' '_').vcf
    ( TEMPDIR=$(pwd) $timer -f "$fmt" vg deconstruct -P $ref \
         $(for i in $(cat $samples); do echo -n ' -A '$i; done) \
         -e -a -t $threads "$gfa" >"$vcf" )
    pigz -p $threads $vcf
done

rm $gfa
