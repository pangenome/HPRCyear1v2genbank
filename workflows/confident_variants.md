
## confident variants

Prepare tools:

```shell
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 120284fda67b3f3a67e8879c6ce7923923224f4f
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-120284fda67b3f3a67e8879c6ce7923923224f4f
```

Download and prepare unreliable regions (info at this [link](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components)):

```shell
mkdir -p /lizardfs/guarracino/HPRC/annotations/unreliable/
cd /lizardfs/guarracino/HPRC/annotations/unreliable/

cut -f 1 -d '#' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/hprcy1v2genbank.genomes.combined.fai | grep 'chm13\|grch38' -v | sort | uniq | while read SAMPLE; do
  wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/APR_08_2022/FINAL_HIFI_BASED/FLAGGER_HIFI_ASM_SIMPLIFIED_BEDS/$SAMPLE/$SAMPLE.hifi.flagger_final.simplified.unreliable_only.bed
done
```

Prepare reliable graphs:

```shell
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-120284fda67b3f3a67e8879c6ce7923923224f4f

mkdir -p /lizardfs/guarracino/HPRC/confident_variants/
cd /lizardfs/guarracino/HPRC/confident_variants/

( seq 1 22; echo X) | while read i; do
  echo chr$i
  
  mkdir -p chr$i
  cd chr$i
  
  PATH_GRAPH_OG_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/chr$i.pan/chr$i.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og.gz
  PREFIX=$(basename $PATH_GRAPH_OG_GZ .og.gz)

  # Take path information
  $RUN_ODGI paths -i <(zcat $PATH_GRAPH_OG_GZ) --list-paths --list-path-start-end > $PREFIX.paths.tsv

  # Take unreliable regions for the paths present in the input graph
  cat /lizardfs/guarracino/HPRC/annotations/unreliable/*.hifi.flagger_final.simplified.unreliable_only.bed | \
    grep -f <(grep chr -v $PREFIX.paths.tsv | cut -f 1) | bedtools sort > $PREFIX.unreliable.bed

  # Take the complement of the unreliable regions, obtaining the reliable regions
  # Add the future subpath names as 4th column in the BED file
  # -L limits the output to solely the chromosomes with records in the input file.
  bedtools complement -i $PREFIX.unreliable.bed -g <(grep chr -v $PREFIX.paths.tsv | cut -f 1,3 | sort) -L | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' > $PREFIX.reliable.subpaths.bed
  
  # Inject the reliable regions as subpaths, using the 4th as (sub)path names
  # Remove only full contigs of which there are reliable/unreliable regions
  $RUN_ODGI inject -i <(zcat $PATH_GRAPH_OG_GZ) -b $PREFIX.reliable.subpaths.bed -o - -t 48 -P | \
    $RUN_ODGI prune -i - -r <(cut -f 1 $PREFIX.reliable.subpaths.bed | sort | uniq) -o - -P | \
    $RUN_ODGI view -i - -g > $PREFIX.reliable.gfa

  cd ..
done
```

Call variants:

```shell
( seq 1 22; echo X) | while read i; do
  PATH_GRAPH_OG_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/chr$i.pan/chr$i.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og.gz
  PREFIX=$(basename $PATH_GRAPH_OG_GZ .og.gz)

  sbatch -p workers -c 24 --wrap "hostname; cd /scratch; \time -v vg-1.36.0 deconstruct -P chm13 -H '#' -e -a -t 24 /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.reliable.gfa | bgzip -c -@ 24 > $PREFIX.reliable.vcf.gz && tabix $PREFIX.reliable.vcf.gz && mv $PREFIX.reliable.vcf.gz* /lizardfs/guarracino/HPRC/confident_variants/chr$i/"
  sbatch -p workers -c 24 --wrap "hostname; cd /scratch; \time -v vg-1.40.0 deconstruct -P chm13 -H '#' -e -a -t 24 /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.reliable.gfa | bgzip -c -@ 24 > $PREFIX.reliable.vg_1_40_0.vcf.gz && tabix $PREFIX.reliable.vg_1_40_0.vcf.gz && mv $PREFIX.reliable.vg_1_40_0.vcf.gz* /lizardfs/guarracino/HPRC/confident_variants/chr$i/"
done
```

[comment]: <> (OLD BACKUP, TO IGNORE)
[comment]: <> (mkdir -p /lizardfs/guarracino/HPRC/confident_variants/)
[comment]: <> (cd /lizardfs/guarracino/HPRC/confident_variants/)
[comment]: <> (&#40; seq 13 15; seq 21 22&#41; | while read i; do)
[comment]: <> (mkdir -p chr$i)
[comment]: <> (cd chr$i)
[comment]: <> (PATH_GRAPH_OG_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/chr$i.pan/chr$i.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og.gz)
[comment]: <> (PREFIX=$&#40;basename $PATH_GRAPH_OG_GZ .og.gz&#41;)
[comment]: <> (odgi paths -i <&#40;zcat $PATH_GRAPH_OG_GZ&#41; -L > $PREFIX.paths.txt)
[comment]: <> (cut -f 1 -d '#' $PREFIX.paths.txt | sort | uniq | grep 'chm13\|grch38' -v | while read SAMPLE; do)
[comment]: <> (echo $SAMPLE)
[comment]: <> (    sbatch -p workers -c 48 ../call_confident_variants.sh $PATH_GRAPH_OG_GZ $PREFIX $SAMPLE chr$i $RUN_ODGI /lizardfs/guarracino/HPRC/confident_variants/chr$i)
[comment]: <> (    #sbatch -p workers -c 12 --wrap "\time -v vg deconstruct -P chm13 -H '#' -e -a -t 12 /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.gfa | sed 's/:[0-9]*-[0-9]*//g' | bgzip -c -@ 48 > /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.vcf.gz && tabix /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.vcf.gz")
[comment]: <> (    #vg deconstruct -P chm13 -H '#' -e -a -t 48 /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.gfa | sed 's/:[0-9]*-[0-9]*//g' | bgzip -c -@ 48 > /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.vcf.gz)
[comment]: <> (    #tabix /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.$SAMPLE.reliable.vcf.gz)
[comment]: <> (done)
[comment]: <> (cd ..)
[comment]: <> (done)
[comment]: <> (#call_confident_variants.sh)
[comment]: <> (#!/bin/bash)
[comment]: <> (PATH_GRAPH_OG_GZ=$1)
[comment]: <> (PREFIX=$2)
[comment]: <> (SAMPLE=$3)
[comment]: <> (CHR=$4)
[comment]: <> (RUN_ODGI=$5)
[comment]: <> (DIR_OUTPUT=$6)
[comment]: <> (grep "^$SAMPLE\|chm13" /lizardfs/erikg/HPRC/year1v2genbank/parts/$CHR.pan.fa.fai | awk -v OFS='\t' '{print&#40;$1,"0",$2&#41;}' > $PREFIX.$SAMPLE.bed)
[comment]: <> ($RUN_ODGI extract -i <&#40;zcat $PATH_GRAPH_OG_GZ&#41; -b $PREFIX.$SAMPLE.bed -p <&#40; cut -f 1 $PREFIX.$SAMPLE.bed | sort | uniq &#41; -o - -t 48 -P | \)
[comment]: <> ($RUN_ODGI unchop -i - -o - -t 48 -P | \)
[comment]: <> ($RUN_ODGI view -i - -g | sed 's/:[0-9]*-[0-9]*//g' > $PREFIX.$SAMPLE.gfa)
[comment]: <> (bedtools subtract \)
[comment]: <> (-a $PREFIX.$SAMPLE.bed \)
[comment]: <> (-b /lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed \)
[comment]: <> (> $PREFIX.$SAMPLE.reliable.bed)
[comment]: <> ($RUN_ODGI chop -i $PREFIX.$SAMPLE.gfa -c 1 -o - -t 48 -P | \)
[comment]: <> ($RUN_ODGI extract -i - -b $PREFIX.$SAMPLE.reliable.bed -p <&#40; cut -f 1 $PREFIX.$SAMPLE.reliable.bed | sort | uniq &#41; -o - -t 48 -P | \)
[comment]: <> ($RUN_ODGI unchop -i - -o - -t 48 -P | \)
[comment]: <> ($RUN_ODGI view -i - -g > $PREFIX.$SAMPLE.reliable.gfa)
[comment]: <> (vg deconstruct -P chm13 -H '#' -e -a -t 48 $PREFIX.$SAMPLE.reliable.gfa | sed 's/:[0-9]*-[0-9]*//g' | bgzip -c -@ 48 > $PREFIX.$SAMPLE.reliable.vcf.gz)
[comment]: <> (tabix $PREFIX.$SAMPLE.reliable.vcf.gz)
[comment]: <> (mv $PREFIX.$SAMPLE* $DIR_OUTPUT)
[comment]: <> (vcf2tsv $PATH_SNV_VCF_GZ -g | sed '1d' | awk '$17 != 0 && $17 != "0|0" && $17 != "."' > $PATH_SNV_VCF_GZ.bed)
