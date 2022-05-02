# build log

## download and preprocessing

Get the URLs of the assemblies:

```shell
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
<Year1_assemblies_v2_genbank.index grep 'chm13\|h38' | awk '{ print $2 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' >refs.urls
<Year1_assemblies_v2_genbank.index grep -v 'chm13\|h38' | awk '{ print $2; print $3 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' >samples.urls
```

Download them:

```shell
mkdir assemblies
cd assemblies
cat ../refs.urls ../samples.urls | parallel -j 4 'wget -q {} && echo got {}'
```

Add a prefix to the reference sequences:

```shell
( fastix -p 'grch38#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) >grch38.fa && samtools faidx grch38.fa ) &
( fastix -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) >chm13.fa && samtools faidx chm13.fa ) &
wait
```

Combine them into a single reference for competitive assignment of sample contigs to chromosome bins:

```
cat chm13.fa grch38.fa >chm13+grch38_full.fa && samtools faidx chm13+grch38_full.fa
```

Remove unplaced contigs from grch38 that are (hopefully) represented in chm13:

```shell
samtools faidx chm13+grch38_full.fa $(cat chm13+grch38_full.fa.fai | cut -f 1 | grep -v _ ) >chm13+grch38.fa && samtools faidx chm13+grch38.fa
```

Unpack our assemblies:

```shell
ls *.gz | grep genbank | while read f; do sbatch -p lowmem -c 16 --wrap 'gunzip '$f' && samtools faidx '$(basename $f .gz); done >unpack.jobids
```

Manually break the only non-acrocentric misjoin:

```shell
samtools faidx HG02080.paternal.f1_assembly_v2_genbank.fa $( ( cat HG02080.paternal.f1_assembly_v2_genbank.fa.fai | cut -f 1 | grep -v HG02080#1#JAHEOW010000073.1 ; echo HG02080#1#JAHEOW010000073.1:0-7208112; echo HG02080#1#JAHEOW010000073.1:7208112-12869124 ) | sort -V ) | sed s/HG02080#1#JAHEOW010000073.1:0-7208112/HG02080#1#JAHEOW010000073.1_a/ | sed s/HG02080#1#JAHEOW010000073.1:7208112-12869124/HG02080#1#JAHEOW010000073.1_b/  >HG02080.paternal.f1_assembly_v2_genbank_split.fa && samtools faidx HG02080.paternal.f1_assembly_v2_genbank_split.fa
# prevents us from using this file in partitioning
pigz HG02080.paternal.f1_assembly_v2_genbank.fa
rm HG02080.paternal.f1_assembly_v2_genbank.fa.fai
```

Then we'll step into the directory below with `cd ..`.

## partitioning by chromosome

Partition the assembly contigs by chromosome by mapping each assembly against the scaffolded references, and then subsetting the graph. Here we use [wfmash](https://github.com/ekg/wfmash) for the mapping:

```shell
dir=approx_mappings
mkdir -p $dir
ref=assemblies/chm13+grch38.fa
aligner=/gnu/store/gsgm077q7krn9i08n15lshsw28paim14-wfmash-0.5.0+96d4426-12/bin/wfmash
for hap in $(cat haps.list);
do
    in=assemblies/$(ls assemblies | grep $hap | grep .fa$)
    out=$dir/$hap.vs.ref.paf
    sbatch -p lowmem -c 16 --wrap "$aligner -t 16 -m -N -s 50000 -p 90 $ref $in >$out" >>partition.jobids
done
```

Collect unmapped contigs and remap them in split mode:

```shell
dir=approx_mappings
ref=assemblies/chm13+grch38.fa
aligner=/gnu/store/gsgm077q7krn9i08n15lshsw28paim14-wfmash-0.5.0+96d4426-12/bin/wfmash  
for hap in $(cat haps.list);
do
    in=assemblies/$(ls assemblies | grep $hap | grep .fa$)
    paf=$dir/$hap.vs.ref.paf
    out=$dir/$hap.unaligned
    comm -23 <(cut -f 1 $in.fai | sort) <(cut -f 1 $paf | sort) >$out.txt
    if [[ $(wc -l $out.txt | cut -f 1 -d\ ) != 0 ]];
    then 
        samtools faidx $in $(tr '\n' ' ' <$out.txt) >$out.fa
        samtools faidx $out.fa
        sbatch -p lowmem -c 16 --wrap "$aligner -t 16 -m -s 50000 -p 90 $ref $out.fa >$out.split.vs.ref.paf" >>partition.jobids
    fi
    echo $hap
done
```

Collect our best mapping for each of our attempted split rescues.

```shell
dir=approx_mappings
ls $dir/*.unaligned.split.vs.ref.paf | while read f;
do
    cat $f | awk '{ print $1,$11,$0 }' | tr ' ' '\t' |  sort -n -r -k 1,2 | awk '$1 != last { print; last = $1; }'
done >$dir/rescues.paf
```

Subset by chromosome:

```shell
dir=approx_mappings
mkdir -p parts
( seq 22; echo X; echo Y; echo M ) | while read i; do awk '$6 ~ "chr'$i'$"' $(ls $dir/*.vs.ref.paf | grep -v unaligned | sort; echo $dir/rescues.paf) | cut -f 1 | sort >parts/chr$i.contigs; done
( seq 22; echo X; echo Y; echo M ) | while read i; do sbatch -p lowmem -c 16 --wrap './collect.sh '$i' >parts/chr'$i'.pan.fa && samtools faidx parts/chr'$i'.pan.fa' ; done >parts.jobids
# make a combined X+Y
cat parts/chrX.pan.fa parts/chrY.pan.fa >parts/chrS.pan.fa && samtools faidx parts/chrS.pan.fa
# make a combined acrocentric input
cat parts/chr{13,14,15,21,22}.pan.fa >parts/chrA.pan.fa && samtools faidx parts/chrA.pan.fa
```

This results in chromosome-specific FASTAs in `parts/chr*.pan.fa`.

## graph building

We now apply [pggb](https://github.com/pangenome/pggb):

```shell
( echo 1 16 2 3 4 5 6 7 8 X 9 10 11 12 13 14 15 17 18 19 20 21 22 Y | tr ' ' '\n') \
    | while read i; do sbatch -p workers -c 48 --wrap 'hostname; cd /scratch &&
    /gnu/store/2mjiai3jvwgi4cdw0cmjzpl3g97lb8lr-pggb-0.2.0+531f85f-1/bin/pggb
        -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -o chr'$i'.pan
        -t 48 -p 98 -s 100000 -n 90 -k 311 -O 0.03 -T 48
        -U -v -L -V chm13:#,grch38:# -Z ; mv /scratch/chr'$i'.pan '$(pwd);
    done >>pggb.jobids
```

Note that, for clarity, the command line given in quotes has been broken across multiple lines, which may cause problems if it is copied and pasted without editing.

A slightly different command line is used for the mitochondria, specifically we set `-s 1000` to improve sensitivity in this short chromosome.

## evaluation

Go in the `evaluation` directory:

```shell
cd evaluation
```

Download and prepare the reference:

```shell
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

run_rtg=/home/tools/RealTimeGenomics/3.12/rtg
$run_rtg format -o GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download the 'truth' set:

```shell
wget -c http://hypervolu.me/~guarracino/HPRC/HPRC_variant_calling_evaluation/HG00438.GRCh38_no_alt.deepvariant.vcf.gz
```

Download the Dipcall confident regions for the HG00438 sample:

```shell
wget -c http://hypervolu.me/~guarracino/HPRC/HPRC_variant_calling_evaluation/HG00438.f1_assembly_v2.dip.bed
```

<!---
```shell
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/original_vcfs/HG00438.GRCh38_no_alt.deepvariant.vcf.gz
```
-->

<!---
Download the easy/hard regions:

```shell
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/dipcall_truth/pggb-burned/GRCh38_notinalldifficultregions.bed.gz
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/dipcall_truth/pggb-burned/GRCh38_alldifficultregions.bed.gz
```
-->

Download the stratification files:

```shell
wget -r -nH --cut-dirs=6 ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38
```

- `-nH` avoids the creation of a directory named after the server name;
- `--cut-dirs=6` allows to put the content in the directory where you launch `wget`. The number 6 is used to filter out
  the 6-th components of the path.


Run the evaluations. For example, for chromosome 20, run:

```shell
grep chr20 HG00438.f1_assembly_v2.dip.bed | bgzip > HG00438.f1_assembly_v2.dip.chr20.bed.gz
tabix -p bed HG00438.f1_assembly_v2.dip.chr20.bed.gz
./vcf_evaluation.sh HG00438 chr20.pan.fa.c3d3224.7748b33.eb1aaa2.smooth.vcf.gz HG00438.f1_assembly_v2.dip.chr20.bed.gz GRCh38_notinalldifficultregions.bed.gz GRCh38_alldifficultregions.bed.gz HG00438_eval_out 16
```

The detailed results will be in `HG00438_eval_out/`.

## confident variants

Prepare tools:

```shell
git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout ed9390a47e6b029a753cbb83b29945eb48ca5c3b
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-ed9390a47e6b029a753cbb83b29945eb48ca5c3b

vg version
#vg version v1.40.0 "Suardi"
```

Prepare reliable graphs:

```shell
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-ed9390a47e6b029a753cbb83b29945eb48ca5c3b

mkdir -p /lizardfs/guarracino/HPRC/confident_variants/
cd /lizardfs/guarracino/HPRC/confident_variants/

# 
( seq 1 22; echo X; echo Y) | while read i; do
  echo chr$i
  
  mkdir -p chr$i
  cd chr$i
  
  PATH_GRAPH_OG_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/chr$i.pan/chr$i.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og.gz
  PREFIX=$(basename $PATH_GRAPH_OG_GZ .og.gz)

  # Take path information
  $RUN_ODGI paths -i <(zcat $PATH_GRAPH_OG_GZ) --list-paths --list-path-start-end > $PREFIX.paths.tsv

  # Take unreliable regions for the paths present in the input graph
  cat /lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/*.hifi.flagger_final.bed | \
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
( seq 1 22; echo X; echo Y) | while read i; do
  PATH_GRAPH_OG_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/chr$i.pan/chr$i.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og.gz
  PREFIX=$(basename $PATH_GRAPH_OG_GZ .og.gz)

  sbatch -p workers -c 24 --wrap "\time -v vg deconstruct -P chm13 -H '#' -e -a -t 24 /lizardfs/guarracino/HPRC/confident_variants/chr$i/$PREFIX.reliable.gfa | bgzip -c -@ 24 > $PREFIX.reliable.vcf.gz && tabix $PREFIX.reliable.vcf.gz && mv $PREFIX.reliable.vcf.gz* /lizardfs/guarracino/HPRC/confident_variants/chr$i/"
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
