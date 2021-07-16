# build log

Get the URLs of the assemblies:

```
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index
<Year1_assemblies_v2_genbank.index grep 'chm13\|h38' | awk '{ print $2 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' >refs.urls
<Year1_assemblies_v2_genbank.index grep -v 'chm13\|h38' | awk '{ print $2; print $3 }' | sed 's%s3://human-pangenomics/working/%https://s3-us-west-2.amazonaws.com/human-pangenomics/working/%g' >samples.urls
```

Download them:

```
mkdir assemblies
cd assemblies
cat ../refs.urls ../samples.urls | parallel -j 4 'wget -q {} && echo got {}'
```

Add a prefix to the reference sequences:

```
( fastix -p 'grch38#' <(zcat GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) >grch38.fa && samtools faidx grch38.fa ) &
( fastix -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) >chm13.fa && samtools faidx chm13.fa ) &
wait
```

Combine them into a single reference for competitive assignment of sample contigs to chromosome bins:

```
cat chm13.fa grch38.fa >chm13+grch38_full.fa && samtools faidx chm13+grch38_full.fa
```

Remove unplaced contigs from grch38 that are (hopefully) represented in chm13:

```
samtools faidx chm13+grch38_full.fa $(cat chm13+grch38_full.fa.fai | cut -f 1 | grep -v _ ) >chm13+grch38.fa && samtools faidx chm13+grch38.fa
```

Unpack our assemblies:

```
ls *.gz | grep genbank | while read f; do sbatch -p lowmem -c 16 --wrap 'gunzip '$f' && samtools faidx '$(basename $f .gz); done >unpack.jobids
```

Manually break the only non-acrocentric misjoin:

```
samtools faidx HG02080.paternal.f1_assembly_v2_genbank.fa $( ( cat HG02080.paternal.f1_assembly_v2_genbank.fa.fai | cut -f 1 | grep -v HG02080#1#JAHEOW010000073.1 ; echo HG02080#1#JAHEOW010000073.1:0-7208112; echo HG02080#1#JAHEOW010000073.1:7208112-12869124 ) | sort -V ) | sed s/HG02080#1#JAHEOW010000073.1:0-7208112/HG02080#1#JAHEOW010000073.1_a/ | sed s/HG02080#1#JAHEOW010000073.1:7208112-12869124/HG02080#1#JAHEOW010000073.1_b/  >HG02080.paternal.f1_assembly_v2_genbank_split.fa && samtools faidx HG02080.paternal.f1_assembly_v2_genbank_split.fa
# prevents us from using this file in partitioning
pigz HG02080.paternal.f1_assembly_v2_genbank.fa
rm HG02080.paternal.f1_assembly_v2_genbank.fa.fai
```

Then we'll step into the directory below with `cd ..`.

Partition the assembly contigs by chromosome by mapping each assembly against the scaffolded references, and then subsetting the graph. Here we use [wfmash](https://github.com/ekg/wfmash) for the mapping:

```
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

```
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

```
dir=approx_mappings
ls $dir/*.unaligned.split.vs.ref.paf | while read f;
do
    cat $f | awk '{ print $1,$11,$0 }' | tr ' ' '\t' |  sort -n -r -k 1,2 | awk '$1 != last { print; last = $1; }'
done >$dir/rescues.paf
```

Subset by chromosome:

```
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

We now apply [pggb](https://github.com/pangenome/pggb):

```
# nuclear chromosomes
( seq 8; echo X; seq 9 22; echo Y ) | while read i; do sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && pggb -t 48 -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -p 98 -s 100000 -n 70 -k 79 -B 10000000 -w 1000000 -G 6733,9929 -P 1,19,39,3,81,1 -v -S -L -o chr'$i'.pan -Z ; mv /scratch/chr'$i'.pan '$(pwd); done >>pggb.jobids
# chrM
( echo M ) | while read i; do sbatch -p lowmem -c 48 --wrap 'cd /scratch && pggb --resume -t 48 -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -Y "#" -p 98 -s 1000 -l 3000 -n 70 -k 79 -B 10000000 -w 10000000 -G 10000 -v -S -L -o chr'$i'.pan -Z ; mv /scratch/chr'$i'.pan '$(pwd); done >>pggb.jobids
# acrocentric / sex chromosome communities
( echo A; echo S ) | while read i; do sbatch -p lowmem -c 48 --wrap 'hostname; cd /scratch && pggb -t 48 -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -p 98 -s 100000 -n 70 -k 79 -B 10000000 -w 1000000 -G 6733,9929 -P 1,19,39,3,81,1 -v -S -L -o chr'$i'.pan -Z ; mv /scratch/chr'$i'.pan '$(pwd); done >>pggb.jobids
```

## evaluation

Go in the `evaluation` directory:

```
cd evaluation
```

Download and prepare the reference:

```
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

run_rtg=/home/tools/RealTimeGenomics/3.12/rtg
$run_rtg format -o GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sdf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download the 'truth' set:

```
wget -c http://hypervolu.me/~guarracino/HPRC_variant_calling_evaluation/HG00438.GRCh38_no_alt.deepvariant.vcf.gz
```

Download the Dipcall confident regions for the HG00438 sample:

```
wget -c http://hypervolu.me/~guarracino/HPRC_variant_calling_evaluation/HG00438.f1_assembly_v2.dip.bed
```

<!---
```
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/original_vcfs/HG00438.GRCh38_no_alt.deepvariant.vcf.gz
```
-->

<!---
Download the easy/hard regions:

```
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/dipcall_truth/pggb-burned/GRCh38_notinalldifficultregions.bed.gz
wget -c https://9a3fe.03c0.data.globus.org/benchmarking/dipcall_truth/pggb-burned/GRCh38_alldifficultregions.bed.gz
```
-->

Download the stratification files:
```
wget -r -nH --cut-dirs=6 ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38
```
- `-nH` avoids the creation of a directory named after the server name;
- `--cut-dirs=6` allows to put the content in the directory where you launch `wget`. The number 6 is used to filter out
  the 6-th components of the path.


Run the evaluations. For example, for chromosome 20, run:

```
grep chr20 HG00438.f1_assembly_v2.dip.bed > HG00438.f1_assembly_v2.dip.chr20.bed
grep union v2.0-GRCh38-stratifications.tsv | grep difficult > v2.0-GRCh38-stratifications.easy_hard.tsv

bash vcf_evaluation.sh HG00438 chr20.pan.fa.c3d3224.7748b33.eb1aaa2.smooth.vcf.gz HG00438.f1_assembly_v2.dip.chr20.bed v2.0-GRCh38-stratifications.easy_hard.tsv HG00438_eval_out 16
```