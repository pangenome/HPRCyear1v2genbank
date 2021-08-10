#!/bin/bash

i=$1
m=$2

/gnu/store/s5a40kp4w9nzjypb029vmbgygy5ycnvi-odgi-0.6.0+84de376-1/bin/odgi untangle \
	-t 16 -i <(zcat /lizardfs/erikg/HPRC/year1v2genbank/wgg.75/chr$i.pan/*.smooth.og.gz) \
	-m $m -r chm13#chr$i -q chm13#chr$i >chr$i.pan.chm13.untangle-m$m.bed
