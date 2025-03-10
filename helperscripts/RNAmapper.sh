#!/usr/bin/env bash

dir=$1
prefix=$2
out=${3-$dir}${prefix}
chrs=$(ls $dir | grep -v "/" | grep -E ".*mut_chr([0-9]+)\\.vcf" | sed -E "s/.*mut_chr([0-9]+)\\.vcf/\\1/" | sort -n -k1)

echo "RNAmapper.py running on files ${dir}${prefix}..."

for chr in $chrs
do
    wt="${dir}${prefix}_wt_chr${chr}.vcf"
    mut="${dir}${prefix}_mut_chr${chr}.vcf"
    ../RNAmapper.py -wt $wt -mut $mut -o $out$prefix -ch $chr &
done

wait
echo "RNAmapper.py completed on ${dir}${prefix}!"