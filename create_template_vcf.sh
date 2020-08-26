indir=$1
vcf_file=$2
vcf_filestart=${vcf_file::-7}

num_words=`bcftools view -h $indir$vcf_filestart.vcf.gz | tail -n 1 | wc -w`

num_inds=`expr $num_words - 9 `

echo $num_inds

indlist=""

for i in $(seq 1 $num_inds);
do
	echo $i
	indlist="${indlist}"'\t1/1'"";
done

echo $indlist

bcftools view -h $indir$vcf_filestart.vcf.gz > $indir$vcf_filestart.metainfo
# Write metainfo lines (all but the last one)
sed -e '$ d' $indir$vcf_filestart.metainfo > $indir$vcf_filestart.template.vcf
# append the GT field format to header
echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $indir$vcf_filestart.template.vcf
# append header line
bcftools view -h $indir$vcf_filestart.vcf.gz | grep '#CHROM' >> $indir$vcf_filestart.template.vcf
# append VCF body
echo -e  "20	2355080	rs566601297	G	A	.	PASS	.	GT"$indlist >> $indir$vcf_filestart.template.vcf

rm $indir$vcf_filestart.metainfo
echo "wrote to " $indir$vcf_filestart.template.vcf


