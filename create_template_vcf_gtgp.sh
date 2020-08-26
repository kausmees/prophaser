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
	indlist="${indlist}"'\t1/1:0,0,1'"";
done

echo $indlist

bcftools view -h $indir$vcf_filestart.vcf.gz > $indir$vcf_filestart.metainfo.gtgp
# Write metainfo lines (all but the last one)
sed -e '$ d' $indir$vcf_filestart.metainfo.gtgp > $indir$vcf_filestart.template_gtgp.vcf
# append the GT and GP fields format to header
echo -e '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> $indir$vcf_filestart.template_gtgp.vcf
echo -e '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">' >> $indir$vcf_filestart.template_gtgp.vcf
# append header line
bcftools view -h $indir$vcf_filestart.vcf.gz | grep '#CHROM' >> $indir$vcf_filestart.template_gtgp.vcf
# append VCF body
echo -e  "20	2355080	rs566601297	G	A	.	PASS	.	GT:GP"$indlist >> $indir$vcf_filestart.template_gtgp.vcf

rm $indir$vcf_filestart.metainfo.gtgp
echo "wrote to " $indir$vcf_filestart.template_gtgp.vcf


