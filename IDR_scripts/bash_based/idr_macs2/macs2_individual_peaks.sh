input_dir="$1"
output_dir="$2"
#init_dir=$(pwd)
mkdir $output_dir
#cd $output_dir
#output_dir=$(pwd)
#cd $init_dir
cd $input_dir
#input_dir=$(pwd)

files=($(ls A*K36me3*.bam))
for file in ${files[@]}
do
    /ngs_share/tools/miniconda3/envs/ATAC/bin/macs2 callpeak --bdg --keep-dup all --name ${file%.bam} --format BAM -g hs --verbose 3 --pvalue 0.01 -t $input_dir/$file --outdir $output_dir
    #qsub -b y -e $output_dir/$file.log /ngs_share/tools/miniconda3/envs/ATAC/bin/macs2 callpeak --bdg --keep-dup all --name ${file%.bam} --format BAM -g hs --verbose 3 --pvalue 0.0001 -t $input_dir/$file --outdir $output_dir
done
