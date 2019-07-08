input_dir="$1"
output_dir="$2"
mkdir $output_dir
cd $input_dir

echo ">>> Performing IDR normalisation "
files_gr1=($(ls AS*K36me3*.narrowPeak))
files_gr2=($(ls AT*K36me3*.narrowPeak))
num_cases=${#files_gr1[@]}
for (( i=0; i<=$num_cases; i++))
do
    echo "> Samples: $input_dir/${files_gr1[$i]} + $input_dir/${files_gr2[$i]}"
    idr --plot --samples $input_dir/${files_gr1[$i]} $input_dir/${files_gr2[$i]} \
    --output-file $output_dir/${files_gr1[$i]%_peaks.narrowPeak}_vs_${files_gr2[$i]%_peaks.narrowPeak}_idr_values.txt
    #--use-best-multisummit-IDR \
done