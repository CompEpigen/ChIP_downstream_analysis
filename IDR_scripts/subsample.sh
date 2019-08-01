ROOT_DIR="/Users/adams/Documents/Heidelberg/IDR"
comparison_file="/Users/adams/Documents/Heidelberg/IDR/comparison_files/comparison_file_subsample.csv"
lines=`tail -n +2 $comparison_file`
for ((i=0;i<10;i++))
do
for line in ${lines[@]}
do
    IFS=';' read -r -a array <<< "$line"
    lab="${array[0]}"
    mark="${array[1]}"
    control="${array[2]}"
    shuf -n 7215 $ROOT_DIR/${lab}/${mark}/${control}/peak1 > $ROOT_DIR/subsample_Plass/${mark}/${control}/peak1_subsample${i}
    shuf -n 3751 $ROOT_DIR/${lab}/${mark}/${control}/peak2 > $ROOT_DIR/subsample_Plass/${mark}/${control}/peak2_subsample${i}
done
done

for ((i=0;i<10;i++))
do
for line in ${lines[@]}
do
    IFS=';' read -r -a array <<< "$line"
    lab="${array[0]}"
    mark="${array[1]}"
    control="${array[2]}"
    idr --samples $ROOT_DIR/subsample_Plass/${mark}/${control}/peak1_subsample${i} $ROOT_DIR/subsample_Plass/${mark}/${control}/peak2_subsample${i} --input-file-type broadPeak --plot --output-file $ROOT_DIR/subsample_Plass/${mark}_output/${control}/subsample${i}_${mark}
done
done
