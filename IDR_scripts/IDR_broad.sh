ls
ROOT_DIR="/Users/adams/Documents/Heidelberg/IDR"
comparison_file="comparison_file_broad.csv"
lines=`tail -n +2 $comparison_file`
for line in ${lines[@]}
do
    IFS=';' read -r -a array <<< "$line"
    lab="${array[0]}"
    mark="${array[1]}"
    control="${array[2]}"
    idr --samples $ROOT_DIR/${lab}/${mark}/${control}/peak1 $ROOT_DIR/${lab}/${mark}/${control}/peak2 --input-file-type broadPeak --plot --output-file $ROOT_DIR/output/${lab}_${mark}_${control}_output_broad
done