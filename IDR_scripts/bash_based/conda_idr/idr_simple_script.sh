BASE_DIR=/ngs_share/scratch/pavlo/gctb/analysis/ChIP/normr
IDR_DIR=$BASE_DIR/idr
for mark in K4me1 K27ac K27me3 K36me3
do
#qsub -b y -e $IDR_DIR/AR_AS_${mark}_idr_cluster.log /opt/miniconda/bin/idr --samples $BASE_DIR/AR_${mark}_normr.peakNarrow $BASE_DIR/AS_${mark}_normr.peakNarrow --output-file $IDR_DIR/AR_AS_${mark}_idr.txt --log-output-file $IDR_DIR/AR_AS_${mark}_idr.log --plot
qsub -b y -e $IDR_DIR/AR_AT_${mark}_idr_cluster.log /opt/miniconda/bin/idr --samples $BASE_DIR/AR_${mark}_normr.peakNarrow $BASE_DIR/AT_${mark}_normr.peakNarrow --output-file $IDR_DIR/AR_AT_${mark}_idr.txt --log-output-file $IDR_DIR/AR_AT_${mark}_idr.log #--plot
#qsub -b y -e $IDR_DIR/GCTB11_KM882_${mark}_idr_cluster.log /opt/miniconda/bin/idr --samples $BASE_DIR/GCTB11_${mark}_normr.peakNarrow $BASE_DIR/KM882_${mark}_normr.peakNarrow --output-file $IDR_DIR/GCTB11_KM882_${mark}_idr.txt --log-output-file $IDR_DIR/GCTB11_KM882_${mark}_idr.log --plot
done