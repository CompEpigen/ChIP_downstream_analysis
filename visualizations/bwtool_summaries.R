library(GenomicRanges)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)


ROOT_DIR="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/"

SAMPLE_WISE<-FALSE
############ read bwtool matrices
marks<-c("WGBS",
        "ATAC", 
        "K4me1", "K4me3",  "K27ac","K27me3", "K36me3",
        "K9me3",
        "H3.3", "G34W", "H3",
        "RNAseq", "CAGE")#[c(1:6,10)]

#marks<-c("G34W")

mt_wt_sample_groups<-list(
        "WT"=c("GCTB11", "GCTB_11", "GCTB12", "GCTB_12", "AO", "AZ","KM882","KM931", "KM1404", "KM1234"),
        ### KM921 has too few reads
        #"WT"=c("GCTB11", "GCTB_11", "GCTB12", "GCTB_12", "AO", "AZ","KM882", "KM931", "KM1404", "KM1234"),
        "MT"=c("AP","AR","AS","AT","AU","AW", "AX2"))
mt_wt_sample_groups$MT<-c(mt_wt_sample_groups$MT, sprintf("GCTB_%s", mt_wt_sample_groups$MT))

#
#wt_mt_gctb_sample_groups<-list(
#        "MSC"=c("KM882","KM921", "KM931", "KM1404"),
#        "GCTB_WT"=c("GCTB11", "AO", "AZ"),
#        "MT"=c("AP","AR","AS","AT","AU","AW", "AX2"))
#wt_mt_gctb_sample_groups$MT<-c(wt_mt_gctb_sample_groups$MT, sprintf("GCTB_%s", wt_mt_gctb_sample_groups$MT))


#mt_wt_sample_groups<-list(
#        "HEK293T-utr"=c("HEK293T-utr"),
#        "HEK293T-G34W"=c("HEK293T-G34W"),
#        "HEK293T-WT"=c("HEK293T-WT"),
#        "KM882"=c("KM882"),
#        "AO"=c("AO"),
#        "AR"=c("AR"),
#        "AS"=c("AS")
#### KM921 has too few reads
##"WT"=c("GCTB11", "GCTB_11", "GCTB12", "GCTB_12", "AO", "AZ","KM882", "KM931", "KM1404", "KM1234"),
#)

if(SAMPLE_WISE){
    sample_groups<-unlist(mt_wt_sample_groups)
    names(sample_groups)<-sample_groups
    sample_groups<-tapply(sample_groups, gsub("GCTB_|GCTB", "", sample_groups), as.list)
    sgn<-names(sample_groups); sgn[grep("^[0-9]+$", sgn)]<-paste0("GCTB", sgn[grep("^[0-9]+$", sgn)]); names(sample_groups)<-sgn
}else{
    sample_groups<-mt_wt_sample_groups
}



mut_wt_cols<-ggsci:::pal_npg()(10)[c(4,1)]
group_colors<-c("WT"=mut_wt_cols[1],"MT"=mut_wt_cols[2])

hek_cols<-ggsci:::pal_npg()(10)[c(4,5,9,1,2,3,7)]
group_colors<-setNames(hek_cols, names(sample_groups))


#seq_depths_chip<-read.csv("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/qc/ChIP_read_counts_raw.csv", sep=";")
seq_depths_chip<-read.csv("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/qc/ChIP_read_counts_duprem.csv", sep=";")

#seq_depths_atac<-read.csv("/ngs_share/scratch/pavlo/gctb/analysis/ATAC/qc/ATAC_read_counts.csv", sep=";")
seq_depths_atac<-as.matrix(c("AR"=1, "AS"=1, "AT"=1, "AO"=NA, "GCTB11"=1, "KM882"=NA, "KM1404"=1, "KM931"=1, "KM921"=1))
seq_depths_atac<-data.frame(Sample=rownames(seq_depths_atac), ATAC=seq_depths_atac)
seq_depths<-merge(seq_depths_chip, seq_depths_atac, by="Sample", all=TRUE)

#ANALYSIS="Basic_annotations"
#ANALYSIS="Bivalent"

#ANALYSIS="ESC_Bivalent"
#ANALYSIS="ESC_Bivalent_unstranded_3clusters"
#ANALYSIS="ESC_Bivalent_3groups"
#ANALYSIS="MSC_Bivalent_3groups"




#ANALYSIS="Bivalent_ESC_TSS"
#ANALYSIS="Bivalent_MSC_TSS"

#ANALYSIS="Bivalent_MSC_TSS_strat"


#ANALYSIS="Bivalent_ESC_full"
#ANALYSIS="Bivalent_MSC_full"
#ANALYSIS="Bivalent_MSC_stringent"

#ANALYSIS="Bivalent_MSC_TSS_sample_wise"


#ANALYSIS="MD_full"

#ANALYSIS="ATAC_diff"
#ANALYSIS="ENCODE_ChromHMM"

#ANALYSIS="ENCODE_ChromHMM_3groups"
#ANALYSIS="ENCODE_ChromHMM_signalNorm"

#ANALYSIS="ATAC_stratified"
#ANALYSIS="ATAC_stratified_3groups"
#ANALYSIS="ATAC_stratified_TSS"


#ANALYSIS="ChIP_union_peaks"
#ANALYSIS="ChIP_union_peaks_5k"

#ANALYSIS="Enhancers"
#ANALYSIS="Superenhancers"


#ANALYSIS="ChromHMM_WT_control_3groups"
#ANALYSIS="TeloCentroAll"
#ANALYSIS="TeloCentroStarts"
#ANALYSIS="TeloCentroEnds"
#ANALYSIS="H3.3_regions"
#ANALYSIS="H3.3_regions_OFFS_0"
#ANALYSIS="H3.3_regions_no_control_OFFS_0_rep2"
#ANALYSIS="H3.3_regions_wt_bg_OFFS_0_rep1_withY"
#ANALYSIS="H3.3_regions_selected"

#ANALYSIS="G34W_regions_old_redone"
#ANALYSIS="G34W_regions_new"
#ANALYSIS="G34W_regions_new_short"
#ANALYSIS="G34W_regions_new_long2"
#ANALYSIS="G34W_regions_new_long2_k9"
#ANALYSIS="G34W_SUZ12_overlap"

ANALYSIS="G34W_regions_final"

#ANALYSIS="CT_genes_new"
#ANALYSIS="CT_genes_full"
#ANALYSIS="CT_genes_TSS"

#ANALYSIS="PRC2_genes_TSS"
#ANALYSIS="PRC2_genes_full"

#ANALYSIS="CGIs"

#ANALYSIS="Repeats_simple"
#ANALYSIS="Repeats_simple_long"
#ANALYSIS="Repeats_TAR1"


#ANALYSIS="DE_genes_TSS"

#ANALYSIS="All_genes_TSS"
#ANALYSIS="All_genes_full"


#ANALYSIS="TFBS_TR4"
#ANALYSIS="TFBS_EZH2"
#ANALYSIS="TFBS_ZNF143"
#ANALYSIS="TFBS_SUZ12"

#ANALYSIS="RI_sites"
#ANALYSIS="RI_sites_binned"


#### for bwtools
ANNOTATION_DIRS<-list()
#ANNOTATION_DIR="/ngs_share/scratch/pavlo/gctb/analysis/integrative/general_annotations/"
#ANNOTATION_DIR="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/basic"
#ANNOTATION_DIR="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/ATAC_stratified"
#ANNOTATION_DIR="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CHMM_ENCODE"
ANNOTATION_DIRS[["ATAC_stratified"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/ATAC_stratified"
ANNOTATION_DIRS[["ATAC_stratified_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/ATAC_stratified_by_TSS"

ANNOTATION_DIRS[["MD_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/PMDs/new/"

ANNOTATION_DIRS[["CGIs"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CGIs"

ANNOTATION_DIRS[["ENCODE_ChromHMM"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CHMM_ENCODE"
ANNOTATION_DIRS[["ChromHMM_WT_control_3groups"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CHMM_WT/WT_merged_control_200"

ANNOTATION_DIRS[["Bivalent_ESC_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/ESC/TSS"
ANNOTATION_DIRS[["Bivalent_MSC_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/TSS"
ANNOTATION_DIRS[["Bivalent_MSC_TSS_strat"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/TSS/strat_up_down/"

ANNOTATION_DIRS[["ESC_Bivalent_unstranded"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/ESC/TSS_unstranded/"
ANNOTATION_DIRS[["ESC_Bivalent_unstranded_3clusters"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/ESC/TSS_unstranded_raw/"

ANNOTATION_DIRS[["Bivalent_ESC_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/ESC/full"
ANNOTATION_DIRS[["Bivalent_MSC_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/full"
ANNOTATION_DIRS[["Bivalent_MSC_stringent"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/TSS_stringent/"

ANNOTATION_DIRS[["Bivalent_MSC_stringent"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/TSS_stringent/"

#ANNOTATION_DIRS[["Bivalent_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent_full"

ANNOTATION_DIRS[["Bivalent_MSC_TSS_sample_wise"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Bivalent/MSC/TSS"

ANNOTATION_DIRS[["TeloCentroAll"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/general_annotations/cytoband/"
ANNOTATION_DIRS[["TeloCentroStarts"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/general_annotations/cytoband/starts"
ANNOTATION_DIRS[["TeloCentroEnds"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/general_annotations/cytoband/ends"

ANNOTATION_DIRS[["ChIP_union_peaks"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/ChIP/union_peaks_macs_fdr001/"
ANNOTATION_DIRS[["ChIP_union_peaks_5k"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/ChIP/union_peaks_macs_fdr001_rand5k/"

ANNOTATION_DIRS[["Enhancers"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/enhancers/"
ANNOTATION_DIRS[["Superenhancers"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/enhancers/super/"

ANNOTATION_DIRS[["H3.3_regions"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/H3_variants_all_control_200bp_0/"
ANNOTATION_DIRS[["H3.3_regions_OFFS_0"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/H3_variants_all_control_200bp_0/"
ANNOTATION_DIRS[["H3.3_regions_no_control_OFFS_0_rep2"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/H3_variants_all_200bp_0_rep1"
ANNOTATION_DIRS[["H3.3_regions_wt_bg_OFFS_0_rep1_withY"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks//H3_variants_all_control_bg_200bp_0_rep1"
ANNOTATION_DIRS[["H3.3_regions_selected"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/"

ANNOTATION_DIRS[["G34W_regions_old"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/Homer_g34w_old"
ANNOTATION_DIRS[["G34W_regions_old_redone"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/Homer_g34w_old"

ANNOTATION_DIRS[["G34W_regions_new"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/H3_variants_all_bg_200bp_0_top5000/"
ANNOTATION_DIRS[["G34W_regions_new_short"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/H3_variants_all_bg_200bp_0_top5000_g34_only/"
ANNOTATION_DIRS[["G34W_regions_new_long2"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/H3_variants_all_bg_200bp_0_rep2_full"
ANNOTATION_DIRS[["G34W_regions_new_long2_k9"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/H3_variants_all_bg_200bp_0_rep2_full"

ANNOTATION_DIRS[["G34W_regions_final"]]="/ngs_share/scratch/pavlo/gctb/analysis/ChIP/G34W/tracks/selected/g34w_final"

ANNOTATION_DIRS[["G34W_SUZ12_overlap"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/G34W/overlaps/SUZ12/"

ANNOTATION_DIRS[["CT_genes"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CT_genes/"
ANNOTATION_DIRS[["CT_genes_new"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CT_genes/"
ANNOTATION_DIRS[["CT_genes_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CT_genes/full"
ANNOTATION_DIRS[["CT_genes_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/CT_genes/TSS"

ANNOTATION_DIRS[["PRC2_genes"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/PRC2_target_genes/"

ANNOTATION_DIRS[["PRC2_genes_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/PRC2_target_genes/Full"
ANNOTATION_DIRS[["PRC2_genes_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/PRC2_target_genes/Full"

ANNOTATION_DIRS[["Repeats_simple"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Simple_repeats"
ANNOTATION_DIRS[["Repeats_simple_long"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Simple_repeats/long"
ANNOTATION_DIRS[["Repeats_TAR1"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/Telomere_repeats/"

ANNOTATION_DIRS[["DE_genes_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/DE_genes/TSS/"

ANNOTATION_DIRS[["All_genes_TSS"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/All_genes_kallisto/top5000/TSS/"
ANNOTATION_DIRS[["All_genes_TES"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/All_genes_kallisto/top5000/TES/"
ANNOTATION_DIRS[["All_genes_full"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/All_genes_kallisto/top5000/full/"


ANNOTATION_DIRS[["TFBS_TR4"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/TFBS/TR4"
ANNOTATION_DIRS[["TFBS_EZH2"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/TFBS/EZH2"
ANNOTATION_DIRS[["TFBS_ZNF143"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/TFBS/ZNF143"
ANNOTATION_DIRS[["TFBS_SUZ12"]]="/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/TFBS/SUZ12"

ANNOTATION_DIRS[["RI_sites"]]="/ngs_share/scratch/pavlo/gctb/analysis/RNA-seq/rMATS/tracks/"
ANNOTATION_DIRS[["RI_sites_binned"]]="/ngs_share/scratch/pavlo/gctb/analysis/RNA-seq/rMATS/tracks/"

        #"/ngs_share/scratch/pavlo/gctb/analysis/integrative/important_annotations/G34W_comb/unique/"

ANNOTATION_DIR<-ANNOTATION_DIRS[[ANALYSIS]]

#<-file.path(ROOT_DIR, "summaries", ANALYSIS)

#OFFSET<-20000
#OFFSET<-10000
#OFFSET<-7500
OFFSET<-5000
#OFFSET<-2500
#OFFSET<-5000
#OFFSET<-1000
#OFFSET<-0
WIN_SIZE<-100
N_BINS<-150
MAX_K<-33
BINNED<-FALSE

ANCHOR_SPECS=c("start", "mid", "end")[2]
#ANCHOR_NAME="TSS"
ANCHOR_NAME="TSS"

BINNED_REG_START_NAME<-"TSS"
BINNED_REG_END_NAME<-"TES"

#BINNED_REG_START_NAME<-"start"
#BINNED_REG_END_NAME<-"end"

if(BINNED){
    n_flank_win<-OFFSET/WIN_SIZE
    nwin<-2*OFFSET/WIN_SIZE + N_BINS
    bin_locs<-c((seq(-n_flank_win, 0.5)-0.5)[-1]*WIN_SIZE, 0, 1/seq(N_BINS)[-N_BINS] ,(seq(0, n_flank_win)-0.5)[-1]*WIN_SIZE)
    #bin_locs<-seq(0, nwin)[-1]
}else{
    nwin<-2*OFFSET/WIN_SIZE
    bin_locs<-(seq(-nwin/2, nwin/2)-0.5)[-1]*WIN_SIZE
}
NORMALIZE_ATAC<-TRUE
NORMALIZE_CHIP<-FALSE

#STATE_RGNS_DIR<-file.path(ROOT_DIR, "summaries", "Bivalent")
#STATE_RGNS_DIR<-file.path(ROOT_DIR, "summaries", "CHMM_ENCODE")
#STATE_RGNS_DIR<-file.path(ROOT_DIR, "summaries", "RNAseq")
#STATE_RGNS_DIR<-file.path(ROOT_DIR, "summaries", "Bivalent")

STATE_RGNS_DIR=file.path("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/", 
        paste(ANALYSIS, if(BINNED) "binned" else "centered", OFFSET, WIN_SIZE, sep="_"))


OUT_DIR=STATE_RGNS_DIR
if(!file.exists(OUT_DIR)){
    dir.create(OUT_DIR)
}

#features<-c("RNAseq_expressed_WT","RNAseq_silent_WT")

#features<-c("RNAseq_silent_WT","RNAseq_present_WT","RNAseq_active_WT")[c(1,3)]
#features<-c("RNAseq_activated","RNAseq_repressed")
#features<-c("RNAseq_full_fully_silent","RNAseq_full_low_expressed", "RNAseq_full_med_expressed", "RNAseq_full_high_expressed")
#features<-c("RNAseq_tss_fully_silent","RNAseq_tss_low_expressed", "RNAseq_tss_med_expressed", "RNAseq_tss_high_expressed",
#        "RNAseq_tes_fully_silent","RNAseq_tes_low_expressed", "RNAseq_tes_med_expressed", "RNAseq_tes_high_expressed") 
#features<-c("ATAC_gained_peaks", "ATAC_lost_peaks")

#chmm_annot<-read.table("/ngs_share/scratch/pavlo/gctb/data/annotations/chromhmm/ENCODE/chmm_15_state_model_states.txt", sep="\t", header=TRUE)
#features<-gsub("/", "_", chmm_annot$MNEMONIC)

#features<-c("bivalent_ESC_cluster_1","bivalent_ESC_cluster_2","bivalent_ESC_cluster_3")

#features<-sprintf("bivalent_ESC_cluster_combined_category_%s", readRDS("/ngs_share/scratch/pavlo/gctb/data/annotations/biv_genes_human_ESCs_combined_categories.RDS"))
#features<-features[-grep("Cluster1_CGINo_gene_body", features)]

#pdf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/summaries/all_atac_peaks_lim.pdf", width=20, height=20)
#layout(matrix(1:42, ncol=length(unlist(sample_groups)), byrow=FALSE))
#layout(cbind(seq(MAX_K), matrix(seq(MAX_K+1,(length(marks)+2)*MAX_K), ncol=length(marks)+1, byrow=TRUE)))
#layout(cbind(seq(length(features)),matrix(seq(length(features)+1,(length(marks)+2)*length(features)), ncol=length(marks)+1, byrow=TRUE)))

#for(feature in features){
#    plot.new()
#    text(0.5,0.5,sprintf("%s", gsub("RNAseq_|_expressed","",feature)), cex=1.5, srt=90)
#}
#BW_DIR="/C010-datasets/Internal/GCTB/ChIP/tracks/counts/simple_bamCoverage/10/bw"
#BW_DIR="/C010-datasets/Internal/GCTB/ATAC-Seq/tracks"
#BW_DIR="/C010-datasets/Internal/GCTB/WGBS/analysis/results/pilot_summarized/covg5_bw"

features<-gsub("\\.bed", "", list.files(ANNOTATION_DIR, pattern="*.bed$"))

if(TRUE){
    BW_DIRS<- c(
            "WGBS"="/C010-datasets/Internal/GCTB/WGBS/analysis/results/pilot_summarized/smoothed/bw",
            "ATAC"="/C010-datasets/Internal/GCTB/ATAC-Seq/tracks",
            #WGBS="/C010-datasets/Internal/GCTB/WGBS/analysis/results/pilot_summarized/covg5_bw",
            #"ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/counts/simple_bamCoverage/10/bw",
            #"ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/signalNormalized/bw",
            #"ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/RPKM/10/bw"
            #"ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/bamcompareRPKM/10/bw",
            "ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/seqDepthNorm/10/bw/",
            #"ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/macs2/input",
           "RNAseq"="/C010-datasets/Internal/GCTB/RNA-Seq/tracks/read_cov/bw/",
           "CAGE"="/C010-datasets/Internal/GCTB/CAGE-Seq/tracks"
    )
}else{
    BW_DIRS<-c("ChIP"="/C010-datasets/Internal/GCTB/ChIP/tracks/epicwl/bw")
}



for(BW_DIR in BW_DIRS){

    bw_files<-gsub("\\.bw", "", list.files(BW_DIR, pattern="bw$"))
    
    #nsamp<-length(bw_files)
    
    ########### extract profile matrices with BWtools
    
    command<-paste(
        sprintf("ANALYSIS=%s\n",ANALYSIS),
        sprintf("ANNOTATION_DIR=%s\n",ANNOTATION_DIR),
        sprintf("BW_DIR=%s\n",BW_DIR),
        sprintf("OUT_DIR=%s\n",OUT_DIR),
        sprintf("UPSTREAM=%d\n",OFFSET),
        sprintf("ANCHOR=%s\n", if(ANCHOR_SPECS=="start") "-starts" else if(ANCHOR_SPECS=="end") "-ends" else ""),
        sprintf("DOWNSTREAM=%d\n",OFFSET),
        sprintf("TILE_AVG=%d\n",WIN_SIZE),
        sprintf("TMP_DIR=%s\n","/ngs_share/tmp"),
        sprintf("REGIONS=\"%s\"\n", paste(features, collapse=" ")),
        sprintf("TOKEN=%s\n", if(BINNED) "binned" else ANCHOR_SPECS),
        sprintf("META=%s\n", N_BINS ),
        "mkdir $TMP_DIR/logs/$ANALYSIS\n",
        "INPUT_BWS=`ls -1 $BW_DIR`\n",
        "for rgn in `echo $REGIONS`\n",
        "do\n",
        "for bwig in `echo $INPUT_BWS`\n",
        "do\n",
        "samp=${bwig%.bw}\n",
            "qsub -N bwtool_${ANALYSIS}_${rgn}_${samp} -l h_vmem=10G,virtual_free=10G -b y \\\n",
            "-e $TMP_DIR/logs/$ANALYSIS/bwtool_${rgn}_${samp}.err -o $TMP_DIR/logs/bwtool_${rgn}_${samp}.log \\\n",
            if(BINNED)
                "/ngs_share/tools/bwtool_patched/bwtool matrix ${UPSTREAM}:${META}:${DOWNSTREAM} \\\n"
            else
                "/ngs_share/tools/bwtool/bwtool matrix ${UPSTREAM}:${DOWNSTREAM} ${ANCHOR} \\\n",
            "$ANNOTATION_DIR/${rgn}.bed $BW_DIR/$bwig \\\n",
            "$OUT_DIR/summary_${rgn}_${TOKEN}_${UPSTREAM}_${DOWNSTREAM}_${TILE_AVG}_${samp}.txt \\\n",
            #if(!BINNED) "-tiled-averages=${TILE_AVG}\n" else "\n",
            "-tiled-averages=${TILE_AVG}\n",
            "done\n",
            "done\n",
     sep="")
    
    system(command)
    #print(command)
    #waitForClusterJobs<-function(analysis_id, lookup_int=10, verbose=TRUE){
    
}
repeat{
    lookup_cmd<-sprintf("qstat -r | grep \"Full jobname\" | grep -e %s_%s", "bwtool", ANALYSIS)
    suppressWarnings({
                running_jobs<-system(lookup_cmd, intern=TRUE)
            })
    if(length(running_jobs)>0){
        cat(sprintf("[SGE jobs:] %d remaining\n", length(running_jobs)))
        Sys.sleep(10)
    }else{
        break;
    }
}
#}

AVERAGES<-TRUE
HEATMAP<-TRUE

AVG_SCATTERPLOTS<-TRUE
AVG_BOXPLOTS<-TRUE
################################################
################ read the extracted matrices
################################################

final_data<-final_heatmap_data<-final_reg_data<-final_swprofiles<-final_rmats<-setNames(as.list(rep(NA, length(features))), features)
region_ids<-feature_grs<-list()
#heatmaps<-list()
for(feature in features){
    print(feature)
    final_data[[feature]]<-final_heatmap_data[[feature]]<-final_reg_data[[feature]]<-final_swprofiles[[feature]]<-final_rmats[[feature]]<-list()
    
    fbed<-read.table(file.path(ANNOTATION_DIR, paste0(feature,".bed")))
    feature_grs[[feature]]<-GRanges(fbed$V1, IRanges(fbed$V2, fbed$V3), "*")#, fbed$V4)
    region_ids[[feature]]<-paste0("region_", seq_along(feature_grs[[feature]]))
    #heatmaps[[feature]]<-list
    
#    if(BINNED){
#        matrix<-read.table(file.path(STATE_RGNS_DIR, sprintf("summary_%s_binned_%d_%d_%d.txt",feature, OFFSET, OFFSET,WIN_SIZE)))
#    }else{
#        #matrix<-read.table(file.path(STATE_RGNS_DIR, sprintf("summary_feature_%d_1500bp.txt",feature)))
#        matrix<-read.table(file.path(STATE_RGNS_DIR, sprintf("summary_%s_center_%d_%d_%d.txt",feature, OFFSET, OFFSET, WIN_SIZE)))
#    }
    
#    nwin<-ncol(matrix)/nsamp
#    
#    indices<-lapply(seq(nsamp)-1, function(si) seq(si*nwin+1,nwin*(si+1)))
#    
#    mat_list<-lapply(indices, function(idx) matrix[,idx])
#    names(mat_list)<-bw_files

    mat_list<-list()
    for(sample in unlist(sample_groups)){
        for(mark in c(marks)){
                mat_file<-file.path(STATE_RGNS_DIR, sprintf(
                                "summary_%s_%s_%d_%d_%d_%s_%s.txt",
                                feature, if(BINNED) "binned" else ANCHOR_SPECS, OFFSET, OFFSET,WIN_SIZE, sample,mark))
                if(file.exists(mat_file)){
                    mat_list[[paste(sample, mark, sep="_")]]<-read.table(mat_file)
                }else{
                    print(paste(mat_file, "does not exist"))
                }
        }
    }
    
    print(lapply(mat_list, dim))
    
    #nwin<-ncol(mat_list[[1]])
    
    nreg<-nrow(mat_list[[1]])
    
    rand_subset<-sort(sample.int(nreg, min(nreg,5000)))
    
    #for(i in seq(length(mat_list))){
    for(mark in c(marks)){
        print(mark)
        
        if(!mark %in% c("WGBS", "RNAseq", "CAGE") & ((mark %in% c("ATAC") & NORMALIZE_ATAC) || NORMALIZE_CHIP)){
            norm_factors<-setNames(seq_depths[[mark]]/min(seq_depths[[mark]], na.rm=TRUE), seq_depths[["Sample"]])
        }else{
            norm_factors<-setNames(rep(1, length(unlist(sample_groups))), unlist(sample_groups))
        }
        sw_profiles<-sg_profiles<-sg_averages<-sg_matrices<-sg_rmats<-setNames(as.list(rep(NA, length(sample_groups))), names(sample_groups))
        sg_vars<-list()
        sg_sems<-list()
        sg_maxs<-sg_mins<-list()
        sg_dfs<-list()
        for(sgn in names(sample_groups)){
            sg<-sample_groups[[sgn]]
            sw_profiles[[sgn]]<-list()
            sg_profiles[[sgn]]<-list()
            sg_rmats[[sgn]]<-list()
            sg_averages[[sgn]]<-list()
            sg_matrices[[sgn]]<-list()
            
            present_samples<-0
            sg_matrix<-matrix(0, ncol=nwin, nrow=nreg)
            for(sample in sg){
                print(sample)
                if(sprintf("%s_%s",sample, mark) %in% names(mat_list)){
                    present_samples<-present_samples+1
                    sg_matrix<-sg_matrix + mat_list[[sprintf("%s_%s",sample, mark)]]/norm_factors[sample]
                    profile<-colMeans(mat_list[[sprintf("%s_%s",sample, mark)]], na.rm=TRUE)
                    rgn_avg<-rowMeans(mat_list[[sprintf("%s_%s",sample, mark)]], na.rm=TRUE)
                    
                    profile<-profile/norm_factors[sample]
                    rgn_avg<-rgn_avg/norm_factors[sample]
                    sg_averages[[sgn]][[sample]]<-rgn_avg
                    sg_profiles[[sgn]][[sample]]<-profile
                    
                    #            plot(bin_locs, profile, 
                    #                    ylab="Read counts", type="l", ylim=c(0,15), main=sprintf("%s_%s",sample, mark))
                }else{
                    #plot(0,0,ylim=c(0,22), main=sprintf("%s_%s",sample, mark))
                    sg_profiles[[sgn]][[sample]]<-rep(NA,nwin)
                    sg_averages[[sgn]][[sample]]<-rep(NA,nreg)
                    
                }
            }
            #browser()
            imat<-do.call("rbind", sg_profiles[[sgn]])
            sw_profiles[[sgn]]<-imat[which(apply(!is.na(imat),1,any)),]
            sg_profiles[[sgn]]<-colMeans(imat, na.rm=TRUE)
            sg_vars[[sgn]]<-apply(imat, 2, sd, na.rm=TRUE)
            sg_sems[[sgn]]<-sg_vars[[sgn]]/sqrt(present_samples)
            
            if(AVERAGES){
                rmat<-do.call("cbind", sg_averages[[sgn]])
                sg_rmats[[sgn]]<-rmat[,which(apply(!is.na(rmat),2,any))]
                sg_averages[[sgn]]<-rowMeans(rmat, na.rm=TRUE)
                #sg_vars[[length(sg_profiles)]]<-apply(imat, 2, sd, na.rm=TRUE)
            }
            
            if(HEATMAP){
                sg_matrix<-sg_matrix/present_samples
                sg_matrices[[sgn]]<-sg_matrix
                
                if(!all(as.numeric(as.matrix(sg_matrix))==0)){
                    #rownames(sg_matrix)<-sprintf("region_%s", seq(nrow(sg_matrix)))
                    rownames(sg_matrix)<-region_ids[[feature]]
                    sg_matrix<-sg_matrix[order(rowSums(sg_matrix), decreasing=TRUE),]
                    sg_df<-cbind(data.frame(Distance_kb=bin_locs/1000), t(sg_matrix))#[rand_subset,]))
                    sg_df<-gather(sg_df, key="Region", value="value", -Distance_kb)
                    sg_df$Group<-rep(sgn, nrow(sg_df))
                    sg_dfs[[length(sg_dfs)+1]]<-sg_df
                }
            }
            
        }
        
        if(HEATMAP){
            if(length(sg_dfs)>0){
                sg_dfs<-do.call("rbind", sg_dfs)
                sg_dfs$Region<-factor(sg_dfs$Region, levels=rev(unique(sg_dfs$Region)))
                sg_dfs$Feature<-rep(feature, nrow(sg_dfs))
                sg_dfs$Mark<-rep(mark, nrow(sg_dfs))
                final_heatmap_data[[feature]][[mark]]<-sg_dfs
            }else{
                print("not found any heatmap data!!!")
            }
        }
        #sw_profiles<-do.call("rbind", sw_profiles)
        final_swprofiles[[feature]][[mark]]<-sw_profiles

        sg_profiles<-do.call("rbind", sg_profiles)
        sg_vars<-do.call("rbind", sg_vars)
        sg_sems<-do.call("rbind", sg_sems)
        rownames(sg_vars)<-paste0("SD_", names(sample_groups))#c("var_WT", "var_G34W")
        rownames(sg_sems)<-paste0("SEM_", names(sample_groups))#c("var_WT", "var_G34W")
        
        rownames(sg_profiles)<-paste0("average_", names(sample_groups))#c("average_WT", "average_G34W")
        final_data[[feature]][[mark]]<-cbind(
                data.frame(
                       #Feature=rep(paste0("Feature_", feature), ncol(sg_profiles)),
                       Feature=rep(feature, ncol(sg_profiles)), 
                Mark=rep(mark, ncol(sg_profiles)),
                Distance=bin_locs),
                as.data.frame(t(sg_profiles)),
                as.data.frame(t(sg_vars)),
                as.data.frame(t(sg_sems))
               )
       if(AVERAGES){
           sg_averages<-do.call("cbind", sg_averages)
           #sg_vars<-do.call("rbind", sg_vars)
           #rownames(sg_vars)<-paste0("var_", names(sample_groups))#c("var_WT", "var_G34W")
           colnames(sg_averages)<-paste0("average_", names(sample_groups))#c("average_WT", "average_G34W")
           final_reg_data[[feature]][[mark]]<-cbind(
                   data.frame(
                           Feature=rep(feature, nrow(sg_averages)), 
                           Mark=rep(mark, nrow(sg_averages)),
                           Region=region_ids[[feature]]),
                   as.data.frame(sg_averages)#,
                   #as.data.frame(t(sg_vars))
           )
           final_rmats[[feature]][[mark]]<-sg_rmats
        }
        #matplot(bin_locs, t(sg_profiles), 
#                ylim=c(0,mark_lims[mark]),
#                ylab="Read counts", type="l", main=mark) #main=sprintf("%s_%s",sample, mark))
    }
}
saveRDS(final_swprofiles, file=file.path(STATE_RGNS_DIR, "final_swprofiles.RDS"))
saveRDS(final_data, file=file.path(STATE_RGNS_DIR, "final_data.RDS"))
saveRDS(final_reg_data, file=file.path(STATE_RGNS_DIR, "final_reg_data.RDS"))
saveRDS(final_rmats, file=file.path(STATE_RGNS_DIR, "final_rmats.RDS"))
saveRDS(final_heatmap_data, file=file.path(STATE_RGNS_DIR, "final_heatmap_data.RDS"))

if(FALSE){
    #ANALYSIS_GROUP<-"H3.3_G34W"
    ANALYSIS_GROUP<-"Bivalent_genes"
    
    STATE_RGNS_DIR2=file.path("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/",
            ANALYSIS_GROUP,
            paste(ANALYSIS, if(BINNED) "binned" else "centered", OFFSET, WIN_SIZE, sep="_"))
    
    final_swprofiles<-readRDS(file=file.path(STATE_RGNS_DIR2, "final_swprofiles.RDS"))
    final_data<-readRDS(file=file.path(STATE_RGNS_DIR2, "final_data.RDS"))
    final_reg_data<-readRDS(file=file.path(STATE_RGNS_DIR2, "final_reg_data.RDS"))
    final_rmats<-readRDS(file=file.path(STATE_RGNS_DIR2, "final_rmats.RDS"))
    final_heatmap_data<-readRDS(file=file.path(STATE_RGNS_DIR2, "final_heatmap_data.RDS"))
}


######################################################################## 
############################################ SIGNIFICANCE
#library(Hotelling)
#library(ICSNP)
#library(FRB)
library(EmpiricalBrownsMethod)
library(Hmisc)

calculateKostCovariance2 <- function(data_matrix) {
    covar_matrix = EmpiricalBrownsMethod:::kostPolyFit(Hmisc:::rcorr(t(data_matrix))$r)
    return(covar_matrix)
}
kost_method<-function(data_matrix, p_values, extra_info=FALSE){
    covar_matrix<-calculateKostCovariance2(data_matrix)
    EmpiricalBrownsMethod:::combinePValues(covar_matrix, p_values, extra_info = extra_info)
}

########################################################################################
#### using averages 
TEST_METHOD=c("averages", "profiles")[1]

comp_groups<-c("WT", "MT")
htsq<-list()
t.test.pvals<-list()
for(feat in features){
    for(mark in marks){
        if(TEST_METHOD=="averages"){
            grp1<-t(final_rmats[[feat]][[mark]][[comp_groups[1]]])
            grp2<-t(final_rmats[[feat]][[mark]][[comp_groups[2]]])
        }else if(TEST_METHOD=="profiles"){
            grp1<-final_swprofiles[[feat]][[mark]][[comp_groups[1]]]
            grp2<-final_swprofiles[[feat]][[mark]][[comp_groups[2]]]
        }
        if(!is.null(nrow(grp1)) && !is.null(nrow(grp2)) && (nrow(grp1)>1) && (nrow(grp2)>1)){
            nna1<-which(colSums(is.na(grp1) | grp1==0)==0)
            nna2<-which(colSums(is.na(grp2) | grp2==0)==0)
            nna<-intersect(nna1,nna2)
            nzsd1<-which(apply(grp1, 2, sd)>0)
            nzsd2<-which(apply(grp2, 2, sd)>0)
            grp1<-grp1[,nna]
            grp2<-grp2[,nna]
            
            print(length(nna))
            
            #ht<-hotelling.test(grp1[,nna], grp2[,nna],shrinkage=TRUE)
            #ht<-HotellingsT2(grp1[,nna], grp2[,nna])
            #ht<-FRBhotellingMM(grp1[,nna], grp2[,nna])
            
            #htsq[[length(htsq)+1]]<-ht
            
            p_vals<-sapply(seq(nna), function(i) t.test(grp1[,i], grp2[,i])$p.value)
            t.test.pvals[[sprintf("%s_%s", feat, mark)]]<-p_vals
            if(!is.null(nrow(grp1)) && !is.null(nrow(grp2)) && (nrow(grp1)+nrow(grp2)>4)){
                kmp<-kost_method(t(rbind(grp1,grp2)), p_vals, extra_info=TRUE)
                htsq[[sprintf("%s_%s", feat, mark)]]<-data.frame(Feature=feat, Mark=mark, KostP=kmp$P_test)
            }else{
                p_sum <- 2 * sum(-log(p_vals))
                df_fisher <- 2 * (ncol(grp1))
                p_fisher <- pchisq(df = df_fisher, q = p_sum, lower.tail = FALSE)
                htsq[[sprintf("%s_%s", feat, mark)]]<-data.frame(Feature=feat, Mark=mark, KostP=p_fisher)
            }
        }else{
            htsq[[sprintf("%s_%s", feat, mark)]]<-data.frame(Feature=feat, Mark=mark, KostP=NA)
        }
    }
}
htsq<-do.call("rbind", htsq)
htsq$p.value<-htsq$KostP

#plotting

######################################################################## 
########################   PROFILE PLOTS

final_df<-do.call("rbind",unlist(final_data, recursive=FALSE))
rownames(final_df)<-NULL

#final_df_avg_wt<-final_df[,-grep(c("_G34W"),colnames(final_df))]
#final_df_avg_mt<-final_df[,-grep(c("_WT"),colnames(final_df))]
#colnames(final_df_avg_wt)[grep("WT", colnames(final_df_avg_wt))]<-colnames(final_df_avg_mt)[grep("G34W", colnames(final_df_avg_mt))]<-c("Mean", "Variance")
#final_df_avg_mt$H3F3A<-rep("G34W", ncol(final_df_avg_mt))
#final_df_avg_wt$H3F3A<-rep("WT", ncol(final_df_avg_wt))
#final_df_clean<-rbind(final_df_avg_wt, final_df_avg_mt)

lll<-list()
for(gn in names(sample_groups)){
    
    final_df_avg_grp<-final_df[,-grep(paste(paste0("_",setdiff(names(sample_groups),gn), collapse="|")),colnames(final_df))]
    colnames(final_df_avg_grp)[grep(gn, colnames(final_df_avg_grp))]<-c("Mean", "SD", "SEM")
    final_df_avg_grp$Group<-rep(gn, nrow(final_df_avg_grp))
    
    if(SAMPLE_WISE){
        final_df_avg_grp$SampleWiseGroup<-final_df_avg_grp$Group
        final_df_avg_grp$Group<-substr(names(unlist(mt_wt_sample_groups))[match(final_df_avg_grp$SampleWiseGroup, unlist(mt_wt_sample_groups))], 1,2)
    }
    
    lll[[gn]]<-final_df_avg_grp
}

final_df_clean<-do.call("rbind", lll)

#final_df_clean$Mean[is.na(final_df_clean$Mean)]<-0

final_df_clean$Distance_kb<-final_df_clean$Distance/1000

if(BINNED){
    final_df_clean$Distance_kb<-factor(final_df_clean$Distance_kb, levels=bin_locs/1000)
}

if(any(grepl("RNA", features))){
    
    final_df_clean$Feature<-gsub("RNAseq_", "", final_df_clean$Feature)
    final_df_clean$Feature<-gsub("tss", "TSS", final_df_clean$Feature)
    final_df_clean$Feature<-gsub("tes", "TES", final_df_clean$Feature)
    final_df_clean$Feature<-gsub("_expressed", "", final_df_clean$Feature)
    
    levs<-gsub("RNAseq_", "", features)
    levs<-gsub("tss", "TSS", levs)
    levs<-gsub("tes", "TES", levs)
    levs<-gsub("_expressed", "", levs)
    final_df_clean$Feature<-factor(final_df_clean$Feature, levs)
    
}else if(any(grepl("bivalent_ESC_cluster_combined_category", features))){
    final_df_clean$Feature<-as.character(final_df_clean$Feature)
    #final_df_clean$Feature<-gsub("Feature_", "", final_df_clean$Feature)
    final_df_clean$Feature<-gsub("bivalent_ESC_cluster_combined_category_", "", final_df_clean$Feature)
    final_df_clean$Feature<-gsub("bivalent_ESC_", "", final_df_clean$Feature)
    final_df_clean$Feature<-factor(final_df_clean$Feature,sort(gsub("bivalent_ESC_cluster_combined_category_|bivalent_ESC_","",features)))
    
}else if(any(grepl("TssAFlnk", features))){
    final_df_clean$Feature<-gsub("Feature_", "", final_df_clean$Feature)
    chmm_annot<-read.table("/ngs_share/scratch/pavlo/gctb/data/annotations/chromhmm/ENCODE/chmm_15_state_model_states.txt", sep="\t", header=TRUE)
    final_df_clean$Feature<-factor(final_df_clean$Feature,gsub("/", "_", chmm_annot$MNEMONIC))
#}else if(any(grepl("H3.3_", features))){
    #final_df_clean$Feature<-gsub("Feature_", "", final_df_clean$Feature)
#    final_df_clean$Feature<-factor(final_df_clean$Feature, levels=c("H3.3_WT_unique","H3.3_MT_unique","G34W_unique"))
}
#}else if(any(grepl("G34W_filt", features))){
#    final_df_clean$Feature<-gsub("Feature_", "", final_df_clean$Feature)
#    lvls<-c("H3.3_WT","H3.3_WT_unique","H3.3_MT","H3.3_MT_unique","H3.3_MT_rep2","H3.3_overlap","G34W","G34W_filt","G34W_all","G34W_unique","G34W_MT_olap", "H3.Y_depleted", "h3.3_background")
#    final_df_clean$§e<-factor(final_df_clean$Feature,lvls)
#}else{
#    final_df_clean$Feature<-gsub("Feature_", "", final_df_clean$Feature)
#    final_df_clean$Feature<-factor(final_df_clean$Feature,features)
#}

#library(tidyr)
#final_df_melt<-gather(final_df_clean,key="Key", value="Value", -Distance, -State, -Mark)
#
#### subset the data
#final_df_clean_subs<-final_df_clean[final_df_clean$Feature %in% c("TssBiv", "BivFlnk", "EnhBiv") & final_df_clean$Mark %in% c("WGBS", "ATAC", "K4me3", "K27me3"),]
#final_df_clean_subs<-final_df_clean[final_df_clean$Feature %in% "H3.Y_depleted",]
#final_df_clean_subs2<-final_df_clean[final_df_clean$Mark %in% "K27me3",]


feat_column<-unlist(lapply(names(final_swprofiles), function(feat) rep(feat,length(final_swprofiles[[feat]]))))
mk_column<-unlist(lapply(names(final_swprofiles), function(feat) lapply(names(final_swprofiles[[feat]]), function(mk) rep(mk,length(bin_locs)))))
final_swprofiles_df<-do.call("rbind", setNames(unlist(unlist(final_swprofiles, recursive=FALSE), recursive=FALSE), NULL))

final_swprofiles_df<-list()
for(feat in names(final_swprofiles)){
    for(mk in marks){
        for(sgn in names(sample_groups)){
            if(is.matrix(final_swprofiles[[feat]][[mk]][[sgn]])){
                final_swprofiles_df[[length(final_swprofiles_df)+1]]<-data.frame(final_swprofiles[[feat]][[mk]][[sgn]], Feature=feat, Mark=mk, Group=sgn)
            }else{
#                browser()
                final_swprofiles_df[[length(final_swprofiles_df)+1]]<-data.frame(as.data.frame(t(final_swprofiles[[feat]][[mk]][[sgn]])), Feature=feat, Mark=mk, Group=sgn)
            }
        }
    }
}
final_swprofiles_df<-do.call("rbind", final_swprofiles_df)
final_swprofiles_df$Sample<-rownames(final_swprofiles_df)
final_swprofiles_df<-gather(final_swprofiles_df, key="Distance", value="Mean", -Feature, -Mark, -Group,-Sample)

final_swprofiles_df$Distance<-setNames(bin_locs, paste0("V", seq(bin_locs)))[final_swprofiles_df$Distance]
final_swprofiles_df$Distance_kb<-final_swprofiles_df$Distance/1000

profile_meta_plot<-function(summarized_df, facet_vars=c("Mark", "Feature"), free_scales=FALSE, col_var="Feature", pval_df=NULL, per_sample=FALSE, binned=FALSE){
    
    library(ggplot2)
    library(ggpubr)
    library(ggsci)
    
    if(per_sample){
        p<-ggplot(summarized_df, aes(x=Distance_kb, y=Mean, color=Group, group=Sample)) 
    }else{
        p<-ggplot(summarized_df, aes(x=Distance_kb, y=Mean, color=Group))
    }
    
    if(BINNED){
        p<-p + geom_line(lwd=1,aes(group=Group))
    }else{
        p<-p + geom_line(lwd=1)
    }
    
    if(!per_sample & !binned){
        p<-p + geom_ribbon(aes(ymin=Mean-SEM, ymax=Mean+SEM, fill=Group),alpha=0.3, lwd=0.25)
    }
    
    p<-p + scale_colour_manual(values=group_colors)
    
    if(binned){
        p<-p + scale_x_discrete(
            #breaks=c(floor(n_flank_win/2), n_flank_win, n_flank_win+N_BINS, N_BINS+n_flank_win+ceiling(n_flank_win/2)),
            breaks=c(levels(summarized_df$Distance_kb)[
                            c(ceiling(n_flank_win/2),n_flank_win+1, N_BINS+n_flank_win, N_BINS+n_flank_win+ceiling(n_flank_win/2))]),
            labels=c(
                    sprintf("%.1f",round(as.numeric(levels(summarized_df$Distance_kb)[ceiling(n_flank_win/2)]))),
                    BINNED_REG_START_NAME, BINNED_REG_END_NAME, 
                    sprintf("%.1f",round(as.numeric(levels(summarized_df$Distance_kb)[N_BINS+n_flank_win+ceiling(n_flank_win/2)])))
            ), 
            expand = c(0, 0), position="bottom")
    }
    if(binned){
        xl<-"distance to anchor, kb"
    }else{
        xl<-"distance to center, kb"
    }
    p<-p + xlab(xl) + ylab("normalized read count")
    if(free_scales){
        #p<-p + facet_wrap(Feature~Mark, scales="free", ncol=length(unique(final_df_clean$Mark)))
         p<-p + facet_wrap(as.formula(sprintf("%s~%s",facet_vars[1],facet_vars[2])), scales="free", 
                 ncol=length(unique(summarized_df[[col_var]])), drop=TRUE)
    }else{
        p<-p + facet_grid(as.formula(sprintf("%s~%s",facet_vars[1],facet_vars[2])), scales="free", drop=TRUE)
    }
    #p<-p + facet_grid(Feature~Mark, scales="free")
    
    if(!is.null(pval_df)){
        pval_df$p_val_formatted=sprintf("P=%1.2g", pval_df$p.value)
        pval_df$Group<-rep(NA, nrow(pval_df))
        p<-p + geom_text(mapping=aes(x=-Inf, y=Inf, label=p_val_formatted), data=pval_df#, 
                #nudge_x=-1, nudge_y=1,
                ,color="black"
                ,hjust=-0.25
                ,vjust=1
        )
    }

    #p<-p + theme_bw(
    #        axis.text.x=element_text(color="black"),
    #        axis.text.y=element_text(color="black"),
    #        )
    p<-p + theme_classic()
    
    return(p)

}
final_df_clean_subs<-final_df_clean
#if(length(htsq_subs)>0) htsq_subs<-htsq else htsq_subs<-NULL
htsq_subs<-htsq
htsq_subs<-NULL


#profiles_mark_subset<-c("K4me3", "K27me3","WGBS")
profiles_mark_subset<-marks



final_df_clean_subs<-final_df_clean[final_df_clean$Mark %in% profiles_mark_subset,]
if(!is.null(htsq_subs)) htsq_subs<-htsq[ htsq$Mark %in% profiles_mark_subset,]
final_df_clean_subs$Mark<-factor(final_df_clean_subs$Mark, levels=profiles_mark_subset)
if(!is.null(htsq_subs)) htsq_subs$Mark<-factor(htsq_subs$Mark, levels=profiles_mark_subset)


pdf(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/profile_summaries_SEM_%s_pval_subs.pdf", 
                paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")),
        width=1*length(unique(final_df_clean_subs$Feature))+3.5, height=1*length(unique(final_df_clean_subs$Mark))+2)# width=20, height=7)
#print(p)
print(profile_meta_plot(final_df_clean_subs, pval_df=htsq_subs, binned=BINNED))
dev.off()


final_swprofiles_df_subs<-final_swprofiles_df[final_swprofiles_df$Mark %in% profiles_mark_subset,]
final_swprofiles_df_subs$Mark<-factor(final_swprofiles_df_subs$Mark, levels=profiles_mark_subset)


pdf(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/profile_summaries_SEM_%s_pval_subs_samplewise.pdf", 
                paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")),
        width=1*length(unique(final_swprofiles_df_subs$Feature))+3.5, height=1*length(unique(final_swprofiles_df_subs$Mark))+2)# width=20, height=7)
#print(p)
print(profile_meta_plot(final_swprofiles_df_subs, pval_df=NULL, per_sample=TRUE))
dev.off()


########

comp_groups<-c("WT","MT")
#final_df_clean_subs$Group<-factor(final_df_clean, levels=c("WT", "MT"))
final_df_clean_subs_diff<-final_df_clean_subs[,c("Feature","Mark","Distance","Mean","Distance_kb", "Group")]
#final_df_clean_subs_diff<-final_df_clean_subs
final_df_clean_subs_diff<-spread(final_df_clean_subs_diff, key="Group", value="Mean")
final_df_clean_subs_diff$Delta<-final_df_clean_subs_diff[[comp_groups[2]]]-final_df_clean_subs_diff[[comp_groups[2]]]

################################################

if(AVG_SCATTERPLOTS || AVG_BOXPLOTS || SEGMENT_CHANGE_PLOTS){
    
    final_reg_df<-do.call("rbind",unlist(final_reg_data, recursive=FALSE))
    rownames(final_reg_df)<-NULL
    
    if(any(grepl("RNA", features))){
        
        final_reg_df$Feature<-gsub("Feature_RNAseq_", "", final_reg_df$Feature)
        final_reg_df$Feature<-gsub("tss", "TSS", final_reg_df$Feature)
        final_reg_df$Feature<-gsub("tes", "TES", final_reg_df$Feature)
        final_reg_df$Feature<-gsub("_expressed", "", final_reg_df$Feature)
        
        levs<-gsub("RNAseq_", "", features)
        levs<-gsub("tss", "TSS", levs)
        levs<-gsub("tes", "TES", levs)
        levs<-gsub("_expressed", "", levs)
        final_reg_df$Feature<-factor(final_reg_df$Feature, levs)
    }else if(any(grepl("bivalent_ESC_cluster_combined_category", features))){
        final_reg_df$Feature<-as.character(final_reg_df$Feature)
        final_reg_df$Feature<-gsub("Feature_", "", final_reg_df$Feature)
        final_reg_df$Feature<-gsub("bivalent_ESC_cluster_combined_category_", "", final_reg_df$Feature)
        final_reg_df$Feature<-gsub("bivalent_ESC_cluster_", "Cluster_", final_reg_df$Feature)
        feature_levels<-gsub("bivalent_ESC_cluster_combined_category_","",features)
        feature_levels<-gsub("bivalent_ESC_cluster_","Cluster_",feature_levels)
        
        final_reg_df$Feature<-factor(final_reg_df$Feature,sort(feature_levels))
        
    }else if(any(grepl("TssAFlnk", features))){
        final_reg_df$Feature<-gsub("Feature_", "", final_reg_df$Feature)
        chmm_annot<-read.table("/ngs_share/scratch/pavlo/gctb/data/annotations/chromhmm/chmm_15_state_model_states.txt", sep="\t", header=TRUE)
        final_reg_df$Feature<-factor(final_reg_df$Feature,gsub("/", "_", chmm_annot$MNEMONIC))
    }else if(any(grepl("G34W_filt", features))){
        final_reg_df$Feature<-gsub("Feature_", "", final_reg_df$Feature)
        lvls<-c("H3.3_WT","H3.3_WT_unique","H3.3_MT","H3.3_MT_rep2", "H3.3_MT_unique","H3.3_overlap","G34W","G34W_filt","G34W_all","G34W_unique","G34W_MT_olap")
        final_reg_df$Feature<-factor(final_reg_df$Feature,lvls)
    }else{
        final_reg_df$Feature<-gsub("Feature_", "", final_reg_df$Feature)
        final_reg_df$Feature<-factor(final_reg_df$Feature,features)
    }
    
}
######################################################################## 
########################   AVERAGE SCATTERPLOTS
if(AVG_SCATTERPLOTS){
    
    region_average_scatterplot<-function(final_reg_df, comp_groups=NULL, comp_marks=NULL, facet_vars=NULL, scatter_fit=TRUE){
       
        
        if(is.null(comp_groups) & is.null(comp_marks)){
            stop("either of comp_groups or comp_marks has to be specified")
        }
        
        CoordSymmetric <- ggproto("CoordSymmetric", CoordCartesian,
                
                range = function(scale_details) {
                    print("adjusting_ranges")
                    cur_x_range<-scale_details$x.range
                    cur_y_range<-scale_details$y.range
                    sym_range<-c(min(cur_x_range[1L], cur_y_range[1L]), max(cur_x_range[2L], cur_y_range[2L]))
                    return(list(x = sym_range, y = sym_range))
                },
                train = function (self, scale_details)
                {
                    train_cartesian <- function(scale_details, limits, name) {
                        if (self$expand) {
                            expand <- ggplot2:::expand_default(scale_details)
                        }
                        else {
                            expand <- c(0, 0)
                        }
                        if (is.null(limits)) {
                            range <- scale_details$dimension(expand)
                        }
                        else {
                            range <- range(scale_details$transform(limits))
                            range <- scales:::expand_range(range, expand[1], expand[2])
                        }
                        out <- scale_details$break_info(range)
                        out$arrange <- scale_details$axis_order()
                        names(out) <- paste(name, names(out), sep = ".")
                        out
                    }
                    trained_ranges<-c(
                            train_cartesian(scale_details$x, self$limits$x, "x"), 
                            train_cartesian(scale_details$y, self$limits$y, "y"))
                    cur_x_range<-trained_ranges$x.range
                    cur_y_range<-trained_ranges$y.range
                    sym_range<-c(min(cur_x_range[1L], cur_y_range[1L]), max(cur_x_range[2L], cur_y_range[2L]))
                    c(
                            train_cartesian(scale_details$x, sym_range, "x"), 
                            train_cartesian(scale_details$y, sym_range, "y"))
                    
                }
        )
        coord_symmetric<-function(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE){
            ggproto(NULL, CoordSymmetric,
                    limits = list(x = xlim, y = ylim),
                    ratio = ratio,
                    expand = expand
            )
        }
        
        if(!is.null(comp_marks)){
            
            lll<-list()
            for(gn in names(sample_groups)){
                
                final_df_avg_grp<-final_reg_df[,-grep(paste(paste0("_",setdiff(names(sample_groups),gn), collapse="|")),colnames(final_reg_df))]
                colnames(final_df_avg_grp)[grep(gn, colnames(final_df_avg_grp))]<-c("Mean")
                final_df_avg_grp$Group<-rep(gn, nrow(final_df_avg_grp))
                
                lll[[gn]]<-final_df_avg_grp
            }
            final_reg_df_clean<-do.call("rbind", lll)
            
        }
        if(!is.null(facet_vars)){
         fvl<-Reduce("*", apply(final_reg_df[facet_vars],2,function(col) length(unique(col))))
        }else{
         fvl<-0
        }
        if(is.null(comp_marks)){
            
            #comp_groups<-c("MSC", "GCTB_WT", "MT")[c(2,3)]
            #comp_groups<-c("WT", "MT")
            
            p<-ggplot(final_reg_df, aes_string(x=sprintf("average_%s", comp_groups[1]), y=sprintf("average_%s", comp_groups[2]))) + geom_point(size=0.5, alpha=0.5)
            
            p<-p+geom_abline(slope=1,intercept=0,linetype="dashed")
            p<-p + coord_symmetric()#coord_fixed(ratio=1)# coord_cartesian(ylim = c(-20, 80)) 
            
            p<-p + facet_wrap(c("Mark","Feature", facet_vars), 
                    ncol=length(unique(final_reg_df$Feature))*(fvl+(fvl==0)), scales="free")
            
            p<-p + theme_bw()
            
        }
            ################ different marks vs each other
        else if(is.null(comp_groups)){
            
            final_reg_df_by_mark<-spread(final_reg_df_clean, key=Mark, value=Mean)
            
            ###
        #    comp_marks<-c("K27me3", "K4me3")
        #    comp_marks<-c("ATAC", "WGBS")
        #    comp_marks<-c("ATAC", "G34W")
        #    comp_marks<-c("K27me3", "G34W")
            
            
            PNG<-FALSE
            
            p<-ggplot(final_reg_df_by_mark, aes_string(x=comp_marks[1], y=comp_marks[2])) + geom_point(size=0.5, alpha=0.5)
            
            p<-p+geom_abline(slope=1,intercept=0,linetype="dashed")
            
            #p<-p + coord_symmetric()#coord_fixed(ratio=1)# coord_cartesian(ylim = c(-20, 80)) 
            
            p<-p + facet_wrap(c("Group","Feature",facet_vars), 
                    ncol=length(unique(final_reg_df$Feature))*(fvl+(fvl==0)+is.null(fvl)))
            
            p<-p + theme_bw()
            
        }
        ######### changes vs each other
        else if(!is.null(comp_groups) & !is.null(comp_marks)){
#            comp_marks<-c("G34W", "K27me3")
#            comp_groups<-c("WT", "MT")
            
            PNG<-FALSE
            
            final_reg_diffs<-final_reg_df_clean %>% unite(Mark_Group, Mark, Group) %>% spread(Mark_Group, Mean)
            
            final_reg_diffs$Diff_mark1<-final_reg_diffs[[sprintf("%s_%s", comp_marks[1], comp_groups[2])]]-final_reg_diffs[[sprintf("%s_%s", comp_marks[1], comp_groups[1])]]
            final_reg_diffs$Diff_mark2<-final_reg_diffs[[sprintf("%s_%s", comp_marks[2], comp_groups[2])]]-final_reg_diffs[[sprintf("%s_%s", comp_marks[2], comp_groups[1])]]
            
            ####
            p<-ggplot(final_reg_diffs, aes(x=Diff_mark1, y=Diff_mark2)) + geom_point(size=0.5, alpha=0.5)
            p<-p+geom_vline(xintercept=0,linetype="dashed")+geom_hline(yintercept=0,linetype="dashed")
            
            p<-p + xlab(sprintf("%s, %s - %s",  comp_marks[1], comp_groups[2], comp_groups[1]))
            p<-p + ylab(sprintf("%s, %s - %s",  comp_marks[2], comp_groups[2], comp_groups[1]))
            
            #p<-p + coord_symmetric()#coord_fixed(ratio=1)# coord_cartesian(ylim = c(-20, 80)) 
            
            p<-p + facet_wrap(c("Feature", facet_vars), ncol=length(unique(final_reg_df$Feature)))
            
            p<-p + theme_bw()
            
            if(scatter_fit){
                
                lm_eqn <- function(df, y, x){
                    formula = as.formula(sprintf('%s ~ %s', y, x))
                    m <- lm(formula, data=df);
                    # formating the values into a summary string to print out
                    # ~ give some space, but equal size and comma need to be quoted
                    eq <- substitute(
                            #italic(target) == a + b %.% italic(input)*","~~italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                            italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
                            list(target = "y",
                                    input = "x",
                                    a = format(as.vector(coef(m)[1]), digits = 2), 
                                    b = format(as.vector(coef(m)[2]), digits = 2), 
                                    r2 = format(summary(m)$r.squared, digits = 3),
                                    # getting the pvalue is painful
                                    pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
                            )
                    )
                    as.character(as.expression(eq));
                }
                
                eq<-plyr:::ddply(final_reg_diffs,c("Feature", facet_vars), lm_eqn, y='Diff_mark1', x='Diff_mark2')
                p<-p+geom_smooth(method = "lm", se = TRUE, color='red', size=1, linetype=2)
                p<-p+geom_text(data=eq, 
                        aes(#x=0.75*max(final_reg_diffs$Diff_mark1),y=0.9*max(final_reg_diffs$Diff_mark2),
                                x=-Inf, y=Inf,
                                label=V1), 
                        vjust=1, hjust=0, color='black', parse=T)
            }
            #
        }
        return(p)
    }
    
    comp_groups<-c("WT", "MT")
    comp_marks<-c("H3.3", "K27me3")
    
    
    PNG<-FALSE
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/region_scatterplots_%s_%s_vs_%s", 
            paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"), comp_groups[1], comp_groups[2])
    if(PNG){
        png(paste0(fn, ".png"),
                width=250*length(unique(final_reg_df$Feature))+200, height=250*length(unique(final_reg_df$Mark)))# width=20, height=7)
    }else{
        pdf(paste0(fn, ".pdf"),
                width=2.5*length(unique(final_reg_df$Feature))+2, height=2.5*length(unique(final_reg_df$Mark)))# width=20, height=7))
    }
    print(region_average_scatterplot(final_reg_df, comp_groups=c("WT", "MT")))
    dev.off()
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/region_scatterplots_%s_vs_%s_%s", comp_marks[1], comp_marks[2],
            paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))
    if(PNG){
        png(paste0(fn, ".png"),
                width=200*length(unique(final_reg_df$Feature))+200, height=200*length(unique(final_reg_df$Group)))# width=20, height=7)
    }else{
        pdf(paste0(fn, ".pdf"), width=1.5*length(unique(final_reg_df$Feature))+2, height=1.5*length(grep("average_", colnames(final_reg_df))))
    }
    print(region_average_scatterplot(final_reg_df, comp_marks=comp_marks))
    dev.off()
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/region_scatterplot_diffs_%s_vs_%s_in_%s_vs_%s_%s", 
            comp_marks[1], comp_marks[2],comp_groups[1], comp_groups[2],
            paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))
    if(PNG){
        png(paste0(fn, ".png"),
                width=200*length(unique(final_reg_df$Feature))+200, height=200)# width=20, height=7)
    }else{
        pdf(paste0(fn, ".pdf"), width=2*length(unique(final_reg_df$Feature))+2, height=2)
    }
    print(region_average_scatterplot(final_reg_df, comp_groups=c("WT", "MT"), comp_marks=comp_marks))
    dev.off()
    
    
}

####################### AVERAGE BOXPLOTS
#AVG_DELTA_REF_FEATURE="5k_random_not_MSC_biv"
AVG_DELTA_REF_FEATURE=NULL

PNG<-FALSE

if(AVG_BOXPLOTS){
    
    region_average_boxplot<-function(final_reg_df, facet_vars=NULL, comp_groups=NULL, ref_feature=NULL, test="wilcox.test", pval_df=NULL){
        
        if(!is.null(facet_vars)){
            fvl<-Reduce("*", apply(final_reg_df[facet_vars],2,function(col) length(unique(col))))
        }else{
            fvl<-0
        }
            if(is.null(comp_groups)){
                
                final_reg_df_bp<-gather(final_reg_df, key="Group", value="Coverage", -Feature,-Mark,-Region)
                final_reg_df_bp$Group<-gsub("average_", "",final_reg_df_bp$Group)
                final_reg_df_bp$Group<-factor(final_reg_df_bp$Group, levels=names(sample_groups))
                #comp_groups<-c("MSC", "GCTB_WT", "MT")[c(2,3)]
                #comp_groups<-c("WT", "MT")
                plot_var_x<-"Group"
                plot_var_y<-"Coverage"
            }else{
                final_reg_df_bp<-final_reg_df
                final_reg_df_bp$Difference<-final_reg_df_bp[[sprintf("average_%s", comp_groups[2])]]-final_reg_df_bp[[sprintf("average_%s", comp_groups[1])]]
                
                if(is.null(pval_df)){
                        
                    if(is.null(ref_feature)){
                        if(test=="t.test"){
                            pvals <- final_reg_df_bp %>% group_by(Mark, Feature) %>%
                                    #filter(n() > 1) %>%
                                    summarize(p.value = (t.test(Difference, mu = 0)$p.value))
                        }else if(test=="wilcox.test"){
                            pvals <- final_reg_df_bp %>% group_by(Mark, Feature) %>%
                                    #filter(n() > 1) %>%
                                    summarize(p.value = (wilcox.test(Difference, mu = 0)$p.value))
                        }
                    }else{
                        
                        pvals<-list()
                        
                        for(mark in unique(final_reg_df_bp[["Mark"]])){
                            for(feat in unique(final_reg_df_bp[["Feature"]])){
                                grp1<-final_reg_df_bp$Difference[final_reg_df_bp$Mark==mark & final_reg_df_bp$Feature==feat]
                                grp2<-final_reg_df_bp$Difference[final_reg_df_bp$Mark==mark & final_reg_df_bp$Feature==ref_feature]
                                
                                if(test=="t.test"){
                                    pval<-t.test(grp1, grp2)$p.value
                                }else if(test=="wilcox.test"){
                                    pval<-wilcox.test(grp1, grp2)$p.value
                                }
                                pvals[[length(pvals)+1]]<-data.frame(Mark=mark, Feature=feat, p.value=pval)
                                }
                            }
                            
                      
                    }
                    pvals<-do.call("rbind", pvals)
                }else{
                    pvals<-pval_df
                }
                plot_var_x<-sprintf("'%s-%s'", comp_groups[2], comp_groups[1])
                plot_var_y<-"Difference"
                
                
                pvals$p_val_formatted=sprintf("P=%1.2g", pvals$p.value)
                
            }
        
        
        #browser()
        
        p<-ggplot(final_reg_df_bp, 
                aes_string(x=plot_var_x, y=plot_var_y)) + geom_boxplot()
        
        #p<-p+geom_abline(slope=1,intercept=0,linetype="dashed")
        
        if(!is.null(comp_groups)){
            p<-p + geom_hline(yintercept=0, linetype=2)
        }
        if(is.null(comp_groups)) {
            p<-p + facet_wrap(c("Mark","Feature", facet_vars), 
                ncol=length(unique(final_reg_df_bp$Feature))*(fvl+(fvl==0)), scales="free")
            p<-p + xlab("Group means")
        }else{
            p<-p + facet_grid(as.formula("Mark~Feature"),
                   scales="free_y")
            p<-p + xlab("Delta of means")
            p<-p + geom_text(mapping=aes(x=-Inf, y=0, label=p_val_formatted), data=pvals#, 
                    #nudge_x=-1, nudge_y=1,
                    ,hjust=-0.5
                    ,vjust=-1
            )
        }
        p<-p + theme_bw()
        
        p
    }
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/region_boxplot_diffs_%s", 
            paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))
    if(PNG){
        png(paste0(fn, ".png"),
                width=250*length(unique(final_reg_df$Feature))+200, height=250*length(unique(final_reg_df$Mark)))# width=20, height=7)
    }else{
        pdf(paste0(fn, ".pdf"),
                width=2.5*length(unique(final_reg_df$Feature))+2, height=2.5*length(unique(final_reg_df$Mark)))# width=20, height=7))
    }
    print(region_average_boxplot(final_reg_df))
    dev.off()
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/region_boxplot_diffs_%s_vs_%s_%s",
            comp_groups[2], comp_groups[1],
            paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))
    if(PNG){
        png(paste0(fn, ".png"),
                width=150*length(unique(final_reg_df$Feature))+100, height=150*length(unique(final_reg_df$Mark)))# width=20, height=7)
    }else{
        pdf(paste0(fn, ".pdf"),
                width=1.5*length(unique(final_reg_df$Feature))+1, height=1.5*length(unique(final_reg_df$Mark)))# width=20, height=7))
    }
    print(region_average_boxplot(final_reg_df, comp_groups=comp_groups, ref_feature=AVG_DELTA_REF_FEATURE, pval_df=htsq))
    dev.off()
    
}


######################################################################## 
########################   HEATMAPS
#
#pdf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/encode_bivalent_ESC_cluster_2__2_5kb_combined_categories.pdf", width=7, height=10)
#for(feat in features){
#    print(feature)
#    for(mark in marks){
#        
#        print(mark)
#        sg_dfs<-final_heatmap_data[[feat]][[mark]]
#        
#        sg_dfs$value_log2<-log2(sg_dfs$value)
#        sg_dfs$value_log2[sg_dfs$value_log2<0]<-0
#        
#        
#        p<-ggplot(sg_dfs, aes(x=Distance_kb, y=Region, fill=value_log2))
#        p<-p+geom_raster()
#        p<-p+scale_fill_gradient(
#                #colours=rev(colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(length(breaksList))),
#                low="white", high="darkblue")
#        p<-p+theme(axis.text.y=element_blank())
#        p<-p+facet_wrap("Group")#, scales="free")
#        p<-p+ggtitle(sprintf("%s,%s", gsub("bivalent_ESC_cluster_combined_category_","",feat), mark))
#        #pdf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/summaries/encode_bivalent_ESC_cluster_2__2_5kb_heatmap_test_biv_log2_K27me3.pdf", width=5, height=10)
#        print(p)
#    #dev.off()
#    }
#}
#dev.off()
#

######################################################## selected marks
K_MEANS_N_CLUSTERS<-7
MIN_CLUSTER_SIZE=20

#plot_group<-"G34W"
plot_marks<-"G34W"
#plot_marks<-"K27me3"
#plot_marks<-"K36me3"
#plot_marks<-c("G34W", "H3.3")
#plot_marks<-c("G34W","H3.3","K36me3", "K27me3", "K27ac")
#plot_marks<-c("G34W","H3.3","K36me3", "K27me3", "K27ac", "K4me3", "K4me1")
#plot_marks<-c("G34W","H3.3","K36me3", "K27me3", "K27ac", "K9me3", "K4me3", "K4me1")

#plot_marks<-c("G34W","H3.3","K36me3", "K27me3", "K27ac", "K4me3", "K4me1", "WGBS", "ATAC")

#plot_marks<- c("K27me3", "K36me3", "G34W")
#plot_marks<- c("K27me3","WGBS", "ATAC")
##plot_marks<- c("K27me3", "K4me3", "K36me3")
#plot_marks<- c("K27me3", "K4me3", "K27ac", "K36me3", "H3.3", "G34W")
#plot_marks<- c("K27me3", "K4me3", "K36me3", "G34W")

#FEATURE_SUBSET<-c(1:length(features))

#FEATURE_SUBSET<-2

clusters<-list()
#cluster_plots<-list()
CLUSTER_PROFILE_PLOTS<-TRUE
CLUSTER_AVERAGE_PLOTS<-FALSE

CLUSTER_PROFILE_FREE_SCALE<-TRUE

#GROUP_AVG_COMPARISONS<-list(c("WT","MT"))
#MARK_AVG_COMPARISONS<-list(c("G34W","K27me3"), c("G34W", "K36me3"), c("K27me3","K36me3"), c("H3.3", "K36me3"), c("G34W","K27ac"), c("K27me3", "K27ac"))


GROUP_AVG_COMPARISONS<-list()
MARK_AVG_COMPARISONS<-list()


PANEL_SPACING1<-0.25#1.25
PANEL_SPACING2<-1


for(FEATURE_SUBSET in seq(features)){
    
    plot_features<-features[FEATURE_SUBSET]
    N_CLUSTERS=setNames(rep(K_MEANS_N_CLUSTERS, length(plot_features)),plot_features)
    PNG<-FALSE
    
    heatmap_plot_data<-setNames(as.list(rep(NA, length(plot_features))), plot_features)
    
    #cluster_plots[[paste(plot_features, collapse="_")]]<-list()
    cluster_plots<-list()
    cluster_summaries<-list()
    cluster_avg_plots<-list()
    cluster_averages<-list()
    
    for(feat in plot_features){
        print(feat)
    #    for(mark in marks){
            
            sg_dfs<-do.call("rbind", final_heatmap_data[[feat]][plot_marks])
            sg_dfs$Mark<-factor(sg_dfs$Mark, levels=plot_marks)
            #sg_dfs_subs<-sg_dfs[sg_dfs$Group%in%plot_group,]
                
            sg_dfs_subs_matr<-sg_dfs %>% unite(Mark_Group, Mark, Group, Distance_kb) %>% spread(key=Mark_Group, value=value)#, "Feature", "Mark", "Group", "Group")
            #common_matr<-do.call("cbind", apply(expand.grid(plot_group, plot_marks), 1,  function(gmk) sg_dfs_subs_matr[sg_dfs_subs_matr$Mark==gmk[1] & sg_dfs_subs_matr$Group==gmk[2],setdiff(colnames(sg_dfs_subs_matr),c("Region","Feature", "Mark", "Group", "Group"))]))
            #common_rgns<-do.call("cbind", apply(expand.grid(plot_group, plot_marks), 1,  function(gmk) sg_dfs_subs_matr[sg_dfs_subs_matr$Mark==gmk[1] & sg_dfs_subs_matr$Group==gmk[2],"Region"]))
            common_matr<-as.matrix(sg_dfs_subs_matr[,setdiff(colnames(sg_dfs_subs_matr),c("Region","Feature", "Mark", "Group", "Group"))])
            regions<-sg_dfs_subs_matr[,"Region"]
            
            common_matr_clean<-common_matr
            common_matr_clean[is.na(common_matr_clean)]<-0
            
            kmc<-kmeans(common_matr_clean, centers=N_CLUSTERS[feat])
            
            clustr_ids<-kmc$cluster
            
            ### order clusters by size
            
            clustr_sizes<-table(clustr_ids)
            clustr_ids<-rank(-clustr_sizes, ties.method="first")[clustr_ids]
            clustr_sizes<-table(clustr_ids)
            valid_clusters<-which(clustr_sizes>=MIN_CLUSTER_SIZE)
            
            names(clustr_ids)<-regions
            
            clusters[[feat]]<-clustr_ids
            
            reg_by_cluster<-lapply(sort(unique(clustr_ids)), function(i) sort(as.integer(gsub("region_", "", unique(regions[clustr_ids==i])))))
            saveRDS(reg_by_cluster, file=file.path(STATE_RGNS_DIR, sprintf("kmeans_clustering_%s_%s.RDS", paste(plot_features, collapse="_"), paste(plot_marks, collapse="_"))))
            
            ranks<-clustr_ids*1000000 + rowSums(common_matr)
            
            sg_dfs$Region<-factor(sg_dfs$Region, levels=regions[order(ranks)])
            
            sg_dfs$Cluster<-clustr_ids[as.character(sg_dfs$Region)]
            
            
            if(CLUSTER_PROFILE_PLOTS){
                cids<-sort(unique(clustr_ids))
                for(cid in cids[valid_clusters]){
                    FEATURE_SUBSET<-
                            TARGET_CLUSTER<-2
                    region_subset<-names(clustr_ids[clustr_ids==cid])
                    
                    fdfcs<-sg_dfs[sg_dfs$Cluster == cid & sg_dfs$Feature==feat,]
                    fdfcs<-group_by(fdfcs, Group, Distance_kb, Feature, Mark, Cluster)
                    fdfcs_sum<-dplyr:::summarize(fdfcs, Mean=mean(value, na.rm=TRUE), Variance=var(value,na.rm=TRUE), SEM=var(value, na.rm=TRUE)/sqrt(n()))
                    
                    #cluster_plots[[paste(plot_features, collapse="_")]][[cid]]<-profile_meta_plot(fdfcs)
                    cluster_plots[[cid]]<-profile_meta_plot(fdfcs_sum, c("Cluster", "Mark"))
                    cluster_summaries[[cid]]<-fdfcs_sum
                }
            }
            
            if(CLUSTER_AVERAGE_PLOTS){
                cids<-sort(unique(clustr_ids))
                for(cid in cids[valid_clusters]){
                    FEATURE_SUBSET<-
                            TARGET_CLUSTER<-2
                    region_subset<-names(clustr_ids[clustr_ids==cid])
                    
                    final_reg_df_subs<-final_reg_df[final_reg_df$Region %in% region_subset & final_reg_df$Feature==feat,]
                    final_reg_df_subs$Cluster<-rep(cid, nrow(final_reg_df_subs))
                    #cluster_plots[[paste(plot_features, collapse="_")]][[cid]]<-plofile_meta_plot(fdfcs)
                    cluster_avg_plots[[cid]]<-region_average_scatterplot(final_reg_df_subs, 
                            #comp_groups=c("WT", "MT"), 
                            comp_marks=c("G34W", "K27me3"),
                            scatter_fit=FALSE)
                    cluster_averages[[cid]]<-final_reg_df_subs
                }
            }
            
            sg_dfs$value_log2<-log2(sg_dfs$value)
            sg_dfs$value_log2[sg_dfs$value_log2<0]<-0
            
            heatmap_plot_data[[feat]]<-sg_dfs
            #dev.off()
    #    }
    }
    heatmap_plot_data<-do.call("rbind", heatmap_plot_data)
    
    
    heatmap_plot_data<-heatmap_plot_data[heatmap_plot_data$Cluster %in% valid_clusters,]
    
    heatmap_plot_data$Group_Mark<-factor(paste(heatmap_plot_data$Group, heatmap_plot_data$Mark, sep="_"), 
            levels=apply(expand.grid(plot_marks, names(sample_groups))[,c(2:1)], 1, paste, collapse="_"))
    heatmap_plot_data$Feature_Cluster<-factor(paste(heatmap_plot_data$Feature, heatmap_plot_data$Cluster, sep="_"),
            levels=apply(expand.grid(plot_features, seq(max(heatmap_plot_data$Cluster))), 1, paste, collapse="_"))
    
    
    if(BINNED){
        heatmap_plot_data$Distance_kb<-factor(heatmap_plot_data$Distance_kb, levels=bin_locs/1000)
    }
    
    ### plot z-scores
    heatmap_plot_data_grp<- group_by(heatmap_plot_data, Group, Feature, Mark) # grouping the data by type
    heatmap_plot_data_grp <- mutate(heatmap_plot_data_grp, value_z = (value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE)) #groupwise standardization 
    heatmap_plot_data_grp$value_log<-log2(heatmap_plot_data_grp$value+1)
    #######
    
    #subplot_marks<-"G34W"
    
    #subplot_marks<-"K4me3"
    sub_plots<-list()
    for(subplot_mark in plot_marks){
    
        #hd2plot<-heatmap_plot_data_grp[heatmap_plot_data_grp$Feature %in% c("H3.3_MT"),]
        hd2plot<-heatmap_plot_data_grp[heatmap_plot_data_grp$Mark %in% subplot_mark, ]
        
        value2plot="value"
        
        QUANT<-0.01
        zlim <- quantile(hd2plot[[value2plot]], c(QUANT,1-QUANT), na.rm=TRUE)
        zmin<-zlim[1]
        zmax<-zlim[2]
        hd2plot[[value2plot]][hd2plot[[value2plot]]>zmax]<-zmax
        hd2plot[[value2plot]][hd2plot[[value2plot]]<zmin]<-zmin
        
        hd_grp2<- group_by(hd2plot,Region) # grouping the data by type
        hd_grp2 <- mutate(hd_grp2, region_sum=sum(value)) #groupwise standardization
        region_stats<-dplyr:::summarize(hd_grp2, region_sum=sum(value), region_mean=mean(value))
        hd2plot$Region_sorted<-factor(hd2plot$Region, levels=region_stats$Region[order(region_stats$region_sum, decreasing=FALSE)])
        
        
        p<-ggplot(hd2plot, aes(x=Distance_kb, y=Region, fill=value))
        p<-p+geom_raster()
        #p<-p+scale_fill_gradientn(
        #        #colours=rev(colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(10))#,
        #        colours=colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(10)#,
        #        #low="white", high="darkblue"
        #)
        #p<-p+scale_fill_gradient(
        #        #colours=rev(colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(10))#,
        #        #colours=colorRampPalette(brewer.pal(n = 7, name ="YlOrRd"))(10)#,
        #        low="blue", high="yellow"
        #)
        
        #p<-p+scale_fill_distiller(palette = "Spectral")
        if(BINNED){
            p<-p + scale_x_discrete(
                    #breaks=c(floor(n_flank_win/2), n_flank_win, n_flank_win+N_BINS, N_BINS+n_flank_win+ceiling(n_flank_win/2)),
                    breaks=c(levels(heatmap_plot_data$Distance_kb)[
                                    c(ceiling(n_flank_win/2),n_flank_win+1, N_BINS+n_flank_win, N_BINS+n_flank_win+ceiling(n_flank_win/2))]),
                    labels=c(
                           sprintf("%.1f",round(as.numeric(levels(heatmap_plot_data$Distance_kb)[ceiling(n_flank_win/2)]))),
                           BINNED_REG_START_NAME,BINNED_REG_END_NAME,
                           sprintf("%.1f",round(as.numeric(levels(heatmap_plot_data$Distance_kb)[N_BINS+n_flank_win+ceiling(n_flank_win/2)])))
                        ), 
                    expand = c(0, 0), position="bottom")
        }else{
            p<-p + scale_x_continuous(expand = c(0, 0))
        }
        brks<-round(unname(signif(0.95*zlim,1)), 1)
        names(brks)<-sprintf("%.1f",brks)
        lims<-unname(zlim)
        #p<-p+scale_fill_distiller(palette = "RdBu", direction=-1, breaks=brks, limits=brks)
        
        col6_scale <-  c("#00007F", "blue", "#007FFF", "cyan",
                              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
        col10_scale <- rev(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))
        
        p<-p+scale_fill_gradientn(colors=col10_scale, breaks=brks, limits=brks, guide=guide_legend(title.position="top",title.hjust =0.5))
        
        
        p<-p+ggpubr:::theme_pubr()
        p<-p+facet_grid(Feature_Cluster~Group_Mark,  scales="free", space="free_y")#nrow=length(plot_features),
        
        p<-p+theme(
                axis.text.x=if(BINNED) element_text(size = 10, angle = 45, vjust = 1, hjust = 1) else element_text(size=8),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                strip.background = element_blank(),
                strip.text.y = element_blank(),
                panel.spacing.x=unit(PANEL_SPACING1, "lines"),
                panel.spacing.y=unit(0.25, "lines"),
                plot.margin=margin(l=15, r=15),
                #legend.position="top",
                legend.title.align=0.5,
                #legend.key.size=unit(0.015, "npc")
                legend.key.width=unit(0.03, "npc")
        
        )

        #p<-p+geom_hline(yintercept=table(), linetype="dashed", 
        #        color = "red", size=1)
        
        #p<-p+ggtitle(sprintf("%s,%s", gsub("bivalent_ESC_cluster_combined_category_","",feat), paste(plot_marks, collapse="_")))
        
        p<-p + ylab(NULL)
        p<-p + xlab(NULL)
        p<-p+labs(fill=sprintf("%s\nnorm. cov.",subplot_mark))
        p<-p+guides(fill=guide_colourbar(title.position="top",title.hjust =0.5))
    #    if(BINNED){
    #        p<-p + xlab("Distance to/from anchor, kb")
    #    }else{
    #        p<-p + xlab("Distance to center, kb")
    #    }
    #    
        
        sub_plots[[subplot_mark]]<-p
    }
    
    #combined <- rbind(tableGrob(t(as.character(1:nc)), theme = ttheme_minimal(), rows = ""), 
    #        cbind(tableGrob(samples[1:nr], theme = ttheme_minimal()), 
    #                arrangeGrob(grobs = sub_plots,ncol=length(sub_plots)), size = "last"), size = "last")
    
    if(BINNED){
        xlab_string<-c("Distance from anchor, kb")
    }else{
        xlab_string<-sprintf("Distance to %s, kb", ANCHOR_NAME)
    }
    #    
    
    combined<-arrangeGrob(grobs = sub_plots,ncol=length(sub_plots), 
            bottom=xlab_string, padding = unit(PANEL_SPACING2, "line"))
    #combined<-do.call("grid.arrange", 
    #        c(as.list(unname(sub_plots)), ncol=length(sub_plots), bottom=xlab_string,  padding = unit(2, "line")))
    
    
    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/heatmaps_sorted_unique_%s_%d_%d_%s_%s_subplot", ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_"), paste(plot_marks, collapse="_"))
    h<-length(plot_features)
    w<-length(unique(heatmap_plot_data$Group))+length(plot_marks)
    if(PNG){
        png(paste0(fn, ".png"), width=100*(w+1), height=700*h)#, pointsize=5)
    }else{
        pdf(paste0(fn, ".pdf"), width=2*(w+1), height=7*h)
    }
    grid.newpage()
    grid.draw(combined)
    dev.off()
    
    if(CLUSTER_PROFILE_PLOTS){
        
        cluster_summaries<-do.call("rbind", cluster_summaries)
        cluster_summaries$Cluster<-as.factor(cluster_summaries$Cluster)
        if(BINNED){
            cluster_summaries$Distance_kb<-factor(cluster_summaries$Distance_kb, levels=bin_locs/1000)
        }
        pcc<-profile_meta_plot(cluster_summaries, facet_vars=c("Mark", "Cluster"), free_scales=CLUSTER_PROFILE_FREE_SCALE, 
                col_var=if(length(plot_marks)>1) "Mark" else "Cluster")
        
        combined_cluster<-arrangeGrob(grobs = cluster_plots,ncol=length(cluster_plots), padding = unit(4, "line"))
        #combined<-do.call("grid.arrange", 
        #        c(as.list(unname(sub_plots)), ncol=length(sub_plots), bottom=xlab_string,  padding = unit(2, "line")))
        fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/cluster_profile_plots_%s_%d_%d_%s_%s_subplot", ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_"), paste(plot_marks, collapse="_"))
        
        if(length(plot_marks)>1){
            h<-length(unique(cluster_summaries$Cluster))
            w<-length(unique(cluster_summaries$Mark))
        }else{
            h<-length(unique(cluster_summaries$Mark))
            w<-length(unique(cluster_summaries$Cluster))
        }
    
        if(PNG){
            png(paste0(fn, ".png"), width=100*(w+1), height=700*h)#, pointsize=5)
        }else{
            pdf(paste0(fn, ".pdf"), width=(1+CLUSTER_PROFILE_FREE_SCALE)*(w+1)+2, height=(1+CLUSTER_PROFILE_FREE_SCALE)*(h+1))
        }
        #grid.newpage()
        #grid.draw(combined_cluster)
        print(pcc)
        dev.off()
    }
    
    if(CLUSTER_AVERAGE_PLOTS){
        
        cluster_averages<-do.call("rbind", cluster_averages)
        
        for (comp_gr in c(list(NULL), GROUP_AVG_COMPARISONS)){
            for (comp_mk in c(list(NULL), MARK_AVG_COMPARISONS)){
                
                if(!(is.null(comp_gr) && is.null(comp_mk))){
                    if(is.null(comp_gr)){
                        plot_subname<-paste(comp_mk, collapse="_vs_")
                        h<-2*length(sample_groups)
                        w<-2*length(unique(cluster_averages$Cluster))
                    }else if(is.null(comp_mk)){
                        plot_subname<-paste(comp_gr, collapse="_vs_")
                        h<-2*length(unique(cluster_averages$Mark))
                        w<-2*length(unique(cluster_averages$Cluster))
                    }else{
                        plot_subname<-sprintf("%s_%s_vs_%s_%s", comp_mk[1], comp_gr[1], comp_mk[2], comp_gr[2])
                        h<-2*length(unique(cluster_averages$Cluster))
                        w<-2#length(unique(cluster_averages$Mark))
                    }
                    
                    pac<-region_average_scatterplot(
                            cluster_averages, comp_groups=comp_gr, comp_marks=comp_mk, facet_vars=c("Cluster"), scatter_fit=TRUE)
                    
                    #combined_cluster_avg<-arrangeGrob(grobs = cluster_avg_plots,ncol=length(cluster_avg_plots), padding = unit(4, "line"))
                    #combined<-do.call("grid.arrange", 
                    #        c(as.list(unname(sub_plots)), ncol=length(sub_plots), bottom=xlab_string,  padding = unit(2, "line")))

                    fn<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/cluster_average_plots_%s_%s_%d_%d_%s_subplot", 
                            plot_subname,ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_"))
                    if(PNG){
                        png(paste0(fn, ".png"), width=100*(w+1), height=700*h)#, pointsize=5)
                    }else{
                        pdf(paste0(fn, ".pdf"), width=1*(w+1), height=h+1)
                    }
                    #grid.newpage()
                    #grid.draw(combined_cluster_avg)
                    print(pac)
                    dev.off()
                }
            }
        }
        
    }
}



#prof_feat<-features[1]
#prof_mark<-"K27me3"
#prof_data<-t(do.call("rbind", final_swprofiles[[prof_feat]][[prof_mark]]))
#
#pdf("/ngs_share/scratch/pavlo/gctb/analysis/integrative/summaries/profile_plot.pdf")
#image(prof_data, axes=FALSE,xlab="",ylab="",srt=45)
#axis(2, at = (1.2/ncol(prof_data))*(1:ncol(prof_data))-0.2, labels=colnames(prof_data),tick=FALSE, las=2)
#dev.off()
#




