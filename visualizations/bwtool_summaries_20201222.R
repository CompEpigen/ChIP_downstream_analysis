#library(GenomicRanges)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(ggrastr)
#library(Hotelling)
#library(ICSNP)
#library(FRB)
library(EmpiricalBrownsMethod)
library(Hmisc)
library(ggpubr)
library(ggsci)

############ read bwtool matrices

args = commandArgs(trailingOnly=TRUE); 
source(args[1])

grps<-unique(GROUPS)
sample_groups<-lapply(grps, function(grp) SAMPLES[GROUPS==grp])
names(sample_groups)<-grps

if(SAMPLE_WISE){
  sample_groups<-unlist(mt_wt_sample_groups)
  names(sample_groups)<-sample_groups
  sample_groups<-as.list(sample_groups)
}

#sample_groups<-mt_wt_sample_groups

col_schemes<-ggsci:::pal_npg()(10)[seq(length(names(sample_groups)))]
group_colors<-setNames(col_schemes, names(sample_groups))

seq_depths<-data.frame()

#### for bwtools

if(BINNED){
  n_flank_win<-OFFSET/WIN_SIZE
  nwin<-2*OFFSET/WIN_SIZE + N_BINS
  bin_locs<-c((seq(-n_flank_win, 0.5)-0.5)[-1]*WIN_SIZE, 0, 1/seq(N_BINS)[-N_BINS] ,(seq(0, n_flank_win)-0.5)[-1]*WIN_SIZE)
}else{
  nwin<-2*OFFSET/WIN_SIZE
  bin_locs<-(seq(-nwin/2, nwin/2)-0.5)[-1]*WIN_SIZE
}

OUT_DIR=getwd()
if(!file.exists(OUT_DIR)){
  dir.create(OUT_DIR)
}

STATE_RGNS_DIR=file.path(OUT_DIR, "summary_files")
dir.create(STATE_RGNS_DIR)

features<-gsub("\\.bed", "", list.files(ANNOTATION_DIR, pattern="*.bed$"))


  
  bw_files<-gsub("\\.bw", "", list.files(BW_DIR, pattern="bw$"))
  
  ########### extract profile matrices with BWtools
  
  bwtool_call<-paste0(
    if(BINNED)
      "/02_code/bwtool/bwtool matrix ${UPSTREAM}:${META}:${DOWNSTREAM} \\\n"
    else
      "/02_code/bwtool/bwtool matrix ${UPSTREAM}:${DOWNSTREAM} ${ANCHOR} \\\n",
    "$ANNOTATION_DIR/${rgn}.bed $BW_DIR/$bwig \\\n",
    "$BWTOOL_OUT_DIR/summary_${rgn}_${TOKEN}_${UPSTREAM}_${DOWNSTREAM}_${TILE_AVG}_${samp}.txt \\\n",
    "-tiled-averages=${TILE_AVG}",sep=""
  )
  
  print(bwtool_call)
  
  command<-paste(
    sprintf("ANALYSIS=%s\n",ANALYSIS),
    sprintf("ANNOTATION_DIR=%s\n",ANNOTATION_DIR),
    sprintf("BW_DIR=%s\n",BW_DIR),
    sprintf("BWTOOL_OUT_DIR=%s\n",STATE_RGNS_DIR),
    sprintf("UPSTREAM=%d\n",OFFSET),
    sprintf("ANCHOR=%s\n", if(ANCHOR_SPECS=="start") "-starts" else if(ANCHOR_SPECS=="end") "-ends" else ""),
    sprintf("DOWNSTREAM=%d\n",OFFSET),
    sprintf("TILE_AVG=%d\n",WIN_SIZE),
    sprintf("TMP_DIR=%s\n","/ngs_share/tmp"),
    sprintf("REGIONS=\"%s\"\n", paste(features, collapse=" ")),
    sprintf("TOKEN=%s\n", if(BINNED) "binned" else ANCHOR_SPECS),
    sprintf("META=%s\n", N_BINS ),
    "mkdir $TMP_DIR/logs/$ANALYSIS\n",
    sprintf("INPUT_BWS=`ls -1 $BW_DIR | grep -v log | grep %s`\n", BW_EXTENSION),
    "for rgn in `echo $REGIONS`\n",
    "do\n",
    "for bwig in `echo $INPUT_BWS`\n",
    "do\n",
    sprintf("samp=${bwig%%.%s}\n", BW_EXTENSION),
    bwtool_call,
    " &\n",
    "done\n",
"done\n",
sep="")
  
print(command)

system(command, intern=FALSE)




repeat{
  if(QSUB){
        lookup_cmd<-sprintf("qstat -r | grep \"Full jobname\" | grep -e %s_%s", "bwtool", ANALYSIS)
  }else{
        lookup_cmd<-sprintf("COLUMNS=10000 ps aux | grep -e %s | grep -v grep | grep -v COLUMNS", ANALYSIS)
  }
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

################################################
################ read the extracted matrices
################################################

final_data<-final_heatmap_data<-final_reg_data<-final_swprofiles<-final_rmats<-setNames(as.list(rep(NA, length(features))), features)
region_ids<-feature_grs<-list()

for(feature in features){
  print(feature)
  final_data[[feature]]<-final_heatmap_data[[feature]]<-final_reg_data[[feature]]<-final_swprofiles[[feature]]<-final_rmats[[feature]]<-list()
  
  fbed<-read.table(file.path(ANNOTATION_DIR, paste0(feature,".bed")))
  #feature_grs[[feature]]<-GRanges(fbed$V1, IRanges(fbed$V2, fbed$V3), "*")
  #region_ids[[feature]]<-paste0("region_", seq_along(feature_grs[[feature]]))
  region_ids[[feature]]<-paste0("region_", seq_along(fbed))
  
  mat_list<-list()
  for(sample in unlist(sample_groups)){
    for(mark in c(marks)){
      mat_file<-file.path(STATE_RGNS_DIR, sprintf(
        "summary_%s_%s_%d_%d_%d_%s%s_%s%s.txt",
        feature, if(BINNED) "binned" else ANCHOR_SPECS, OFFSET, OFFSET,WIN_SIZE, BW_PREFIX, sample,mark,BW_SUFFIX))
      if(file.exists(mat_file)){
        mat_list[[paste(sample, mark, sep="_")]]<-read.table(mat_file)
      }else{
        print(paste(mat_file, "does not exist"))
      }
    }
  }
  
  print(lapply(mat_list, dim))
  
  nreg<-nrow(mat_list[[1]])
  
  rand_subset<-sort(sample.int(nreg, min(nreg,5000)))
  
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
          
        }else{
          sg_profiles[[sgn]][[sample]]<-rep(NA,nwin)
          sg_averages[[sgn]][[sample]]<-rep(NA,nreg)
          
        }
      }
      imat<-do.call("rbind", sg_profiles[[sgn]])
      sw_profiles[[sgn]]<-imat[which(apply(!is.na(imat),1,any)),,drop=FALSE]
      sg_profiles[[sgn]]<-colMeans(imat, na.rm=TRUE)
      sg_vars[[sgn]]<-apply(imat, 2, sd, na.rm=TRUE)
      sg_sems[[sgn]]<-sg_vars[[sgn]]/sqrt(present_samples)
      
      if(AVERAGES){
        rmat<-do.call("cbind", sg_averages[[sgn]])
        sg_rmats[[sgn]]<-rmat[,which(apply(!is.na(rmat),2,any))]
        sg_averages[[sgn]]<-rowMeans(rmat, na.rm=TRUE)
      }
      
      if(HEATMAP){
        sg_matrix<-sg_matrix/present_samples
        sg_matrices[[sgn]]<-sg_matrix
        
        if(!all(as.numeric(as.matrix(sg_matrix))==0)){
          rownames(sg_matrix)<-region_ids[[feature]]
          sg_matrix<-sg_matrix[order(rowSums(sg_matrix), decreasing=TRUE),]
          sg_df<-cbind(data.frame(Distance_kb=bin_locs/1000), t(sg_matrix))
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
    final_swprofiles[[feature]][[mark]]<-sw_profiles
    
    sg_profiles<-do.call("rbind", sg_profiles)
    sg_vars<-do.call("rbind", sg_vars)
    sg_sems<-do.call("rbind", sg_sems)
    rownames(sg_vars)<-paste0("SD_", names(sample_groups))
    rownames(sg_sems)<-paste0("SEM_", names(sample_groups))
    
    rownames(sg_profiles)<-paste0("average_", names(sample_groups))
    final_data[[feature]][[mark]]<-cbind(
      data.frame(
        Feature=rep(feature, ncol(sg_profiles)), 
        Mark=rep(mark, ncol(sg_profiles)),
        Distance=bin_locs),
      as.data.frame(t(sg_profiles)),
      as.data.frame(t(sg_vars)),
      as.data.frame(t(sg_sems))
    )
    if(AVERAGES){
      sg_averages<-do.call("cbind", sg_averages)
      colnames(sg_averages)<-paste0("average_", names(sample_groups))
      final_reg_data[[feature]][[mark]]<-cbind(
        data.frame(
          Feature=rep(feature, nrow(sg_averages)), 
          Mark=rep(mark, nrow(sg_averages)),
          Region=region_ids[[feature]]),
        as.data.frame(sg_averages)
      )
      final_rmats[[feature]][[mark]]<-sg_rmats
    }
  }
}
saveRDS(final_swprofiles, file=file.path(OUT_DIR, "final_swprofiles.RDS"))
saveRDS(final_data, file=file.path(OUT_DIR, "final_data.RDS"))
saveRDS(final_reg_data, file=file.path(OUT_DIR, "final_reg_data.RDS"))
saveRDS(final_rmats, file=file.path(OUT_DIR, "final_rmats.RDS"))
saveRDS(final_heatmap_data, file=file.path(OUT_DIR, "final_heatmap_data.RDS"))

if(FALSE){
  ANALYSIS_GROUP<-"Bivalent_genes"
  
  STATE_RGNS_DIR2=file.path(OUT_DIR,
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
      print(length(nna))
      
      if(length(nna)>1){
        nzsd1<-which(apply(grp1, 2, sd)>0)
        nzsd2<-which(apply(grp2, 2, sd)>0)
        grp1<-grp1[,nna]
        grp2<-grp2[,nna]
        
        p_vals<-sapply(seq(nna), function(i) t.test(grp1[,i], grp2[,i])$p.value)
        t.test.pvals[[sprintf("%s_%s", feat, mark)]]<-p_vals
      }else{
        t.test.pvals[[sprintf("%s_%s", feat, mark)]]<-c(1.0)
      }
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


######################################################################## 
########################   PROFILE PLOTS

final_df<-do.call("rbind",unlist(final_data, recursive=FALSE))
rownames(final_df)<-NULL

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

final_df_clean$Distance_kb<-final_df_clean$Distance/1000

if(BINNED){
  final_df_clean$Distance_kb<-factor(final_df_clean$Distance_kb, levels=bin_locs/1000)
}


feat_column<-unlist(lapply(names(final_swprofiles), function(feat) rep(feat,length(final_swprofiles[[feat]]))))
mk_column<-unlist(lapply(names(final_swprofiles), function(feat) lapply(names(final_swprofiles[[feat]]), function(mk) rep(mk,length(bin_locs)))))
final_swprofiles_df<-do.call("rbind", setNames(unlist(unlist(final_swprofiles, recursive=FALSE), recursive=FALSE), NULL))

final_swprofiles_df<-list()
for(feat in names(final_swprofiles)){
  for(mk in marks){
    for(sgn in names(sample_groups)){
      if(is.matrix(final_swprofiles[[feat]][[mk]][[sgn]])){
        final_swprofiles_df[[length(final_swprofiles_df)+1]]<-data.frame(final_swprofiles[[feat]][[mk]][[sgn]], Feature=feat, Mark=mk, Group=sgn, Sample=rownames(final_swprofiles[[feat]][[mk]][[sgn]]))
      }else{
        final_swprofiles_df[[length(final_swprofiles_df)+1]]<-data.frame(as.data.frame(t(final_swprofiles[[feat]][[mk]][[sgn]])), Feature=feat, Mark=mk, Group=sgn, Sample=rownames(final_swprofiles[[feat]][[mk]][[sgn]]))
      }
    }
  }
}
final_swprofiles_df<-do.call("rbind", final_swprofiles_df)
final_swprofiles_df<-gather(final_swprofiles_df, key="Distance", value="Mean", -Feature, -Mark, -Group,-Sample)

final_swprofiles_df$Distance<-setNames(bin_locs, paste0("V", seq(bin_locs)))[final_swprofiles_df$Distance]
final_swprofiles_df$Distance_kb<-final_swprofiles_df$Distance/1000

profile_meta_plot<-function(summarized_df, facet_vars=c("Mark", "Feature"), free_scales=FALSE, col_var="Feature", pval_df=NULL, per_sample=FALSE, binned=FALSE, error_measure="SD"){
  
  summarized_df$Error<-summarized_df[[error_measure]]
  
  if(per_sample){
    p<-ggplot(summarized_df, aes(x=Distance_kb, y=Mean, color=Group, group=Sample)) 
  }else{
    p<-ggplot(summarized_df, aes(x=Distance_kb, y=Mean, color=Group))
  }
  
  if(BINNED){
    p<-p + geom_line(lwd=1, aes(group=Group), alpha=1 - 0.5*per_sample)
  }else{
    p<-p + geom_line(lwd=1, alpha= 1 - 0.5*per_sample)
  }
  
  if(per_sample){
    p<-p + geom_point(aes(shape=Sample),size=1)
  }
  
  if(!per_sample & !binned){
    p<-p + geom_ribbon(aes(ymin=Mean-Error, ymax=Mean+Error, fill=Group),alpha=0.3, lwd=0.25)
  }
  
  p<-p + scale_colour_manual(values=group_colors)
  
  if(binned){
    p<-p + scale_x_discrete(
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
    p<-p + facet_wrap(as.formula(sprintf("%s~%s",facet_vars[1],facet_vars[2])), scales="free", 
                      ncol=length(unique(summarized_df[[col_var]])), drop=TRUE)
  }else{
    p<-p + facet_grid(as.formula(sprintf("%s~%s",facet_vars[1],facet_vars[2])), scales="free", drop=TRUE)
  }
  
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
  
  p<-p + theme_classic()
  
  return(p)
}

final_df_clean_subs<-final_df_clean
htsq_subs<-htsq
htsq_subs<-NULL

profiles_mark_subset<-marks

final_df_clean_subs<-final_df_clean[final_df_clean$Mark %in% profiles_mark_subset,]
if(!is.null(htsq_subs)) htsq_subs<-htsq[ htsq$Mark %in% profiles_mark_subset,]
final_df_clean_subs$Mark<-factor(final_df_clean_subs$Mark, levels=profiles_mark_subset)
if(!is.null(htsq_subs)) htsq_subs$Mark<-factor(htsq_subs$Mark, levels=profiles_mark_subset)


pdf(file.path(OUT_DIR, sprintf("profile_summaries_SEM_%s_pval_subs.pdf", 
                               paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))),
    width=1*length(unique(final_df_clean_subs$Feature))+3.5, height=1*length(unique(final_df_clean_subs$Mark))+2)# width=20, height=7)
print(profile_meta_plot(final_df_clean_subs, pval_df=htsq_subs, binned=BINNED))
dev.off()


final_swprofiles_df_subs<-final_swprofiles_df[final_swprofiles_df$Mark %in% profiles_mark_subset,]
final_swprofiles_df_subs$Mark<-factor(final_swprofiles_df_subs$Mark, levels=profiles_mark_subset)


pdf(file.path(OUT_DIR,sprintf("profile_summaries_SEM_%s_pval_subs_samplewise.pdf", 
                              paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"))),
    width=1*length(unique(final_swprofiles_df_subs$Feature))+3.5, height=1*length(unique(final_swprofiles_df_subs$Mark))+2)# width=20, height=7)
print(profile_meta_plot(final_swprofiles_df_subs, pval_df=NULL, per_sample=TRUE))
dev.off()


######## all final region-level summaries
final_reg_df<-do.call("rbind",unlist(final_reg_data, recursive=FALSE))
rownames(final_reg_df)<-NULL

final_df_clean_subs_diff<-final_df_clean_subs[,c("Feature","Mark","Distance","Mean","Distance_kb", "Group")]
final_df_clean_subs_diff<-spread(final_df_clean_subs_diff, key="Group", value="Mean")
final_df_clean_subs_diff$Delta<-final_df_clean_subs_diff[[comp_groups[2]]]-final_df_clean_subs_diff[[comp_groups[2]]]

################################################

if(AVG_SCATTERPLOTS || AVG_BOXPLOTS || SEGMENT_CHANGE_PLOTS){
  
final_reg_df$Feature<-gsub("Feature_", "", final_reg_df$Feature)
final_reg_df$Feature<-factor(final_reg_df$Feature,features)
}
  
######################################################################## 
########################   AVERAGE SCATTERPLOTS
if(AVG_SCATTERPLOTS){
  
  region_average_scatterplot<-function(final_reg_df, comp_groups=NULL, comp_marks=NULL, facet_vars=NULL, scatter_fit=TRUE, scales_type="free"){
    
    
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
      
      p<-ggplot(final_reg_df, aes_string(x=sprintf("average_%s", comp_groups[1]), y=sprintf("average_%s", comp_groups[2]))) + geom_point_rast(size=0.5, alpha=0.5)
      
      p<-p+geom_abline(slope=1,intercept=0,linetype="dashed")
      p<-p + coord_symmetric()
      
      p<-p + facet_wrap(c("Mark","Feature", facet_vars), 
                        ncol=length(unique(final_reg_df$Feature))*(fvl+(fvl==0)), scales=scales_type)
      
      p<-p + theme_bw()
      
    }
    ################ different marks vs each other
    else if(is.null(comp_groups)){
      
      final_reg_df_by_mark<-spread(final_reg_df_clean, key=Mark, value=Mean)
      
      p<-ggplot(final_reg_df_by_mark, aes_string(x=comp_marks[1], y=comp_marks[2])) + geom_point_rast(size=0.5, alpha=0.5)
      
      p<-p+geom_abline(slope=1,intercept=0,linetype="dashed")
      
      p<-p + facet_wrap(c("Group","Feature",facet_vars), 
                        ncol=length(unique(final_reg_df$Feature))*(fvl+(fvl==0)+is.null(fvl)), scales=scales_type)
      
      p<-p + theme_bw()
      
    }
    ######### changes vs each other
    else if(!is.null(comp_groups) & !is.null(comp_marks)){
      
      final_reg_diffs<-final_reg_df_clean %>% unite(Mark_Group, Mark, Group) %>% spread(Mark_Group, Mean)
      
      final_reg_diffs$Diff_mark1<-final_reg_diffs[[sprintf("%s_%s", comp_marks[1], comp_groups[2])]]-final_reg_diffs[[sprintf("%s_%s", comp_marks[1], comp_groups[1])]]
      final_reg_diffs$Diff_mark2<-final_reg_diffs[[sprintf("%s_%s", comp_marks[2], comp_groups[2])]]-final_reg_diffs[[sprintf("%s_%s", comp_marks[2], comp_groups[1])]]
      
      p<-ggplot(final_reg_diffs, aes(x=Diff_mark1, y=Diff_mark2)) + geom_point_rast(size=0.5, alpha=0.5)
      p<-p+geom_vline(xintercept=0,linetype="dashed")+geom_hline(yintercept=0,linetype="dashed")
      
      p<-p + xlab(sprintf("%s, %s - %s",  comp_marks[1], comp_groups[2], comp_groups[1]))
      p<-p + ylab(sprintf("%s, %s - %s",  comp_marks[2], comp_groups[2], comp_groups[1]))
      
      p<-p + facet_wrap(c("Feature", facet_vars), ncol=length(unique(final_reg_df$Feature)))
      
      p<-p + theme_bw()
      
      if(scatter_fit){
        
        lm_eqn <- function(df, y, x){
          formula = as.formula(sprintf('%s ~ %s', y, x))
          m <- lm(formula, data=df);
          eq <- substitute(
            italic(r)^2~"="~r2*","~~p~"="~italic(pvalue), 
            list(target = "y",
                 input = "x",
                 a = format(as.vector(coef(m)[1]), digits = 2), 
                 b = format(as.vector(coef(m)[2]), digits = 2), 
                 r2 = format(summary(m)$r.squared, digits = 3),
                 pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=1)
            )
          )
          as.character(as.expression(eq));
        }
        
        eq<-plyr:::ddply(final_reg_diffs,c("Feature", facet_vars), lm_eqn, y='Diff_mark1', x='Diff_mark2')
        p<-p+geom_smooth(method = "lm", se = TRUE, color='red', size=1, linetype=2)
        p<-p+geom_text(data=eq, 
                       aes(
                         x=-Inf, y=Inf,
                         label=V1), 
                       vjust=1, hjust=0, color='black', parse=T)
      }
    }
    return(p)
  }
  
  
  fn<-file.path(OUT_DIR, sprintf("region_scatterplots_%s_%s_vs_%s", 
                                 paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_"), comp_groups[1], comp_groups[2]))
  if(PNG){
    png(paste0(fn, ".png"),
        width=250*length(unique(final_reg_df$Feature))+200, height=250*length(unique(final_reg_df$Mark)))# width=20, height=7)
  }else{
    pdf(paste0(fn, ".pdf"),
        width=2.5*length(unique(final_reg_df$Feature))+2, height=2.5*length(unique(final_reg_df$Mark)))# width=20, height=7))
  }
  print(region_average_scatterplot(final_reg_df, comp_groups=comp_groups))
  dev.off()
  
  
  if(!is.null(comp_marks)){
    
    fn<-file.path(OUT_DIR, sprintf("region_scatterplots_%s_vs_%s_%s", comp_marks[1], comp_marks[2],
                                   paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")))
    if(PNG){
      png(paste0(fn, ".png"),
          width=200*length(unique(final_reg_df$Feature))+200, height=200*length(unique(final_reg_df$Group)))# width=20, height=7)
    }else{
      pdf(paste0(fn, ".pdf"), width=1.5*length(unique(final_reg_df$Feature))+2, height=1.5*length(grep("average_", colnames(final_reg_df))))
    }
    print(region_average_scatterplot(final_reg_df, comp_marks=comp_marks))
    dev.off()
    
    fn<-file.path(OUT_DIR, sprintf("region_scatterplot_diffs_%s_vs_%s_in_%s_vs_%s_%s", 
                                   comp_marks[1], comp_marks[2],comp_groups[1], comp_groups[2],
                                   paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")))
    if(PNG){
      png(paste0(fn, ".png"),
          width=200*length(unique(final_reg_df$Feature))+200, height=200)# width=20, height=7)
    }else{
      pdf(paste0(fn, ".pdf"), width=2*length(unique(final_reg_df$Feature))+2, height=4)
    }
    print(region_average_scatterplot(final_reg_df, comp_groups=comp_groups, comp_marks=comp_marks))
    dev.off()
    
  }
}


####################### AVERAGE BOXPLOTS

if(AVG_BOXPLOTS){
  
  ################### define the routine
  
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
      plot_var_x<-"Group"
      plot_var_y<-"Coverage"
    }else{
      final_reg_df_bp<-final_reg_df
      final_reg_df_bp$Difference<-final_reg_df_bp[[sprintf("average_%s", comp_groups[2])]]-final_reg_df_bp[[sprintf("average_%s", comp_groups[1])]]
      
      if(is.null(pval_df)){
        
        if(is.null(ref_feature)){
          if(test=="t.test"){
            pvals <- final_reg_df_bp %>% group_by_(.dots=c("Mark", "Feature", facet_vars)) %>%
              dplyr::summarize(p.value = (t.test(Difference, mu = 0)$p.value))
          }else if(test=="wilcox.test"){
            pvals <- final_reg_df_bp %>% group_by_(.dots=c("Mark", "Feature", facet_vars)) %>%
              dplyr::summarize(p.value = (wilcox.test(Difference, mu = 0)$p.value))
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
          pvals<-do.call("rbind", pvals)
        }
      }else{
        pvals<-pval_df
      }
      plot_var_x<-sprintf("'%s-%s'", comp_groups[2], comp_groups[1])
      plot_var_y<-"Difference"
      
      
      pvals$p_val_formatted=sprintf("P=%1.2g", pvals$p.value)
      
    }
    
    p<-ggplot(final_reg_df_bp, 
              aes_string(x=plot_var_x, y=plot_var_y)) + geom_boxplot()
    
    if(!is.null(comp_groups)){
      p<-p + geom_hline(yintercept=0, linetype=2)
    }
    
    p<-p + facet_wrap(c("Mark","Feature", facet_vars), 
                      ncol=length(unique(final_reg_df_bp$Feature))*(fvl+(fvl==0)), scales="free")
    
    if(is.null(comp_groups)) {
      p<-p + xlab("Group means")
    }else{
      p<-p + xlab("Delta of means")
      p<-p + geom_text(mapping=aes(x=-Inf, y=0, label=p_val_formatted), data=pvals#, 
                       ,hjust=-0.5
                       ,vjust=-1
      )
    }
    
    p<-p + theme_bw()
    
    p
  }
  
  
  ################### generate plots
  
  fn<-file.path(OUT_DIR, sprintf("region_boxplot_diffs_%s", 
                                 paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")))
  if(PNG){
    png(paste0(fn, ".png"),
        width=250*length(unique(final_reg_df$Feature))+200, height=250*length(unique(final_reg_df$Mark)))
  }else{
    pdf(paste0(fn, ".pdf"),
        width=2.5*length(unique(final_reg_df$Feature))+2, height=2.5*length(unique(final_reg_df$Mark)))
  }
  print(region_average_boxplot(final_reg_df))
  dev.off()
  
  fn<-file.path(OUT_DIR, sprintf("region_boxplot_diffs_%s_vs_%s_%s",
                                 comp_groups[2], comp_groups[1],
                                 paste(ANALYSIS, c("centered", "binned")[BINNED+1], OFFSET, WIN_SIZE, sep="_")))
  if(PNG){
    png(paste0(fn, ".png"),
        width=150*length(unique(final_reg_df$Feature))+100, height=150*length(unique(final_reg_df$Mark)))
  }else{
    pdf(paste0(fn, ".pdf"),
        width=1.5*length(unique(final_reg_df$Feature))+1, height=1.5*length(unique(final_reg_df$Mark)))
  }
  print(region_average_boxplot(final_reg_df, comp_groups=comp_groups, ref_feature=AVG_DELTA_REF_FEATURE, pval_df=htsq))
  dev.off()
}


######################################################## selected marks

cluster_marks<-plot_marks[seq(length(plot_marks))]
sort_marks<-plot_marks[1]

clusters<-list()

feature_seq<-seq(features)

for(FEATURE_SUBSET in list(feature_seq)){
  
  plot_features<-features[FEATURE_SUBSET]
  N_CLUSTERS=setNames(rep(K_MEANS_N_CLUSTERS, length(plot_features)),plot_features)
  
  heatmap_plot_data<-setNames(as.list(rep(NA, length(plot_features))), plot_features)
  
  cluster_plots<-list()
  cluster_summaries<-list()
  cluster_avg_plots<-list()
  cluster_averages<-list()
  
  for(feat in plot_features){
    print(feat)
    
    sg_dfs<-do.call("rbind", final_heatmap_data[[feat]][plot_marks])
    sg_dfs$Mark<-factor(sg_dfs$Mark, levels=plot_marks)
    
    sg_dfs_subs_matr<-sg_dfs %>% unite(Mark_Group, Mark, Group, Distance_kb) %>% spread(key=Mark_Group, value=value)#, "Feature", "Mark", "Group", "Group")
    
    common_matr<-as.matrix(sg_dfs_subs_matr[,setdiff(colnames(sg_dfs_subs_matr),c("Region","Feature", "Mark", "Group", "Group"))])
    regions<-sg_dfs_subs_matr[,"Region"]
    
    common_matr_clean<-common_matr
    common_matr_clean[is.na(common_matr_clean)]<-0
    
    cl_marks_col_ids<-grep(paste(paste0("^", cluster_marks,"_"), collapse="|"), colnames(common_matr_clean))
    kmc<-kmeans(common_matr_clean[,cl_marks_col_ids,drop=FALSE], centers=N_CLUSTERS[feat])
    
    clustr_ids<-kmc$cluster
    
    ### order clusters by size
    
    clustr_sizes<-table(clustr_ids)
    clustr_ids<-rank(-clustr_sizes, ties.method="first")[clustr_ids]
    clustr_sizes<-table(clustr_ids)
    valid_clusters<-which(clustr_sizes>=MIN_CLUSTER_SIZE)
    clustr_sizes<-as.numeric(clustr_sizes)
    
    names(clustr_ids)<-regions
    
    clusters[[feat]]<-clustr_ids
    
    reg_by_cluster<-lapply(sort(unique(clustr_ids)), function(i) sort(as.integer(gsub("region_", "", unique(regions[clustr_ids==i])))))
    saveRDS(reg_by_cluster, file=file.path(OUT_DIR, sprintf("kmeans_clustering_%s_%s.RDS", feat, paste(plot_marks, collapse="_"))))
    
    #### sort the regions 
    
    sort_col_ids<-grep(paste(paste0("^", sort_marks,"_"), collapse="|"), colnames(common_matr))
    
    ranks<-clustr_ids*1000000 + rowSums(common_matr[,sort_col_ids,drop=FALSE])
    
    sg_dfs$Region<-factor(sg_dfs$Region, levels=regions[order(ranks)])
    
    sg_dfs$Cluster<-clustr_ids[as.character(sg_dfs$Region)]
    
    sg_dfs$Cluster_description<-sprintf("Cluster %d, %d, %s", 
                                        sg_dfs$Cluster, clustr_sizes[sg_dfs$Cluster], scales::percent(clustr_sizes[sg_dfs$Cluster]/sum(clustr_sizes)))
    
    if(CLUSTER_PROFILE_PLOTS){
      cids<-sort(unique(clustr_ids))
      for(cid in cids[valid_clusters]){
        region_subset<-names(clustr_ids[clustr_ids==cid])
        
        fdfcs<-sg_dfs[sg_dfs$Cluster == cid & sg_dfs$Feature==feat,]
        fdfcs<-group_by(fdfcs, Group, Distance_kb, Feature, Mark, Cluster, Cluster_description)
        fdfcs_sum<-dplyr:::summarize(fdfcs, Mean=mean(value, na.rm=TRUE), SD=sd(value,na.rm=TRUE), SEM=var(value, na.rm=TRUE)/sqrt(n()))
        
        cluster_plots[[cid]]<-profile_meta_plot(fdfcs_sum, c("Cluster", "Mark"))
        cluster_summaries[[length(cluster_summaries)+1]]<-fdfcs_sum
      }
    }
    
    if(CLUSTER_AVERAGE_PLOTS || CLUSTER_AVERAGE_BOXPLOTS){
      cids<-sort(unique(clustr_ids))
      for(cid in cids[valid_clusters]){
        region_subset<-names(clustr_ids[clustr_ids==cid])
        
        final_reg_df_subs<-final_reg_df[final_reg_df$Region %in% region_subset & final_reg_df$Feature==feat,]
        final_reg_df_subs$Cluster<-rep(cid, nrow(final_reg_df_subs))
        cluster_avg_plots[[length(cluster_avg_plots)+1]]<-region_average_scatterplot(final_reg_df_subs, 
                                                                                     comp_marks=comp_marks,
                                                                                     scatter_fit=FALSE)
        cluster_averages[[length(cluster_averages)+1]]<-final_reg_df_subs
      }
    }
    sg_dfs$value_log2<-log2(sg_dfs$value)
    sg_dfs$value_log2[sg_dfs$value_log2<0]<-0
    
    heatmap_plot_data[[feat]]<-sg_dfs
  }
  heatmap_plot_data<-do.call("rbind", heatmap_plot_data)
  
  heatmap_plot_data<-heatmap_plot_data[heatmap_plot_data$Cluster %in% valid_clusters,]
  
  heatmap_plot_data$Group_Mark<-factor(paste(heatmap_plot_data$Group, heatmap_plot_data$Mark, sep="_"), 
                                       levels=apply(expand.grid(plot_marks, names(sample_groups))[,c(2:1)], 1, paste, collapse="_"))
  heatmap_plot_data$Feature_Cluster<-factor(paste(heatmap_plot_data$Feature, heatmap_plot_data$Cluster, sep="_"),
                                            levels=apply(expand.grid(seq(max(heatmap_plot_data$Cluster)), plot_features)[,c(2:1)], 1, paste, collapse="_"))
  
  
  if(BINNED){
    heatmap_plot_data$Distance_kb<-factor(heatmap_plot_data$Distance_kb, levels=bin_locs/1000)
  }
  
  ### plot z-scores
  heatmap_plot_data_grp<- group_by(heatmap_plot_data, Group, Feature, Mark) # grouping the data by type
  heatmap_plot_data_grp <- mutate(heatmap_plot_data_grp, value_z = (value-mean(value,na.rm=TRUE))/sd(value,na.rm=TRUE)) #groupwise standardization 
  heatmap_plot_data_grp$value_log<-log2(heatmap_plot_data_grp$value+1)
  #######
  
  sub_plots<-list()
  for(subplot_mark in plot_marks){
    
    hd2plot<-heatmap_plot_data_grp[heatmap_plot_data_grp$Mark %in% subplot_mark, ]
    
    value2plot="value"
    
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
    p<-p+geom_tile_rast()
    
    if(BINNED){
      p<-p + scale_x_discrete(
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
    brks<-round(unname(zlim),1)
    names(brks)<-sprintf("%.1f",brks)
    lims<-unname(zlim)
    
    print("======Breaks calculation=======")
    print(zlim)
    print(brks)
    
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
      panel.spacing.x=unit(PANEL_SPACING1, "lines"),
      panel.spacing.y=unit(0.25, "lines"),
      plot.margin=margin(l=15, r=15),
      legend.title.align=0.5,
      legend.key.width=unit(0.03, "npc")
      
    )
    
    p<-p + ylab(NULL)
    p<-p + xlab(NULL)
    p<-p+labs(fill=sprintf("%s\nnorm. cov.",subplot_mark))
    p<-p+guides(fill=guide_colourbar(title.position="top",title.hjust =0.5))
    
    sub_plots[[subplot_mark]]<-p
  }
  
  
  if(BINNED){
    xlab_string<-c("Distance from anchor, kb")
  }else{
    xlab_string<-sprintf("Distance to %s, kb", ANCHOR_NAME)
  }
  
  combined<-arrangeGrob(grobs = sub_plots,ncol=length(sub_plots), 
                        bottom=xlab_string, padding = unit(PANEL_SPACING2, "line"))
  
  
  fn<-file.path(OUT_DIR, sprintf("heatmaps_sorted_unique_%s_%d_%d_%s_%s_subplot", 
                                 ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_"), paste(plot_marks, collapse="_")))
  h<-length(plot_features)
  w<-length(unique(heatmap_plot_data$Group))+length(plot_marks)
  if(PNG){
    png(paste0(fn, ".png"), width=100*(w+1), height=700*h)#, pointsize=5)
  }else{
    pdf(paste0(fn, ".pdf"), width=2*(w+1), height=6+1*h)
  }
  grid.newpage()
  grid.draw(combined)
  dev.off()
  
  if(CLUSTER_PROFILE_PLOTS){
    
    cluster_summaries<-do.call("rbind", cluster_summaries)
    cluster_summaries$Cluster<-as.factor(cluster_summaries$Cluster)
    cluster_summaries$Feature_Cluster_description<-factor(paste(cluster_summaries$Feature, cluster_summaries$Cluster_description, sep=": "),
                                                          levels=apply(expand.grid(unique(heatmap_plot_data$Cluster_description), plot_features)[,c(2:1)], 1, paste, collapse=": "))
    
    if(BINNED){
      cluster_summaries$Distance_kb<-factor(cluster_summaries$Distance_kb, levels=bin_locs/1000)
    }
    
    
    pcc<-profile_meta_plot(cluster_summaries, facet_vars=c("Mark", "Feature_Cluster_description"), 
                           free_scales=CLUSTER_PROFILE_FREE_SCALE, 
                           col_var=if(length(plot_marks)==1) "Mark" else "Feature_Cluster_description", error_measure=ERROR_MEASURE)
    
    combined_cluster<-arrangeGrob(grobs = cluster_plots,ncol=length(cluster_plots), padding = unit(4, "line"))
    
    fn<-file.path(OUT_DIR, sprintf("cluster_profile_plots_%s_%d_%d_%s_%s_subplot", 
                                   ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_"), paste(plot_marks, collapse="_")))
    
    if(length(plot_marks)==1){
      h<-length(unique(cluster_summaries$Feature_Cluster_description))
      w<-length(unique(cluster_summaries$Mark))
    }else{
      h<-length(unique(cluster_summaries$Mark))
      w<-length(unique(cluster_summaries$Feature_Cluster_description))
    }
    if(PNG){
      png(paste0(fn, ".png"), width=100*(w+1), height=700*h)#, pointsize=5)
    }else{
      pdf(paste0(fn, ".pdf"), width=(1+CLUSTER_PROFILE_FREE_SCALE)*(w+1)+2, height=(1+CLUSTER_PROFILE_FREE_SCALE)*(h+1))
    }
    print(pcc)
    dev.off()
  }
  
  if(CLUSTER_AVERAGE_PLOTS || CLUSTER_AVERAGE_BOXPLOTS){
    
    cluster_averages<-do.call("rbind", cluster_averages)
    
    for (comp_gr in c(list(NULL), GROUP_AVG_COMPARISONS)){
      
      if(!is.null(comp_gr)){
        
        plot_subname<-paste(comp_gr, collapse="_vs_")
        h<-2*length(unique(cluster_averages$Mark))
        w<-2*length(unique(cluster_averages$Cluster))*length(unique(cluster_averages$Feature))
        
        if(CLUSTER_AVERAGE_BOXPLOTS){
          pac<-region_average_boxplot(
            cluster_averages, comp_groups=comp_gr, facet_vars=c("Cluster"))
          
          fn<-file.path(OUT_DIR, sprintf("cluster_average_boxplots_%s_%s_%d_%d_%s_subplot", 
                                         plot_subname, ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_")))
          if(PNG){
            png(paste0(fn, ".png"), width=100*(w+1), height=700*h)
          }else{
            pdf(paste0(fn, ".pdf"), width=1*(w+1), height=h+1)
          }
          print(pac)
          dev.off()
          
        }
      }
      
      if(CLUSTER_AVERAGE_PLOTS){
        
        for (comp_mk in c(list(NULL), MARK_AVG_COMPARISONS)){
          
          if(!(is.null(comp_gr) && is.null(comp_mk))){
            if(is.null(comp_gr)){
              plot_subname<-paste(comp_mk, collapse="_vs_")
              h<-2*length(sample_groups)
              w<-2*length(unique(cluster_averages$Cluster))*length(unique(cluster_averages$Feature))
            }else if(is.null(comp_mk)){
              plot_subname<-paste(comp_gr, collapse="_vs_")
              h<-2*length(unique(cluster_averages$Mark))
              w<-2*length(unique(cluster_averages$Cluster))*length(unique(cluster_averages$Feature))
            }else{
              plot_subname<-sprintf("%s_%s_vs_%s_%s", comp_mk[1], comp_gr[1], comp_mk[2], comp_gr[2])
              h<-2*length(unique(cluster_averages$Cluster))
              w<-2*length(unique(cluster_averages$Feature))
            }
            
            pac<-region_average_scatterplot(
              cluster_averages, comp_groups=comp_gr, comp_marks=comp_mk, 
              facet_vars=c("Cluster"), 
              scatter_fit=TRUE)
            
            fn<-file.path(OUT_DIR, sprintf("cluster_average_plots_%s_%s_%d_%d_%s_subplot", 
                                           plot_subname, ANALYSIS, OFFSET, WIN_SIZE, paste(plot_features, collapse="_")))
            if(PNG){
              png(paste0(fn, ".png"), width=100*(w+1), height=700*h)
            }else{
              pdf(paste0(fn, ".pdf"), width=1*(w+1), height=h+1)
            }
            print(pac)
            dev.off()
            
          }
        }
      }
    }
  }
}