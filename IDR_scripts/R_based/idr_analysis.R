CALLER<-"macs2"
INPUT=c("wo_input", "h3_input", "AO_input")[1]
STAT=c("pval","fdr")[2]
CUTOFF=c("0.01","0.05","0.1","0.2")[1]

### loading peak files
if(CALLER=="peakranger"){
    #token="_peaks_region"
    #token="_ranger_narrow_peaks_region"
    token="_bcp_broad_peaks_region"
    pf_extension<-rep("bed", length(marks))
    names(pf_extension)<-marks
    signal_column<-5
    
}else if(CALLER=="macs2"){
    token="_peaks"
    pf_extension<-c(rep("narrowPeak", 3), rep("gappedPeak", 4))
    names(pf_extension)<-marks
    #pf_extension<-rep("bed", length(marks))
    #pf_extension<-"gappedPeak"
    #pf_extension<-"broadPeak"
    signal_column<-5
    
}

PEAK_OUT_DIR<-sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/peaks/%s/%s/%s%s/",CALLER, INPUT, STAT, CUTOFF)


#### idr analysis using the general-purpose implementation

library(idr)

#mu <- 5
#sigma <- 1.3
#rho <- 0.5
#p <- 0.5

#
#qval_mat<- -log10(FDR-FDR*signal_matrices[["AR_vs_AS_K27me3"]]/1000)
#
#
#qval_mat<- -log10(signal_matrices[["AR_vs_AS_K27me3"]])

signal_marices<-readRDS(PEAK_OUT_DIR,"signal_matrices.RDS")

###### estimation

mu <- 1
sigma <- 1
rho <- 0.1
p <- 0.5

idr.res<-list()
for(pfn in names(signal_matrices)){
    
    qval_mat<- signal_matrices[[pfn]]
    
    if(CALLER %in% c("peakranger")){
        qval_mat<- -log10(qval_mat)
    }
    
    #idr.out <- est.IDR(qval_mat, mu, sigma, rho, p, eps=0.0001, max.ite=20)
    #idr.res[[pfn]]<-idr.out
    
    
    rank.x <- rank(qval_mat[,1])
    rank.y <- rank(qval_mat[,2])
    uv <- get.correspondence(rank.x, rank.y, seq(0.01, 0.99, by=1/200))
    # plot correspondence curve on the scale of percentage
    
    pdf(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/H33/%s_%s_%s.pdf",CALLER,token,pfn))
    
    plot(uv$psi$t, uv$psi$value, xlab="t", ylab="psi", xlim=c(0, max(uv$psi$t)),
            ylim=c(0, max(uv$psi$value)), cex.lab=2)
    lines(uv$psi$smoothed.line, lwd=4)
    abline(coef=c(0,1), lty=3)
    
    # plot the derivative of correspondence curve on the scale of percentage
    plot(uv$dpsi$t, uv$dpsi$value, xlab="t", ylab="psi ", xlim=c(0, max(uv$dpsi$t)),
            ylim=c(0, max(uv$dpsi$value)), cex.lab=2)
    lines(uv$dpsi$smoothed.line, lwd=4)
    abline(h=1, lty=3)
    # plot correspondence curve on the scale of the number of observations
    plot(uv$psi.n$t, uv$psi.n$value, xlab="t", ylab="psi", xlim=c(0, max(uv$psi.n$t)),
            ylim=c(0, max(uv$psi.n$value)), cex.lab=2)
    lines(uv$psi.n$smoothed.line, lwd=4)
    abline(coef=c(0,1), lty=3)
    
    
    # plot the derivative of correspondence curve on the scale of the number
    # of observations
    plot(uv$dpsi.n$t, uv$dpsi.n$value, xlab="t", ylab="psi ", xlim=c(0, max(uv$dpsi.n$t)),
            ylim=c(0, max(uv$dpsi.n$value)), cex.lab=2)
    lines(uv$dpsi.n$smoothed.line, lwd=4)
    abline(h=1, lty=3)
    # If the rank lists consist of a large number of ties at the bottom
    # (e.g. the same low value is imputed to the list for the observations
    # that appear on only one list), it may be desirable to plot only
    # observations before hitting the ties. Then it can be plotted using the
    # following
    plot(uv$psi$t[1:uv$psi$jump.point], uv$psi$value[1:uv$psi$jump.point], xlab="t",
            ylab="psi", xlim=c(0, max(uv$psi$t[1:uv$psi$jump.point])),
            ylim=c(0, max(uv$psi$value[1:uv$psi$jump.point])), cex.lab=2)
    lines(uv$psi$smoothed.line, lwd=4)
    abline(coef=c(0,1), lty=3)

    dev.off()

}
###### IDR analysis using the ENCODE scripts
#[half.width]: Set this to -1 (NULL) if you want to use the reported peak width in the peak files.
#[overlap.ratio]: fractional bp overlap (ranges from 0 to 1) between peaks in replicates to be considered as overlapping peaks. IMPORTANT: This parameter has not been tested fully. It is recommended to set this to 0.
#[is.broadpeak]: Is the peak file format narrowPeak or broadPeak. Set to F if it is narrowPeak/regionPeak or T if it is broadPeak.
#[ranking.measure] is the ranking measure to use. It can take only one of the following values signal.value, p.value or q.value.

source("/ngs_share/tools/idr_ENCODE/R/functions_idr_ENCODE.R")
half.width <- NULL
overlap.ratio <- 0
#sig.value <- "p.value"
sig.value <- "signal.value"


genome<-read.table("/ngs_share/data/genomes/Hsapiens/hg19/seq/hg19.chrom.sizes")
colnames(genome)<-c("UCSC_seqlevel","UCSC_seqlength")
idx <- which(!grepl("random|hap|chrUn", genome$UCSC_seqlevel))
genome <- genome[idx,1:2]
chr.size<-genome

require(GenomicRanges)
replicate_pairs<-list(c("AR", "AS"), c("AS", "AT"), c("AR", "AT"), c("GCTB11", "KM882"), c("AS","AO"), c("AR","AO"))#[c(1,5:6)]
dir.create(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/idr_ENCODE/%s%s/%s/%s%s/",
                            CALLER,token,INPUT,STAT,CUTOFF), recursive=TRUE)
for(mark in marks){
    for(rp in replicate_pairs){
        
        if(grepl("h3_input", BED_DIR)){
            out_bed_file1<-file.path(BED_DIR, sprintf("%s_%s_vs_%s_H3%s.%s", rp[1], mark, rp[1], token, pf_extension[mark]))
            out_bed_file2<-file.path(BED_DIR, sprintf("%s_%s_vs_%s_H3%s.%s", rp[2], mark, rp[2], token, pf_extension[mark]))
            
        }else{
            out_bed_file1<-file.path(BED_DIR, sprintf("%s_%s%s.%s", rp[1], mark, token, pf_extension[mark]))
            out_bed_file2<-file.path(BED_DIR, sprintf("%s_%s%s.%s", rp[2], mark, token, pf_extension[mark]))
        }
        
       
       print(mark)
       print(rp)
       if(file.exists(out_bed_file1) && file.exists(out_bed_file2)){
            
            #            print(sprintf("reading %s",out_bed_file1 ))
            #            
            #            bf<-read.table(out_bed_file, sep="\t", comment.char="#")
            #            print(dim(bf))
            #            bf<-bf[bf$V3-bf$V2>0,]
            ######## code from http://ccg.vital-it.ch/var/sib_april15/cases/landt12/idr.html#hide2
            is.broadpeak <- mark %in% marks[4:7]
            is.broadpeak<-FALSE
            is.gappedpeak <- mark %in% marks[4:7]
            #is.gappedpeak<-FALSE
            is.peakranger<-CALLER %in% c("peakranger")
            
            rep1 <- process.narrowpeak(paste(out_bed_file1, sep=""), chr.size, 
                    half.width=half.width, summit="offset", broadpeak=is.broadpeak, gappedpeak=is.gappedpeak, peakranger=is.peakranger)
            rep2 <- process.narrowpeak(paste(out_bed_file2, sep=""), chr.size, 
                    half.width=half.width, summit="offset", broadpeak=is.broadpeak, gappedpeak=is.gappedpeak, peakranger=is.peakranger)
            sapply(rep1, function(rr) print(head(rr)))
            sapply(rep2, function(rr) print(head(rr)))
            
            uri.output <- compute.pair.uri(rep1$data.cleaned, rep2$data.cleaned, 
                    sig.value1=sig.value, sig.value2=sig.value, overlap.ratio=overlap.ratio)
            em.output <- fit.em(uri.output$data12.enrich, fix.rho2=T)
            idr.local <- 1-em.output$em.fit$e.z
            IDR <- c()
            o <- order(idr.local)
            IDR[o] <- cumsum(idr.local[o])/c(1:length(o))
            idr_output <- data.frame(chr1=em.output$data.pruned$sample1[, "chr"],
                    start1=em.output$data.pruned$sample1[, "start.ori"],
                    stop1=em.output$data.pruned$sample1[, "stop.ori"],
                    sig.value1=em.output$data.pruned$sample1[, "sig.value"],   
                    chr2=em.output$data.pruned$sample2[, "chr"],
                    start2=em.output$data.pruned$sample2[, "start.ori"],
                    stop2=em.output$data.pruned$sample2[, "stop.ori"],
                    sig.value2=em.output$data.pruned$sample2[, "sig.value"],
                    idr.local=1-em.output$em.fit$e.z, IDR=IDR)
            
            write.table(idr_output, 
                    sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/idr_ENCODE/%s%s/%s/%s%s/idr_overlapped_peaks_%s%s_%s_%s_vs_%s_%s.txt",
                            CALLER,token,INPUT,STAT,CUTOFF,CALLER,token,mark,rp[1],rp[2],sig.value), sep="\t", quote=F)
            
            filtered_peaks <- idr_output[idr_output[,10]<=0.01,]
            dim(filtered_peaks) # get the number of peaks
            
            ez.list <- get.ez.tt.all(em.output, uri.output$data12.enrich$merge1, uri.output$data12.enrich$merge2)
            
            pdf(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/idr_ENCODE/%s%s/%s/%s%s/idr_analysis_%s%s_%s_%s_vs_%s_%s.pdf",
                            CALLER,token,INPUT,STAT,CUTOFF,CALLER,token,mark,rp[1],rp[2],sig.value), width=8, height=3)
                par(mar=c(5,5,0,0.5), mfrow = c(1,3), oma=c(5,0,2,0))
                idr_output$col[idr_output[,10]<=0.01]="black"
                idr_output$col[idr_output[,10]>=0.01]="red"
                plot(log(idr_output[,4]),log(idr_output[,8]),col=idr_output[,11], pch=19, 
                        xlab="log(signal) Rep1", ylab="log(signal) Rep2")
                legend("topleft", c("IDR=>0.01","IDR<=0.01"), col=c("red","black"), pch=19, 
                        bty="n", lty=c(1,1), lwd=c(2,2))
                plot(rank(-idr_output[,4]),rank(-idr_output[,8]),col=idr_output[,11], pch=19, 
                        xlab="Peak rank Rep1", ylab="Peak rank Rep2")
                legend("topleft", c("IDR=>0.01","IDR<=0.01"), col=c("red","black"), pch=19, 
                        bty="n", lty=c(1,1), lwd=c(1,1))
                plot(ez.list$IDR, ylab="IDR", xlab="num of significant peaks")
            dev.off()
            
        }else{
            print(sprintf("%s does not exist",out_bed_file1))
        }
    }
}

################################################
######################### IDR processing
#################################################

IDR_RES_DIR<-"/C010-datasets/Internal/GCTB/ChIP/idr"
#IDR_RES_DIR<-"/C010-datasets/Internal/GCTB/ChIP/idr_pval_e-4"


marks<-c("H3", "H3.3", "G34W")

IDR_THRESHOLD<-0.05
library(GenomicRanges)
gr_all<-gr_passing<-list()
total_peaks<-passing_peaks<-list()
for (mark in marks[7]){
    for(rp in replicate_pairs[1]){
    #idr_tab<-read.table(file.path(IDR_RES_DIR, sprintf("AS_%s_vs_AR_%s_idr_values.txt", mark, mark)),
    #        sep="\t", header=FALSE)
    
    idr_tab<-read.table(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/idr_ENCODE/macs2/wo_input/stringent/gapped/idr_overlapped_peaks_%s%s_%s_%s_vs_%s.txt",
            CALLER,token,mark,rp[1],rp[2]),sep="\t", header=TRUE)
    
    total_peaks[[mark]]<-nrow(idr_tab)
    #gr_all[[mark]]<-GRanges(seqnames=idr_tab$V1, IRanges(idr_tab$V17, idr_tab$V18), strand="*")
    gr_all[[mark]]<-GRanges(seqnames=idr_tab$chr1, IRanges(idr_tab$start1, idr_tab$stop1), strand="*")
        
   # passing<-idr_tab[which(idr_tab$V5>as.integer(-125*log2(IDR_THRESHOLD))),]
     passing<-idr_tab[which(idr_tab$IDR<IDR_THRESHOLD),]
        
    passing_peaks[[mark]]<-nrow(passing)
    #gr_passing[[mark]]<-GRanges(seqnames=passing$V1, IRanges(passing$V17, passing$V18), strand="*")
    gr_passing[[mark]]<-GRanges(seqnames=passing$chr1, IRanges(passing$start1, passing$stop1), strand="*")
    
    rtracklayer::export.bed(gr_passing[[mark]], 
            sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/%s_IDR_peaks.bed", mark))
    
}
}


OFFSET<-500

olaps_h3_h33<-findOverlaps(gr_passing[["H3.3"]]+OFFSET, gr_passing[["H3"]]+OFFSET)

olaps_h33_g34w<-findOverlaps(gr_passing[["G34W"]]+OFFSET, gr_passing[["H3.3"]]+OFFSET)

olaps_h3_g34w<-findOverlaps(gr_passing[["G34W"]]+OFFSET, gr_passing[["H3"]]+OFFSET)


library(ChIPPeakAnno)

pdf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/idr_venn.pdf")
makeVennDiagram(gr_passing, 
        NameOfPeaks=marks, ignore.strand=TRUE)
dev.off()

h3_h33_common<-reduce(c(gr_passing[["H3"]],gr_passing[["H3.3"]]))
g34w_unique<-setdiff(gr_passing[["G34W"]], h3_h33_common)
g34w_unique<-g34w_unique[width(g34w_unique)>100]


rtracklayer::export.bed(g34w_unique, "/ngs_share/scratch/pavlo/gctb/analysis/ChIP/H33/g34w_unique_peaks.bed")


library(EnsDb.Hsapiens.v75) ##(hg19)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

for(mark in marks[7]){

pa<-annotatePeak(gr_passing[[mark]], 
        tssRegion = c(-3000, 3000), 
        TxDb = NULL, level = "transcript", 
        assignGenomicAnnotation = TRUE, 
        genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
        annoDb = NULL, addFlankGeneInfo = FALSE, 
        flankDistance = 5000, sameStrand = FALSE, 
        ignoreOverlap = FALSE, ignoreUpstream = FALSE, ignoreDownstream = FALSE, overlap = "TSS", verbose = TRUE)

pdf(sprintf("/ngs_share/scratch/pavlo/gctb/analysis/ChIP/IDR/%s_passing_unique_annotation_pie.pdf", mark))
plotAnnoPie(pa)
dev.off()

}




