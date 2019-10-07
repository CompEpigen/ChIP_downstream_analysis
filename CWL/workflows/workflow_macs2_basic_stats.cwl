cwlVersion: v1.0
class: Workflow
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  treatment_bam: File[]
  control_bam: File[]
  genome_size: string
  broad: 
    type: boolean
    default: False
  qvalue: 
    type: float
    default: 0.05
  pvalue: float?                #if specified, MACS2 will use p in stead of q-value
  format_tag:
    type: string
    default: "AUTO"             #MACS2: "BAMPE" or "BEDPE" for paired-end data
  output:
    type: string
    default: "results.npz"
  corMethod:
    type: string
    default: "spearman"
  whatToPlot:
    type: string
    default: "heatmap"

outputs:
  output_npz:
    type: File
    outputSource: multiBamSummary/output_npz
  
  output_plotCor:
    type: File
    outputSource: plotCorrelation/output_plotCor

  output_plotPCA:
    type: File
    outputSource: plotPCA/output_plotPCA

steps:
  macs2:
    run: ../tools/MACS2.cwl
    scatter: [treatment_bam, control_bam]
    scatterMethod: "dotproduct"
    in:
      treatment_bam: treatment_bam
      control_bam: control_bam
      broad: broad
      format_tag: format_tag
      qvalue: qvalue
      pvalue: pvalue
      genome_size: genome_size
    out: [peak_bed]

  samtools_index:
    run: ../tools/samtools_index.cwl
    scatter: bam_sorted
    in:
      bam_sorted: 
        source: [treatment_bam, control_bam]
        linkMerge: merge_flattened
    out: [bam_sorted_indexed]
  
  multiBamSummary:
    run: ../tools/multiBamSummary.cwl
    in:
      bamfiles: [samtools_index/bam_sorted_indexed]
      output: output
    out: [output_npz]

  plotCorrelation:
    run: ../tools/plotCorrelation.cwl
    in:
      corData: multiBamSummary/output_npz
      corMethod: corMethod
      whatToPlot: whatToPlot
    out: [output_plotCor]
  
  plotPCA:
    run: ../tools/plotPCA.cwl
    in:
      corData: multiBamSummary/output_npz
    out: [output_plotPCA]