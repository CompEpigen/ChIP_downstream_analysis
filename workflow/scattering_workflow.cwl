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
  #output_basename: string

outputs:
  peak_bed:
    type: File[]
    outputSource: macs2/peak_bed


steps:
  macs2:
    run: /Users/adams/Documents/Heidelberg/CWL/tools/MACS2.cwl
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