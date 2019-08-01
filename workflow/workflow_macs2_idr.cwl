cwlVersion: v1.0
class: Workflow
requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  treatment_bam_1: File
  treatment_bam_2: File
  control_bam_1: File
  control_bam_2: File
  genome_size: string
  broad: 
    type: boolean
    default: False
  qvalue: float?
  pvalue: float?
  tagformat: string?           #MACS2: "BAMPE" or "BEDPE" for paired-end data
  filetype: string?            #IDR: narrowPeak, broadPeak or bed
  plot: boolean?             #IDR
  output_basename: string?

outputs:
  idr_values:
    type: File
    outputSource: idr/peaks_idr_scores
  std_error:
    type: File
    outputSource: idr/stderr_out

steps:
  macs2_1:
    run: /Users/adams/Documents/Heidelberg/CWL/tools/MACS2_1.cwl
    in:
      treatment_bam_1: treatment_bam_1
      control_bam_1: control_bam_1
      broad: broad
      format_tag: tagformat
      qvalue: qvalue
      pvalue: pvalue
      genome_size: genome_size
    out:
      [peak1_bed]

  macs2_2:
    run: /Users/adams/Documents/Heidelberg/CWL/tools/MACS2_2.cwl
    in:
      treatment_bam_2: treatment_bam_2
      control_bam_2: control_bam_2
      broad: broad
      format_tag: tagformat
      qvalue: qvalue
      pvalue: pvalue
      genome_size: genome_size
    out:
      [peak2_bed]

  idr:
    run: /Users/adams/Documents/Heidelberg/CWL/tools/IDR.cwl
    in:
      samples:
        - macs2_1/peak1_bed
        - macs2_2/peak2_bed
      inputFileType: filetype
      output_basename: output_basename
    out:
      - peaks_idr_scores
      - stderr_out               #or [peaks_idr_scores, stderr_out]