cwlVersion: v1.0
class: Workflow
requirements:
 MultipleInputFeatureRequirement: {}

inputs:
 treatment_bam_1: File
 treatment_bam_2: File
 control1: File
 control2: File
 gsize: string
 isbroad: 
  type: boolean
  default: False
 threshq: float?
 threshp: float?
 tagformat: string?           #MACS2: "BAMPE" or "BEDPE" for paired-end data
 filetype: string?            #IDR: narrowPeak, broadPeak or bed
 isplot: boolean?             #IDR
 outp: string?

outputs:
 idr_values:
  type: File
  outputSource: idr/peaks_idr_scores
 std_error:
  type: File
  outputSource: idr/stderr_out

steps:
 macs2_1:
  run: MACS2_1.cwl
  in:
   treatment_bam_1: treatment_bam_1
   control_bam_1: control1
   broad: isbroad
   format_tag: tagformat
   qvalue: threshq
   pvalue: threshp
   genome_size: gsize
  out:
   [peak1_bed]

 macs2_2:
  run: MACS2_2.cwl
  in:
   treatment_bam_2: treatment_bam_2
   control_bam_2: control2
   broad: isbroad
   format_tag: tagformat
   qvalue: threshq
   pvalue: threshp
   genome_size: gsize
  out:
   [peak2_bed]

 idr:
  run: IDR.cwl
  in:
   samples:
    - macs2_1/peak1_bed
    - macs2_2/peak2_bed
   inputFileType: filetype
   outputFile: outp
  out:
   - peaks_idr_scores
   - stderr_out               #or [peaks_idr_scores, stderr_out]

