cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}       #stderr
  StepInputExpressionRequirement: {}    #for ValueForm
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: genomicpariscentre/macs2:2.1.0.20140616

baseCommand: ["macs2", "callpeak"]

arguments:
  - valueFrom: "--nomodel"
    position: 3
  - valueFrom: "all"
    prefix: "--keep-dup"
    position: 2
  - valueFrom: $(inputs.treatment_bam.nameroot)
    prefix: "--name"
    position: 100

inputs:
  treatment_bam:
    type: File
    inputBinding:
      position: 1
      prefix: "--treatment"

  control_bam:
    type: File
    inputBinding:
      position: 2
      prefix: "--control"

  genome_size:
    doc: can be "mm", "hs", "ce", "dm", or the total number of genomic bp 
    type: string
    inputBinding:
      position: 3
      prefix: "--gsize"

  format_tag:
    doc: can be "ELAND", "BED", "ELANDMULTI", "ELANDEXPORT", "ELANDMULTIPET", "SAM", "BAM", "BOWTIE", "BAMPE" or "BEDPE"
    #if the format is specified as "BAMPE" or "BEDPE", MACS2 will process the files as paired-end data
    type: string
    inputBinding:
      position: 5
      prefix: "--format"

  broad:
    type: boolean
    inputBinding:
      position: 4
      prefix: "--broad"    #if you add it in the command line its true, else false

  qvalue:
    type: float
    inputBinding:
      position: 6
      prefix: "--qvalue"   #minimum q-value

  pvalue:
    type: float?
    inputBinding:
      position: 7
      prefix: "--pvalue"   #if specified, MACS2 will use p in stead of q-value


outputs:
  peak_bed:    
    type: File
    outputBinding:
      glob: "*Peak"

  peak_xls:
    type: File
    outputBinding:
      glob: "*_peaks.xls"