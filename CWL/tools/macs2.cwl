class: CommandLineTool
cwlVersion: v1.2

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
hints:
  - class: ResourceRequirement
    ramMin: 2000
    ramMax: 50000
    coresMin: 2
  - class: DockerRequirement
    dockerPull: 'genomicpariscentre/macs2:2.1.0.20140616'

baseCommand:
  - macs2
  - callpeak
inputs:
  - id: bam
    type: File
    inputBinding:
      position: 101
      prefix: '--treatment'
  - id: control
    type: File?
    inputBinding:
      position: 2
      prefix: '--control'
    doc: |
      The control or input file (e.g. IgG)
  - id: broad
    type: boolean?
    inputBinding:
      position: 3
      prefix: '--broad'
    doc: |
      MACS will call borad peaks. It is often suggested to turn off 
      model building and local lambda computation
  - id: genome_size
    type: long
    inputBinding:
      position: 3
      prefix: '--gsize'
    doc: |
      mappable genome size
  - default: true
    id: is_paired_end
    type: boolean
    doc: |
      if true,format is automatically set to BAMPE. This  will guess the 
      effective genome size by itself and will ignore --extsize
  - id: qvalue
    type: float?
    inputBinding:
      position: 3
      prefix: '--qvalue'
    doc: |
      minimum FDR cutoff
  - id: pvalue
    type: float?
    inputBinding:
      position: 3
      prefix: --pvalue
    doc:  |
      Pvalue cutoff for peak detection. DEFAULT: not set. In ENCODE ChIP-seq workflow: 0.01
  - id: nomodel
    type: boolean?
    inputBinding:
      position: 6
      prefix: '--nomodel'
    doc: |
      if true, will skip the model building where shift and extension size are calculated
      this is especially suggested if there is no input/control data
      If this is set true, arguments -shift and -extsize need to be set manually 
  - id: extsize
    type: int?
    inputBinding:
      position: 7
      prefix: '--extsize'
    doc: |
      if --nomodel flag is set: reads are extended by this value in 5'->3' direction
      in ENCODE ATAC seq: 150 bp is used
      in ENCODE ChIP seq: fragment length is used
  - id: shiftsize
    type: int?
    inputBinding:
      position: 2
      prefix: '--shift'
    doc: |
      if --nomodel flag is set: number of bp for shifting. 
      in ENCODE ATAC seq: $extsize/2
      in ENCODE ChIP seq: 0
  - id: broad_cutoff
    type: float?
    inputBinding:
      position: 2
      prefix: '--broad-cutoff'
outputs:
  - id: narrowpeak_file
    doc: outputs a bed file format containing the peak positions
    type: 'File[]'
    outputBinding:
      glob: '*Peak'
  - id: peaks_xls
    type: File
    outputBinding:
      glob: '*_peaks.xls'
  - id: bed_file
    type: File?
    outputBinding:
      glob: '*summits.bed'
arguments:
  - position: 1
    prefix: '--format'
    valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return "BAMPE";
        }
        else {
          return "BAM";
        }
      }
  - position: 2
    prefix: '--keep-dup'
    valueFrom: all
  - position: 100
    prefix: '--name'
    valueFrom: $(inputs.bam.nameroot + "_macs2")
