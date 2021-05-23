class: CommandLineTool
cwlVersion: v1.2
id: compute_GC_bias

requirements:
  - class: DockerRequirement
    dockerPull: 'kerstenbreuer/deeptools:3.1.1'
  - class: InlineJavascriptRequirement

baseCommand:
  - computeGCBias
  
inputs:
  - id: bam
    type: File
    secondaryFiles:
      - .bai
    inputBinding:
      position: 100
      prefix: '--bamfile'
  - id: effective_genome_size
    type: long
    inputBinding:
      position: 1
      prefix: '--effectiveGenomeSize'
  - id: genome
    type: File
    inputBinding:
      position: 2
      prefix: '--genome'
    doc: >-
      genome in 2/bit format. Download genome data here:
      http://hgdownload.cse.ucsc.edu/gbdb/
  - id: fragment_length
    type: int?
    inputBinding:
      position: 3
      prefix: '-l'
  - id: region
    type: string?
    inputBinding:
      position: 4
      prefix: '--region'
    doc: 
outputs:
  - id: GCbiasCheck
    doc: will be a txt file that can be used to normalize for GC bias if necessary
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + "_GCbias.txt")
  - id: GCplot
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + "_gc.png")
label: GCbiasCorrection
arguments:
  - position: 4
    prefix: '-freq'
    valueFrom: $(inputs.bam.nameroot + "_GCbias.txt")
  - position: 5
    prefix: '--biasPlot'
    valueFrom: $(inputs.bam.nameroot + "_gc.png")
  - position: 6
    prefix: --plotFileFormat
    valueFrom: 'png'
