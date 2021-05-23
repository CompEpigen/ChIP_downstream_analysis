class: CommandLineTool
cwlVersion: v1.0
id: gcbiascorrection

requirements:
  - class: DockerRequirement
    dockerPull: 'kerstenbreuer/deeptools:3.1.1'
  - class: InlineJavascriptRequirement

baseCommand:
  - correctGCBias
inputs:
  - id: bam
    type: File
    inputBinding:
      position: 100
      prefix: '--bamfile'
    secondaryFiles:
      - .bai
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
  - id: GCbias_file
    type: File
    inputBinding:
      position: 3
      prefix: '--GCbiasFrequenciesFile'
  - id: region
    type: string?
    inputBinding:
      position: 4
      prefix: '--region' 
outputs:
  - id: GCcorrectedBam
    doc: will be a txt file that can be used to normalize for GC bias if necessary
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + "_GCcorrected.bam")
label: GCbiasCorrection
arguments:
  - position: 4
    prefix: '-o'
    valueFrom: $(inputs.bam.nameroot + "_GCcorrected.bam")
