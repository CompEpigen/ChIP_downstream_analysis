cwlVersion: v1.0
class: CommandLineTool

requirements:
   - class: InlineJavascriptRequirement
hints:
   - class: DockerRequirement
     dockerPull: 'biowardrobe2/homer:v0.0.2'

baseCommand:
  - findPeaks

inputs:
  - id: tag_directory
    type: Directory
    inputBinding:
      position: 0
    doc: HOMER specific tag directory created by makeTagDir
  - id: input_tag_directory
    type: 
      - 'null'
      - Directory
    inputBinding:
      position: 2
      prefix: -i
    doc: Background signal must also be covnerted to a Tag Directory
  - id: style
    type: string
    inputBinding:
      position: 5
      prefix: '-style'
    doc: |
      what specific analysis is done:
      histone, factor, groseq, tss dnase, super(super enhancers) 
  - id: size
    type: int?
    inputBinding:
      position: 8
      prefix: '-size'
    doc: |
      set the peak width
      for histone mode: default: 500
      for "factor" mode: size is calculared from autocorrelation analysis
  - id: minDist_between_peaks
    type: int?
    inputBinding:
      position: 7
      prefix: '-minDist'
    doc: |
      Discard reads smaller then

arguments:
  - valueFrom: $(inputs.tag_directory.basename + "_HOMERpeaks.txt")
    prefix: -o
    position: 3


outputs:
  - id: peaks
    doc: |
      peak informations will be outputted to peaks.txt file.
      In contrast to the "-o auto" argumment in HOMER software, 
      the name of the output file will be always peaks.txt, no matter if style is factor or histone
    type: File
    outputBinding:
      glob: $(inputs.tag_directory.basename + "_HOMERpeaks.txt")

