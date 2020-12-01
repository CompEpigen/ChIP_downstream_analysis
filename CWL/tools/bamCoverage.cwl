cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}       #stderr
  StepInputExpressionRequirement: {}    #for valueFrom

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["bamCoverage"]
arguments:
  - valueFrom: $(inputs.bam.nameroot).bigwig
    prefix: --outFileName
    position: 100

  - valueFrom: "bigwig"
    prefix: --outFileFormat
    position: 101

inputs:
  bam:
    type: File
    inputBinding:
      position: 1
      prefix: --bam

  binSize:
    type: int
    inputBinding:
      position: 2
      prefix: --binSize
  
  ignoreForNormalization:
    type: string?
    inputBinding:
      position: 3
      prefix: --ignoreForNormalization

  fragmentLength:
    type: int
    inputBinding:
      position: 4
      prefix: --extendReads

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: "*filt.bigwig"
