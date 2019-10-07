cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}       #stderr
  StepInputExpressionRequirement: {}    #for valueFrom

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["computeMatrix", "scale-regions"]

inputs:
  regions_bed:
    type: File[]
    inputBinding:
      position: 1
      prefix: --regionsFileName

  scoreFileName:                                 #bigWigfiles
    type: File[]
    inputBinding:
      position: 2
      prefix: --scoreFileName

  outFileName:
    type: string
    default: "computeMatrix.gz"
    inputBinding:
      position: 3
      prefix: --outFileName

outputs:
  matrix_gzip:
    type: File
    outputBinding:
      glob: "*Matrix.gz"