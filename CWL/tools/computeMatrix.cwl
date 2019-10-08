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
  scale_regions_or_use_reference_point:
    type: boolean
    default: false
    valueFrom: |
      ${
        if(self){
          return("scale-regions")
        }
        else{
          return("reference-point")
        }
      }
    inputBinding:
      position: 1
  
  before_region_start_length:
    type: int?
    inputBinding:
      position: 2
      prefix: "--beforeRegionStartLength"

  after_region_start_length:
    type: int?
    inputBinding:
      position: 2
      prefix: "--afterRegionStartLength"

  reference_point:
    type: string?
    inputBinding:
      position: 2
      prefix: "referencePoint"

  regions_bed:
    type: File[]
    inputBinding:
      position: 10
      prefix: --regionsFileName

  scoreFileName:                                 #bigWigfiles
    type: File[]
    inputBinding:
      position: 10
      prefix: --scoreFileName

  outFileName:
    type: string
    default: "computeMatrix.gz"
    inputBinding:
      position: 10
      prefix: --outFileName

outputs:
  matrix_gzip:
    type: File
    outputBinding:
      glob: "*Matrix.gz"