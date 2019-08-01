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

baseCommand: ["plotPCA"]

arguments:
  - valueFrom: $(inputs.corData.nameroot + "_plotPCA.pdf")
    prefix: --plotFile
    position: 100

inputs:
  corData:
    type: File
    inputBinding:
      position: 1
      prefix: --corData
  
outputs:
  output_plotPCA:
    type: File
    outputBinding:
      glob: "*plotPCA.pdf"