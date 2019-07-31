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

baseCommand: ["plotProfile"]

inputs:
  matrixFile:
    type: File 
    inputBinding:
      position: 1
      prefix: --matrixFile
  
  outFileName:
    type: string
    default: "Profile_plot.pdf"
    inputBinding:
      position: 2
      prefix: --outFileName

outputs:
  plotProfile_pdf:
    type: File
    outputBinding:
      glob: "*plot.pdf"
