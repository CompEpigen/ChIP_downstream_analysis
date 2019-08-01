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

baseCommand: ["plotCorrelation"]

arguments:
  - valueFrom: $(inputs.corData.nameroot + "_plotCorrelation.pdf")
    prefix: --plotFile
    position: 100

inputs:
  corData:
    type: File
    inputBinding:
      position: 1
      prefix: --corData
  
  corMethod:
    type: string                  #spearman or pearson
    inputBinding:
      position: 2
      prefix: --corMethod
  
  whatToPlot:
    type: string                  #heatmap or scatterplot
    inputBinding:
      position: 3
      prefix: --whatToPlot

outputs:
  output_plotCor:
    type: File
    outputBinding:
      glob: "*plotCorrelation.pdf"