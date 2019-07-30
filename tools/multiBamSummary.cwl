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

baseCommand: ["multiBamSummary", "bins"]

arguments:
  - valueFrom: $(inputs.bamfiles.nameroot + inputs.output)
    prefix: --outFileName
    position: 100

inputs:
  bamfiles:
    type: File[]
    secondaryFiles: .bai
    inputBinding:
      position: 1
      prefix: --bamfiles
  
  output:
    type: string

outputs:
  output_npz:
    type: File
    outputBinding:
      glob: "*.npz"