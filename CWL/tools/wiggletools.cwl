cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: libddocker/wiggletools-1.2
    dockerOutputDirectory: /opt
  ShellCommandRequirement: {}    

baseCommand: ["wiggletools", "scale", "10000"]
arguments:
  - valueFrom: "|"
    position: 2
  - valueFrom: "sort"
    position: 3
  - valueFrom: "-k1,1"
    position: 4
  - valueFrom: "-k2,2n"
    position: 5
  - valueFrom: ">"
    position: 6
  - valueFrom: $(inputs.bw_file.basename).bg
    position: 7  

inputs:
  - id: bw_file
    type: File
    inputBinding:
      position: 1

outputs:
  - id: bg_file
    type: File
    outputBinding:
      glob: "*.bg"         