cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: libddocker/wiggletools-1.2
    dockerOutputDirectory: /opt
  ShellCommandRequirement: {}    

baseCommand: ["wiggletools", "scale"]
arguments:
  - valueFrom: "|"
    position: 3
  - valueFrom: "sort"
    position: 4
  - valueFrom: "-k1,1"
    position: 5
  - valueFrom: "-k2,2n"
    position: 6
  - valueFrom: ">"
    position: 7
  - valueFrom: $(inputs.bw_file.basename).bg
    position: 8  

inputs:
  - id: bw_file
    type: File
    inputBinding:
      position: 2
  - id: scaling_factor
    type: int
    inputBinding:
      position: 1

outputs:
  - id: bg_file
    type: File
    outputBinding:
      glob: "*.bg"         