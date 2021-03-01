cwlVersion: v1.0
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerPull: zavolab/bedgraphtobigwig:4
    dockerOutputDirectory: /opt

baseCommand: ["bedGraphToBigWig"]
arguments:
  - valueFrom: $(inputs.bg_file.basename).bw
    position: 3  


inputs:
  - id: bg_file
    type: File
    inputBinding:
      position: 1
  - id: size_file
    type: File
    inputBinding:
      position: 2

outputs:
  - id: bw_file
    type: File
    outputBinding:
      glob: "*.bw"    