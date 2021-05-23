#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'

requirements:
   - class: InlineJavascriptRequirement
hints:
   - class: DockerRequirement
     dockerPull: 'quay.io/biocontainers/ucsc-bedgraphtobigwig:377--h0b8a92a_2'

baseCommand:
  - bedGraphToBigWig
inputs:
  - id: bedGraph
    type: File
    inputBinding:
      position: 0
    doc: bedgraphs are larger than compressed bigwig files
  - id: genome
    type: File
    inputBinding:
      position: 1
    doc: reference genome

arguments:
  - valueFrom: $(inputs.bedGraph.nameroot).bw
    position: 2


outputs:
  - id: bigwig
    type: File
    outputBinding:
      glob: $(inputs.bedGraph.nameroot).bw
