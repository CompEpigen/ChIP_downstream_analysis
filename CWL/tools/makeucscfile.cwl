#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - makeUCSCfile
inputs:
  - id: tag_directory
    type: Directory
    doc: 'HOMER specific tag directory created by makeTagDir'
    inputBinding:
      position: 0
  - id: style
    type:
      - string?
    inputBinding:
      position: 5
      prefix: '-style'
    doc: |
      what specific analysis is done:
      chip-seq experiment is default
  - id: fragment_length
    type:
      - "null"
      - int
      - string
    inputBinding:
      position: 6
      prefix: "-fragLength"
    doc: |
      Set fragment size.
      By default is estimated from autocorrelation
      Possible values:
        "#" - int value to be used as fragment size
        "given" - use read lengths
  - id: resolution
    type: int?
    inputBinding:
      position: 7
      prefix: '-res'
    doc: |
      overall resolution in bp, default is 1
      can be useful to reduce UCSC file size
  - id: bigwig creation
    type: File?
    inputBinding:
      position: 8
      prefix: '-bigWig'
    doc: |
      This option requires bedGraphToBigWig script installation and access
      as file give the chromosome site as input

arguments: 
  - valueFrom: $(inputs.tag_directory.basename).bedGraph
    prefix: -o
    position: 3

outputs:
  - id: bedGraph
    doc: visualization file of the coverage track
    type: File
    outputBinding:
      glob: "*.bedGraph*"

requirements:     
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biowardrobe2/homer:v0.0.2'
