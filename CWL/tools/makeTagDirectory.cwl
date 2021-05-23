class: CommandLineTool
cwlVersion: v1.0

requirements:
  - class: InitialWorkDirRequirement
    listing:
      - |
        ${
          return [
            {"class": "Directory",
             "basename": "default",
             "listing": [inputs.bam_file],
             "writable": true}
          ]
        }
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biowardrobe2/homer:v0.0.2'

baseCommand:
  - makeTagDirectory

inputs:
  - id: bam_file
    type: File
    doc: 'Alignment file, BAM'
  - id: fragment_size
    type:
      - int
      - string
      - 'null'
    inputBinding:
      position: 5
      prefix: '-fragLength'
    doc: |
      Set fragment size.
      By default is estimated as if it was single end ChIP-Seq experiment.
      Possible values:
        "#" - int value to be used as fragment size
        "given" - use read lengths
        "pe" - calculate from paired end read coordinates
  - id: max_length
    type: int?
    inputBinding:
      position: 8
      prefix: '-maxlen'
    doc: |
      Discard reads bigger then
  - id: min_length
    type: int?
    inputBinding:
      position: 7
      prefix: '-minlen'
    doc: |
      Discard reads smaller then
outputs:
  - id: output_tag_folder
    doc: Tag directory
    type: Directory
    outputBinding:
      glob: '$(inputs.bam_file.basename.split(''.'')[0])'
arguments:
  - position: 0
    valueFrom: '$(inputs.bam_file.basename.split(''.'')[0])'
  - position: 0
    valueFrom: $("default/" + inputs.bam_file.basename)