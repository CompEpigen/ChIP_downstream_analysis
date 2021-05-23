class: CommandLineTool
cwlVersion: v1.0
id: pos2bed

requirements:
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biowardrobe2/homer:v0.0.2'

baseCommand:
  - pos2bed.pl
inputs:
  - id: peak_txt_file
    type: File
    doc: the peak file is in a HOMER specific txt file
    inputBinding:
      position: 0

stdout: $(inputs.peak_txt_file.nameroot).bed

outputs:
  - id: BED_file
    doc: | 
      will output a bed file with peaks that can be 
      visualized with a Genome Browser
    type: stdout
