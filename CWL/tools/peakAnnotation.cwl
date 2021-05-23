class: CommandLineTool
cwlVersion: v1.0
id: peakAnnotation

requirements:
  - class: InlineJavascriptRequirement
hints:
  - class: DockerRequirement
    dockerPull: 'biowardrobe2/homer:v0.0.2'

baseCommand:
  - annotatePeaks.pl
inputs:
  - id: peak_txt_file
    type: File
    doc: the peak file is in a HOMER specific txt file
    inputBinding:
      position: 0
  - id: genome
    type: File
    doc: reference genome where the peaks will be mapped to
    inputBinding:
      position: 1
  - id: gtf_file
    type: File
    inputBinding:
      position: 2
      prefix: '-gtf'

stdout: $(inputs.peak_txt_file.basename)_annotated.txt

outputs:
  - id: annotated_peaks
    doc: | 
      will output a bed file with peaks that can be 
      visualized with a Genome Browser
    type: stdout