#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
inputs:
  - id: bw_dir
    type: Directory
  - id: size_file
    type: File    
  - id: output_dir_name
    type: string

outputs: 
#  bw_file:
#    type: File[]
#    outputSource: bedGraphToBigWig/bw_file
#  bg_file:
#    type: File[]
#    outputSource: wiggletools/bg_file
  out_dir:
    type: Directory
    outputSource: output_to_dir/out_dir    

steps:
  - id: files_in_dir
    run:
      class: ExpressionTool
      requirements: { InlineJavascriptRequirement: {} }
      inputs:
        dir: Directory
      expression: '${return {"files": inputs.dir.listing};}'
      outputs:
        files: File[]
    in:
      dir: bw_dir
    out: [files]

  - id: wiggletools
    run: "../tools/wiggletools.cwl"
    in:
      - id: bw_file
        source: files_in_dir/files
    scatter: bw_file
    out:
      - id: bg_file

  - id: bedGraphToBigWig
    run: "../tools/bedGraphToBigWig.cwl"
    in:
      - id: bg_file
        source: wiggletools/bg_file
      - id: size_file
        source: size_file
    out:
      - id: bw_file
    scatter: bg_file

  - id: output_to_dir
    run:
      class: ExpressionTool
      requirements: { InlineJavascriptRequirement: {} }
      inputs:
        fs: File[]
        output_dir_name: string
      expression: '${
                    return {"out_dir": {
                            "class": "Directory", 
                            "basename": inputs.output_dir_name,
                            "listing": inputs.fs
                            }};}'
      outputs:
        out_dir: Directory
    in:
      - id: fs
        source: bedGraphToBigWig/bw_file
      - id: output_dir_name
        source: output_dir_name
    out: 
      - id: out_dir

     



