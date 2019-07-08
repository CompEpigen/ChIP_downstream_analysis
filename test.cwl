#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["samtools", "view"] 
inputs:
  filename:
    type: string
    inputBinding:
      position: 1
outputs: []
