cwlVersion: v1.0
class: CommandLineTool
requirements:
 InlineJavascriptRequirement: {}       #stderr
 StepInputExpressionRequirement: {}    #for ValueForm

baseCommand: ["idr"]

#arguments:
  #- valueFrom:  $(inputs.samples.nameroot + "_idrValue.txt")
    #prefix: "--ouput-file"
    #position: 4

inputs:
 samples:
  type: File[]
  inputBinding:
   position: 1
   prefix: --samples

 inputFileType:
  type: string
  default: "narrowPeak"
  inputBinding:
   position: 2
   prefix: --input-file-type

 plot:
  type: boolean?
  inputBinding:
   position: 3
   prefix: --plot

 outputFile:
  type: string
  default: "idrValue.txt"
  inputBinding:
   position: 4
   prefix: --output-file

 idrThreshold:
  type: float?                      
  inputBinding:
   position: 5
   prefix: --idr-threshold
    
#stderr: stderr_output.txt           #capture standard error

outputs:
 stderr_out:
  type: stderr

 peaks_idr_scores:
  type: File
  outputBinding:
   glob: $(inputs.outputFile)