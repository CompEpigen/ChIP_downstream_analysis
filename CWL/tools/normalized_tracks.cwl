cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
hints:
  - class: ResourceRequirement
    ramMin: 20000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'kerstenbreuer/deeptools:3.1.1'

baseCommand:
  - bamCoverage

inputs:
  - id: bam
    type: File
    inputBinding:
      position: 100
      prefix: '--bam'
    doc: bam file as input; needs bai index file in the same directory
    secondaryFiles:
      - .bai
  - id: effective_genome_size
    type: long
    inputBinding:
      position: 10
      prefix: '--effectiveGenomeSize'
    doc: >
      the effectively mappable genome size
      see: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
  - id: fragment_size
    type: int?
    doc: mean library fragment size; used to extend the reads
  - id: ignoreForNormalization
    type: string?
    inputBinding:
      position: 10
      prefix: '--ignoreForNormalization'
    default: chrX chrY chrM
    doc: |
      List of space-delimited chromosome names that shall be ignored 
      when calculating the scaling factor. 
      Specify as space-delimited string. 
      Default: "chrX chrY chrM"
  - id: is_paired_end
    type: boolean
    doc: 'if false, reads are extended by fragment_size'
  - id: spike_in_count
    type: long?
    inputBinding:
      position: 10
      prefix: '--scaleFactor'
      valueFrom: |
        ${
          if( self == null ){
            return null
          }
          else{
            return (1.0 / parseFloat(self)) 
          }
        }
    doc: 'number of reads aligned to the spike in genome, optional'
  - id: normalization
    type: string?
    inputBinding:
      position: 9
      prefix: '--normalizeUsing'
      valueFrom: |
        ${ 
          if( inputs.spike_in_count ){
            return null
          }
          else{
            return (self) 
          }
        }
    doc: 'method to normalize number of reads per bin RPKM, CPM, BPM, RPGC'
  - id: bin_size
    type: int
    default: 10
    inputBinding:
      position: 4
      prefix: '--binSize'
outputs:
  - id: bigwig
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + "_" + inputs.normalization + ".bigwig")
arguments:
  - position: 1
    prefix: '--extendReads'
    valueFrom: |
      ${
        if ( inputs.is_paired_end ){
           return null;
        }
        else {
          return inputs.fragment_size;
        }
      }
  - position: 10
    prefix: '--outFileName'
    valueFrom: $(inputs.bam.nameroot + "_" + inputs.normalization + ".bigwig")
  - position: 10
    prefix: '--outFileFormat'
    valueFrom: bigwig
