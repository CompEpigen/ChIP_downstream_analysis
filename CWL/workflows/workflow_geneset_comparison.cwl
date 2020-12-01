
cwlVersion: v1.0
class: Workflow
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  regions_bed: File[]
  bam: File[]
  binSize:
    type: int 
    default: 10
  ignoreForNormalization: string?           #for example ignore chrX
  fragmentLength:
    type: int
    default: 200

outputs:
  bigwig:
    type: File[]
    outputSource: bamCoverage/bigwig

  matrix_gzip:
    type: File
    outputSource: computeMatrix/matrix_gzip

  plotProfile_pdf:
    type: File
    outputSource: plotProfile/plotProfile_pdf

  plotHeatmap_pdf:
    type: File
    outputSource: plotHeatmap/plotHeatmap_pdf

steps:
  samtools_index:
    run: ../tools/samtools_index.cwl
    scatter: bam_sorted
    in:
      bam_sorted: bam
    out: [bam_sorted_indexed]

  bamCoverage:                              #bam to bigwig
    run: ../tools/bamCoverage.cwl
    scatter: bam
    in:
      bam: samtools_index/bam_sorted_indexed
      binSize: binSize
      fragmentLength: fragmentLength
    out: [bigwig]

  computeMatrix:
    run: ../tools/computeMatrix.cwl
    in:
      regions_bed: regions_bed
      scoreFileName: bamCoverage/bigwig
    out: [matrix_gzip]

  plotProfile:
    run: ../tools/plotProfile.cwl
    in:
      matrixFile: computeMatrix/matrix_gzip
    out: [plotProfile_pdf]


  plotHeatmap:
    run: ../tools/plotHeatmap.cwl
    in:
      matrixFile: computeMatrix/matrix_gzip
    out: [plotHeatmap_pdf]