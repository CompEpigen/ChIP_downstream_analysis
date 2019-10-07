
cwlVersion: v1.0
class: Workflow
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  regions_bed: File[]
  bigwigs: File[]

outputs:
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