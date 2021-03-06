
cwlVersion: v1.0
class: Workflow
requirements:
  MultipleInputFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  regions_bed: File[]
  bigwigs: File[]
  per_group: boolean
  scale_regions_or_use_reference_point: boolean
  before_region_start_length: int?
  after_region_start_length: int?
  reference_point: string?

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
      scoreFileName: bigwigs
      scale_regions_or_use_reference_point: scale_regions_or_use_reference_point
      before_region_start_length: before_region_start_length
      after_region_start_length: after_region_start_length
      reference_point: reference_point
    out: [matrix_gzip]

  plotProfile:
    run: ../tools/plotProfile.cwl
    in:
      matrixFile: computeMatrix/matrix_gzip
      per_group: per_group
    out: [plotProfile_pdf]


  plotHeatmap:
    run: ../tools/plotHeatmap.cwl
    in:
      matrixFile: computeMatrix/matrix_gzip
    out: [plotHeatmap_pdf]