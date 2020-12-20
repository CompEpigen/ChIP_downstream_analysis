class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - Rscript
inputs:
  - id: SAMPLE_WISE
    type: boolean
  - id: marks
    type: string
  - id: SAMPLES
    type: string
  - id: GROUPS
    type: string
  - id: ANALYSIS
    type: string
  - id: BW_DIR
    type: Directory
  - id: ANNOTATION_DIR
    type: Directory
  - id: OFFSET
    type: int
  - id: WIN_SIZE
    type: int
  - id: N_BINS
    type: int
  - id: MAX_K
    type: int
  - id: BINNED
    type: boolean
  - id: ANCHOR_SPECS
    type:
      type: enum
      symbols:
        - start
        - mid
        - end
      name: ANCHOR_SPECS
  - id: ANCHOR_NAME
    type: string
  - id: BINNED_REG_START_NAME
    type: string
  - id: BINNED_REG_END_NAME
    type: string
  - id: NORMALIZE_ATAC
    type: boolean
  - id: NORMALIZE_CHIP
    type: boolean
  - id: BW_EXTENSION
    type: string
  - id: BW_PREFIX
    type: string?
  - id: BW_SUFFIX
    type: string?
  - id: AVERAGES
    type: boolean
  - id: HEATMAP
    type: boolean
  - id: AVG_SCATTERPLOTS
    type: boolean
  - id: AVG_BOXPLOTS
    type: boolean
  - id: comp_groups
    type: string
  - id: TEST_METHOD
    type:
      type: enum
      symbols:
        - averages
        - profiles
      name: TEST_METHOD
  - id: comp_marks
    type: string
  - id: PNG
    type: boolean
  - id: AVG_DELTA_REF_FEATURE
    type: string?
  - id: K_MEANS_N_CLUSTERS
    type: int
  - id: MIN_CLUSTER_SIZE
    type: int
  - id: plot_marks
    type: string
  - id: CLUSTER_PROFILE_PLOTS
    type: boolean
  - id: CLUSTER_AVERAGE_PLOTS
    type: boolean
  - id: CLUSTER_AVERAGE_BOXPLOTS
    type: boolean
  - id: CLUSTER_PROFILE_FREE_SCALE
    type: boolean
  - id: GROUP_AVG_COMPARISONS
    type: string
  - id: MARK_AVG_COMPARISONS
    type: string
  - id: PANEL_SPACING1
    type: float
  - id: PANEL_SPACING2
    type: float
  - id: ERROR_MEASURE
    type:
      type: enum
      symbols:
        - SD
        - SEM
      name: ERROR_MEASURE
  - id: QUANT
    type: float
outputs:
  - id: example_out
    type: stdout
  - id: cluster_average_plots
    type: 'File[]'
    outputBinding:
      glob: cluster_average_*.pdf
  - id: cluster_profile_plots
    type: File
    outputBinding:
      glob: cluster_profile_plots_*.pdf
  - id: heatmaps_sorted
    type: File
    outputBinding:
      glob: heatmaps_sorted_*.pdf
  - id: profile_summaries
    type: 'File[]'
    outputBinding:
      glob: profile_summaries_*.pdf
  - id: region_boxplot
    type: 'File[]'
    outputBinding:
      glob: region_boxplot_*.pdf
  - id: region_scatterplot
    type: 'File[]'
    outputBinding:
      glob: region_scatterplot_*.pdf
  - id: final_data
    type: File
    outputBinding:
      glob: final_data.RDS
  - id: final_heatmap_data
    type: File
    outputBinding:
      glob: final_heatmap_data.RDS
  - id: final_reg_data
    type: File
    outputBinding:
      glob: final_reg_data.RDS
  - id: final_rmats
    type: File
    outputBinding:
      glob: final_rmats.RDS
  - id: final_swprofiles
    type: File
    outputBinding:
      glob: final_swprofiles.RDS
  - id: kmeans_clustering
    type: 'File[]'
    outputBinding:
      glob: kmeans_clustering_*.RDS
arguments:
  - position: 0
    prefix: ''
    valueFrom: /02_code/MainScript.R
  - position: 0
    valueFrom: ./Parameters.R
requirements:
  - class: DockerRequirement
    dockerPull: 'allybuck/chip_seq_analysis_20201220:latest'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: Parameters.R
        entry: |
          SAMPLE_WISE="$(inputs.SAMPLE_WISE)"=="true"
          marks=c($(inputs.marks))
          SAMPLES=c($(inputs.SAMPLES))
          GROUPS=c($(inputs.GROUPS))
          ANALYSIS="$(inputs.ANALYSIS)"
          BW_DIR="$(inputs.BW_DIR.path)" 
          ANNOTATION_DIR="$(inputs.ANNOTATION_DIR.path)" 
          OFFSET=$(inputs.OFFSET)
          WIN_SIZE=$(inputs.WIN_SIZE) 
          N_BINS=$(inputs.N_BINS)
          MAX_K=$(inputs.MAX_K)
          BINNED="$(inputs.BINNED)"=="true"
          ANCHOR_SPECS="$(inputs.ANCHOR_SPECS)"
          ANCHOR_NAME="$(inputs.ANCHOR_NAME)"
          BINNED_REG_START_NAME="$(inputs.BINNED_REG_START_NAME)"
          BINNED_REG_END_NAME="$(inputs.BINNED_REG_END_NAME)"
          NORMALIZE_ATAC="$(inputs.NORMALIZE_ATAC)"=="true"
          NORMALIZE_CHIP="$(inputs.NORMALIZE_CHIP)"=="true"
          BW_EXTENSION="$(inputs.BW_EXTENSION)"
          BW_PREFIX="$(inputs.BW_PREFIX)"
          BW_SUFFIX="$(inputs.BW_SUFFIX)"
          QSUB="$(inputs.QSUB)"=="true"
          AVERAGES="$(inputs.AVERAGES)"=="true"
          HEATMAP="$(inputs.HEATMAP)"=="true"
          AVG_SCATTERPLOTS="$(inputs.AVG_SCATTERPLOTS)"=="true"
          AVG_BOXPLOTS="$(inputs.AVG_BOXPLOTS)"=="true"
          comp_groups=c($(inputs.comp_groups))
          TEST_METHOD="$(inputs.TEST_METHOD)"
          comp_marks=c($(inputs.comp_marks))
          PNG="$(inputs.PNG)"=="true"
          AVG_DELTA_REF_FEATURE=if("$(inputs.AVG_DELTA_REF_FEATURE)"=="") NULL else "$(inputs.AVG_DELTA_REF_FEATURE)"
          K_MEANS_N_CLUSTERS=$(inputs.K_MEANS_N_CLUSTERS)
          MIN_CLUSTER_SIZE=$(inputs.MIN_CLUSTER_SIZE)
          plot_marks=c($(inputs.plot_marks))
          CLUSTER_PROFILE_PLOTS="$(inputs.CLUSTER_PROFILE_PLOTS)"=="true"
          CLUSTER_AVERAGE_PLOTS="$(inputs.CLUSTER_AVERAGE_PLOTS)"=="true"
          CLUSTER_AVERAGE_BOXPLOTS="$(inputs.CLUSTER_AVERAGE_BOXPLOTS)"=="true"
          CLUSTER_PROFILE_FREE_SCALE="$(inputs.CLUSTER_PROFILE_FREE_SCALE)"=="true"
          GROUP_AVG_COMPARISONS=list(c($(inputs.GROUP_AVG_COMPARISONS)))
          MARK_AVG_COMPARISONS=list(c($(inputs.MARK_AVG_COMPARISONS)))
          PANEL_SPACING1=$(inputs.PANEL_SPACING1)
          PANEL_SPACING2=$(inputs.PANEL_SPACING2)
          ERROR_MEASURE="$(inputs.ERROR_MEASURE)"
          QUANT=$(inputs.QUANT)
        writable: false
      - entryname: MainScript.R
        entry: |
          args = commandArgs(trailingOnly=TRUE);
          source(args[1])
          print(test_param)
        writable: false
  - class: InlineJavascriptRequirement
stdout: output.txt
