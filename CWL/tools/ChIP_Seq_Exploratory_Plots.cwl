class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - Rscript
inputs:
  - id: ANALYSIS
    type: string
    doc: Analysis name
  - id: MARKS
    type: string
    doc: >-
      Token, that will be looked at for several samples; typically similar to
      antibodies in ChIP or assay type like ATAC seq
  - id: SAMPLES
    type: string
    doc: |-
      Token that correspond to samples; 
      It is not recommended that marks and sample names include _ or -
  - id: GROUPS
    type: string
    doc: >-
      Assignment of samples to groups; the parameter groups has to have the same
      length as "samples"
  - id: SAMPLE_WISE
    type: boolean
    default: False
    doc: Each sample will be treated as a seperated group
  - id: BW_DIR
    type: Directory
    doc: >-
      Directory with bigWig files; each file should be named:
      <BW_PREFIX>_<SampleToken>_<MarkToken>_<BW_SUFFIX>.<BW_EXTENSION>
  - id: BW_EXTENSION
    type: string
    default: "bw"
  - id: BW_PREFIX
    type: string?
    default: null
  - id: BW_SUFFIX
    type: string?
    default: null
  - id: ANNOTATION_DIR
    type: Directory
    doc: >-
      Directory with BED files: tab-separated file with columns chr  start  end 
      strand
  - id: BINNED
    type: boolean
    default: False
    doc: >-
      If True all regions from BED files are aligned at start and end (e.g. gene
      TSSs and TESs), and split into an equal number of meta-bins of different
      sizes.

      If False the regions are aligned using pivots defined in ANCHOR_SPECS
      (e.g. peak summits or gene TSSs) and split into unit windows of fixed
      size.
  - id: OFFSET
    type: int
    default: 2500
    doc: >-
      If BINNED is False: the size of the window upstream and downstream of the
      region center for which data extraction and plotting is performed
  - id: BIN_SIZE
    type: int
    default: 100
    doc: 'If BINNED is False: the size of the unit window'
  - id: N_BINS
    type: int
    default: 150
    doc: 'If BINNED is True: the number of meta windows'
  - id: ANCHOR_SPECS
    type:
      type: enum
      symbols:
        - start
        - mid
        - end
      name: ANCHOR_SPECS
    doc: 'If BINNED is False: point of alignment (start, mid, end)'
  - id: ANCHOR_NAME
    type: string
    default: "TSS"
    doc: 'If BINNED is False: name used for plotting for annotating the anchor point'
  - id: BINNED_REG_START_NAME
    type: string
    default: "TSS"
    doc: 'If BINNED is True: name of the starting point'
  - id: BINNED_REG_END_NAME
    type: string
    default: "TES"
    doc: 'If BINNED is True: name of the end point'
  - id: NORMALIZE_ATAC
    type: boolean
    doc: Normalization of ATAC data is currently not supported
  - id: NORMALIZE_CHIP
    type: boolean
    doc: Normalization of ChIP data is currently not supported
  - id: QSUB
    type: boolean
    doc: Submit bwtool jobs to a cluster; currently not supported
  - id: AVERAGES
    type: boolean
    default: True
    doc: Calculate average statistics per region
  - id: HEATMAP
    type: boolean
    default: True
    doc: Generate heatmaps of the signal
  - id: AVG_SCATTERPLOTS
    type: boolean
    default: True
    doc: Generate scatterplots with average values per region
  - id: AVG_BOXPLOTS
    type: boolean
    default: True
    doc: Generate boxplots with average values per region
  - id: COMP_GROUPS
    type: string
    doc: For average scatterplots and average boxplots sample groups to compare
  - id: COMP_MARKS
    type: string
    doc: For average scatterplots and average boxplots marks to compare
  - id: TEST_METHOD
    type:
      type: enum
      symbols:
        - averages
        - profiles
      name: TEST_METHOD
    doc: >-
      Method for testing significance; supported methods: averages: t-test with
      average values per region
  - id: PNG
    type: boolean
    doc: Generate PNGs instead of PDFs
  - id: AVG_DELTA_REF_FEATURE
    type: string?
    default: null
    doc: >-
      If defined and more than one BED file is given, the name of the BED file
      to be used as reference for comparison in average boxplots
  - id: K_MEANS_N_CLUSTERS
    type: int
    default: 1
    doc: 'for Heatmaps: number of clusters for k-Means clustering of regions'
  - id: MIN_CLUSTER_SIZE
    type: int
    default: 0
    doc: >-
      for Heatmaps: minimal number of regions in a cluster used for plotting.
      Smaller clusters will be omitted.
  - id: PLOT_MARKS
    type: string
    doc: 'for Heatmaps: mark tokens used for plotting'
  - id: CLUSTER_PROFILE_PLOTS
    type: boolean
    default: True
    doc: >-
      for Heatmaps: create a matrix of separate profile plot panels; one panel
      for each cluster
  - id: CLUSTER_PROFILE_FREE_SCALE
    type: boolean
    default: True
    doc: >-
      if CLUSTER_PROFILE_PLOTS is True: whether the scale will be dynamic or
      equal in all panels
  - id: CLUSTER_AVERAGE_PLOTS
    type: boolean
    default: True
    doc: 'for Heatmaps: create separate average scatter plots for each cluster'
  - id: CLUSTER_AVERAGE_BOXPLOTS
    type: boolean
    default: True
    doc: 'for Heatmaps: create separate average boxplots for each cluster'
  - id: HEATMAP_PANEL_SPACING
    type: float
    default: 1
    doc: 'for Heatmaps: distance between heatmap panels'
  - id: ERROR_MEASURE
    type:
      type: enum
      symbols:
        - SD
        - SEM
      name: ERROR_MEASURE
    doc: Measure of uncertainty for profile plots
  - id: QUANT
    type: float
    default: 0.01
    doc: >-
      for Heatmaps: upper (1-QUANT) and lower quantiles of the signal
      distribution that will be trimmed during plotting
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
    prefix: '--vanilla'
    valueFrom: /02_code/MainScript.R
  - position: 0
    valueFrom: ./Parameters.R
requirements:
  - class: DockerRequirement
    dockerPull: 'compepigen/enrichment_visualizations:latest'
  - class: InitialWorkDirRequirement
    listing:
      - entryname: Parameters.R
        entry: |
          SAMPLE_WISE="$(inputs.SAMPLE_WISE)"=="true"
          marks=c($(inputs.MARKS))
          SAMPLES=c($(inputs.SAMPLES))
          GROUPS=c($(inputs.GROUPS))
          ANALYSIS="$(inputs.ANALYSIS)"
          BW_DIR="$(inputs.BW_DIR.path)" 
          ANNOTATION_DIR="$(inputs.ANNOTATION_DIR.path)" 
          OFFSET=$(inputs.OFFSET)
          WIN_SIZE=$(inputs.BIN_SIZE) 
          N_BINS=$(inputs.N_BINS)
          BINNED="$(inputs.BINNED)"=="true"
          ANCHOR_SPECS="$(inputs.ANCHOR_SPECS)"
          ANCHOR_NAME="$(inputs.ANCHOR_NAME)"
          BINNED_REG_START_NAME="$(inputs.BINNED_REG_START_NAME)"
          BINNED_REG_END_NAME="$(inputs.BINNED_REG_END_NAME)"
          NORMALIZE_ATAC="$(inputs.NORMALIZE_ATAC)"=="true"
          NORMALIZE_CHIP="$(inputs.NORMALIZE_CHIP)"=="true"
          BW_EXTENSION="$(inputs.BW_EXTENSION)"
          BW_PREFIX=if("$(inputs.BW_PREFIX)" == "null") "" else "$(inputs.BW_PREFIX)"
          BW_SUFFIX=if("$(inputs.BW_SUFFIX)" == "null") "" else "$(inputs.BW_PREFIX)"
          QSUB="$(inputs.QSUB)"=="true"
          AVERAGES="$(inputs.AVERAGES)"=="true"
          HEATMAP="$(inputs.HEATMAP)"=="true"
          AVG_SCATTERPLOTS="$(inputs.AVG_SCATTERPLOTS)"=="true"
          AVG_BOXPLOTS="$(inputs.AVG_BOXPLOTS)"=="true"
          comp_groups=c($(inputs.COMP_GROUPS))
          TEST_METHOD="$(inputs.TEST_METHOD)"
          comp_marks=c($(inputs.COMP_MARKS))
          PNG="$(inputs.PNG)"=="true"
          AVG_DELTA_REF_FEATURE=if("$(inputs.AVG_DELTA_REF_FEATURE)"=="" or "$(inputs.AVG_DELTA_REF_FEATURE)"== "null") NULL else "$(inputs.AVG_DELTA_REF_FEATURE)"
          K_MEANS_N_CLUSTERS=$(inputs.K_MEANS_N_CLUSTERS)
          MIN_CLUSTER_SIZE=$(inputs.MIN_CLUSTER_SIZE)
          plot_marks=c($(inputs.PLOT_MARKS))
          CLUSTER_PROFILE_PLOTS="$(inputs.CLUSTER_PROFILE_PLOTS)"=="true"
          CLUSTER_AVERAGE_PLOTS="$(inputs.CLUSTER_AVERAGE_PLOTS)"=="true"
          CLUSTER_AVERAGE_BOXPLOTS="$(inputs.CLUSTER_AVERAGE_BOXPLOTS)"=="true"
          CLUSTER_PROFILE_FREE_SCALE="$(inputs.CLUSTER_PROFILE_FREE_SCALE)"=="true"
          GROUP_AVG_COMPARISONS=list(c($(inputs.COMP_GROUPS)))
          MARK_AVG_COMPARISONS=list(c($(inputs.COMP_MARKS)))
          PANEL_SPACING1=$(inputs.HEATMAP_PANEL_SPACING)
          PANEL_SPACING2=$(inputs.HEATMAP_PANEL_SPACING)
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
