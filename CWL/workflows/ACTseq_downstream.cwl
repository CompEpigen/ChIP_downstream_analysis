class: Workflow
cwlVersion: v1.2
label: peakCaller.cwl

requirements:
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
  - class: ScatterFeatureRequirement
inputs:
  - id: bam_file
    type: File
    secondaryFiles:
       - pattern: .bai
    doc: |
      bam file needs to be indexed.
      The secondary file should be of the format: file.bam.bai
      if yet not indexed use command: samtools -b file.bam
  - id: style
    type: string
    doc: |
      This parameter indicates if data is from histone marks or TFs for HOMER findPeaks
      possible values: factor or histone
  - id: 2bit_genome
    type: File?
    doc: |
      only necessary if GC bias should be computed 
      proper version for ref genome can be downloaded here: https://hgdownload.soe.ucsc.edu/
  - id: ignore_regions
    type: string
    doc: 'regions to be ignored for GC Bias check and for GC bias correction'
  - id: is_paired_end
    type: boolean
  - id: normalization_method
    type: string[]
    doc: | 
      specify normalization methods for signal tracks or give spike in count
      possibilities: RPKM, CPM, BPM, RPGC (multiple are allowed)
  - id: genome_fasta
    type: File
    doc: | 
      for annotation required. For hg19, hg38, mm10 download from: https://www.gencodegenes.org/
  - id: genome_gtf
    type: File
    doc: | 
      for annotation required. For hg19, hg38, mm10 download from: https://www.gencodegenes.org/
  - id: ignoreForNormalization
    type: string?
    doc: |
      the genomic regions that should be ignoredin bamCOmpare bigwig creation
  - id: spike_in_count
    type: long?
    doc: |
      nneds
  - id: bin_size
    type: int?
    doc: |
      bin size for signal tracks calculated with bamCoverage. default: 10
      can speed up the process significantly
  - id: min_length
    type: int?
    doc: |
      Discard reads smaller than during HOMER makeTagDirectory
  - id: max_length
    type: int?
    doc: | 
      Discard reads larger than during HOMER makeTagDirectory
  - id: computeGCbias
    type: boolean?
    default: false
    doc: |
      Calculate the GC bias and visualize it in 2 plots (#reads over GC fraction)
      Needs VERY much computing time (up to hours)!
  - id: correct_for_GCbias
    type: boolean?
    default: false
  - id: effective_genome_size
    type: long
    doc: | 
      is used for deeptools computations
      take it from:
      https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
  - id: fragment_length
    type: int?
    doc: |
      If specified, this is used for HOMER makeTagDirectory and bamCoverage to extend the reads
      As default the fragment length is estimated in HOMER
      In bamCoverage reads are only extended if data is single end
  - id: broad
    type: boolean?
    doc: |
      better for many histone modifications. default: false
      look up which which histone marks are broad: https://www.encodeproject.org/chip-seq/histone/
  - id: extsize
    type: int?
    doc: 'if --nomodel flag is set: reads are extended by this during macs2 peak calling'
  - id: shiftsize
    type: int?
    doc 'by this value reads are shifted in macs2 if --nomodel flag is set'
  - id: pvalue
    type: float?
    doc: 'adjust the pvalue for macs2 peak calling'
    default: 0.05
  - id: qvalue
    type: float?
    doc: 'for MACS2'
    default: 0.05
  - id: nomodel
    type: boolean?
    doc: |
      will skip the model building during macs2,
      which may be recommended if there is no control as it is in ACT seq
    default: true
  - id: broad_cutoff
    type: float?
  - id: peak_width
    type: int?
    doc: 'for HOMER findPeaks'
  - id: minDist_between_peaks
    type: int?
outputs:
  - id: output_tag_folder
    outputSource:
      - make_tag_directory/output_tag_folder
    doc: |
      contains tsv files for further analyses and quality control txt files
    type: Directory
  - id: GCplot
    outputSource:
      - compute_GC_bias/GCplot
    type: File?
  - id: GCbiasCheck
    outputSource:
      - compute_GC_bias/GCbiasCheck
    type: File?
  - id: narrowpeak_file
    doc: 'is format specific for MACS2 that is like a bed file with additional information'
    outputSource:
      - macs2/narrowpeak_file
    type: 'File[]'
  - id: peaks_xls
    outputSource:
      - macs2/peaks_xls
    type: File
  - id: summits_bed_file
    doc: |
      does only contain the summits of the called peaks. 
      For annotation, IDR, differential binding analysis... use .narrowPeak files
    outputSource:
      - macs2/bed_file
    type: File?
  - id: HOMER_bed_file
    outputSource:
      - pos2bed/BED_file
    type: File
  - id: bigwig
    outputSource:
      - normalized_tracks/bigwig
    type: File[]
  - id: annotated_peaks
    outputSource:
      - annotation/annotated_peaks
    type: File[]
steps:
  - id: make_tag_directory
    in:
      - id: bam_file
        source: bam_file
      - id: fragment_size
        source: fragment_length
      - id: min_length
        source: min_length
      - id: max_length
        source: max_length
    out:
      - id: output_tag_folder
    run: ./../tools/makeTagDirectory.cwl
  - id: homer_findpeaks
    in:
      - id: tag_directory
        source: make_tag_directory/output_tag_folder
      - id: style
        source: style
      - id: size
        source: peak_width
      - id: minDist_between_peaks
        source: minDist_between_peaks
    out:
      - id: peaks
    run: ./../tools/homer_findpeaks.cwl
  - id: pos2bed
    in:
      - id: peak_txt_file
        source: homer_findpeaks/peaks
    out:
      - id: BED_file
    run: ./../tools/pos2bed.cwl
  - id: compute_GC_bias
    in:
      - id: bam
        source: bam_file
      - id: genome
        source: 2bit_genome
      - id: computeGCbias
        source: computeGCbias
      - id: effective_genome_size
        source: effective_genome_size
    out:
      - id: GCbiasCheck
      - id: GCplot
    run: ./../tools/gcBias_check.cwl
    label: GCbiasCorrection
    when: $(inputs.computeGCbias)
  - id: gcbiascorrection
    in:
      - id: bam
        source: bam_file
      - id: genome
        source: 2bit_genome
      - id: GCbias_file
        source: compute_GC_bias/GCbiasCheck
      - id: correct_for_GCbias
        source: correct_for_GCbias
      - id: region
        source: ignore_regions
      - id: effective_genome_size
        source: effective_genome_size
    out:
      - id: GCcorrectedBam
    run: ./../tools/gcbiascorrection.cwl
    label: GCbiasCorrection
    when: $(inputs.correct_for_GCbias)
  - id: samtools_index
    in:
      - id: correct_for_GCbias
        source: correct_for_GCbias
      - id: bam_sorted
        source: gcbiascorrection/GCcorrectedBam
    out:
      - id: bam_sorted_indexed
    run: ./../tools/samtools_index.cwl 
    when: $(inputs.correct_for_GCbias)
  - id: macs2
    in:
      - id: bam
        source:
          - samtools_index/bam_sorted_indexed
          - bam_file
        pickValue: first_non_null
      - id: broad
        source: broad
      - id: genome_size
        source: effective_genome_size
      - id: extsize
        source: extsize
      - id: is_paired_end
        source: is_paired_end
      - id: nomodel
        source: nomodel
      - id: broad_cutoff
        source: broad_cutoff
    out:
      - id: narrowpeak_file
      - id: peaks_xls
      - id: bed_file
    run: ./../tools/macs2.cwl
  - id: normalized_tracks
    scatter: normalization
    in:
      - id: bam
        source:
          - samtools_index/bam_sorted_indexed
          - bam_file
        pickValue: first_non_null
      - id: effective_genome_size
        source: effective_genome_size
      - id: ignoreForNormalization
        source: ignoreForNormalization
      - id: is_paired_end
        source: is_paired_end
      - id: spike_in_count
        source: spike_in_count
      - id: normalization
        source: normalization_method
      - id: bin_size
        source: bin_size
    out:
      - id: bigwig
    run: ./../tools/normalized_tracks.cwl
  - id: annotation
    scatter: peak_txt_file
    in:
      - id: peak_txt_file
        source: macs2/narrowpeak_file
      - id: genome
        source: genome_fasta
      - id: gtf_file
        source: genome_gtf
    out:
      - id: annotated_peaks
    run: ./../tools/peakAnnotation.cwl
