{
    "$graph": [
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "genomicpariscentre/macs2:2.1.0.20140616",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 10000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "macs2",
                "callpeak"
            ],
            "arguments": [
                {
                    "valueFrom": "--nomodel",
                    "position": 3
                },
                {
                    "valueFrom": "all",
                    "prefix": "--keep-dup",
                    "position": 2
                },
                {
                    "valueFrom": "$(inputs.treatment_bam.nameroot)",
                    "prefix": "--name",
                    "position": 100
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--broad"
                    },
                    "id": "#MACS2.cwl/broad"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--control"
                    },
                    "id": "#MACS2.cwl/control_bam"
                },
                {
                    "doc": "can be \"ELAND\", \"BED\", \"ELANDMULTI\", \"ELANDEXPORT\", \"ELANDMULTIPET\", \"SAM\", \"BAM\", \"BOWTIE\", \"BAMPE\" or \"BEDPE\"",
                    "type": "string",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "--format"
                    },
                    "id": "#MACS2.cwl/format_tag"
                },
                {
                    "doc": "can be \"mm\", \"hs\", \"ce\", \"dm\", or the total number of genomic bp",
                    "type": "string",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--gsize"
                    },
                    "id": "#MACS2.cwl/genome_size"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 7,
                        "prefix": "--pvalue"
                    },
                    "id": "#MACS2.cwl/pvalue"
                },
                {
                    "type": "float",
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--qvalue"
                    },
                    "id": "#MACS2.cwl/qvalue"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--treatment"
                    },
                    "id": "#MACS2.cwl/treatment_bam"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*Peak"
                    },
                    "id": "#MACS2.cwl/peak_bed"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    },
                    "id": "#MACS2.cwl/peak_xls"
                }
            ],
            "id": "#MACS2.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/deeptools:3.1.1",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "multiBamSummary",
                "bins"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.bamfiles.nameroot + inputs.output)",
                    "prefix": "--outFileName",
                    "position": 100
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "secondaryFiles": ".bai",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--bamfiles"
                    },
                    "id": "#multiBamSummary.cwl/bamfiles"
                },
                {
                    "type": "string",
                    "id": "#multiBamSummary.cwl/output"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.npz"
                    },
                    "id": "#multiBamSummary.cwl/output_npz"
                }
            ],
            "id": "#multiBamSummary.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/deeptools:3.1.1",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "plotCorrelation"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.corData.nameroot + \"_plotCorrelation.pdf\")",
                    "prefix": "--plotFile",
                    "position": 100
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--corData"
                    },
                    "id": "#plotCorrelation.cwl/corData"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--corMethod"
                    },
                    "id": "#plotCorrelation.cwl/corMethod"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--whatToPlot"
                    },
                    "id": "#plotCorrelation.cwl/whatToPlot"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*plotCorrelation.pdf"
                    },
                    "id": "#plotCorrelation.cwl/output_plotCor"
                }
            ],
            "id": "#plotCorrelation.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/deeptools:3.1.1",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "plotPCA"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.corData.nameroot + \"_plotPCA.pdf\")",
                    "prefix": "--plotFile",
                    "position": 100
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--corData"
                    },
                    "id": "#plotPCA.cwl/corData"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*plotPCA.pdf"
                    },
                    "id": "#plotPCA.cwl/output_plotPCA"
                }
            ],
            "id": "#plotPCA.cwl"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "listing": [
                        "$(inputs.bam_sorted)"
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "kerstenbreuer/samtools:1.7",
                    "class": "DockerRequirement"
                },
                {
                    "coresMin": 1,
                    "ramMin": 20000,
                    "class": "ResourceRequirement"
                }
            ],
            "baseCommand": [
                "samtools",
                "index"
            ],
            "arguments": [
                {
                    "valueFrom": "-b",
                    "position": 1
                }
            ],
            "inputs": [
                {
                    "doc": "sorted bam input file",
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#samtools_index.cwl/bam_sorted"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "secondaryFiles": ".bai",
                    "outputBinding": {
                        "glob": "$(inputs.bam_sorted.basename)"
                    },
                    "id": "#samtools_index.cwl/bam_sorted_indexed"
                }
            ],
            "id": "#samtools_index.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#main/broad"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/control_bam"
                },
                {
                    "type": "string",
                    "default": "spearman",
                    "id": "#main/corMethod"
                },
                {
                    "type": "string",
                    "default": "AUTO",
                    "id": "#main/format_tag"
                },
                {
                    "type": "string",
                    "id": "#main/genome_size"
                },
                {
                    "type": "string",
                    "default": "results.npz",
                    "id": "#main/output"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#main/pvalue"
                },
                {
                    "type": "float",
                    "default": 0.05,
                    "id": "#main/qvalue"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/treatment_bam"
                },
                {
                    "type": "string",
                    "default": "heatmap",
                    "id": "#main/whatToPlot"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/multiBamSummary/output_npz",
                    "id": "#main/output_npz"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plotCorrelation/output_plotCor",
                    "id": "#main/output_plotCor"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plotPCA/output_plotPCA",
                    "id": "#main/output_plotPCA"
                }
            ],
            "steps": [
                {
                    "run": "#MACS2.cwl",
                    "scatter": [
                        "#main/macs2/treatment_bam",
                        "#main/macs2/control_bam"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": "#main/broad",
                            "id": "#main/macs2/broad"
                        },
                        {
                            "source": "#main/control_bam",
                            "id": "#main/macs2/control_bam"
                        },
                        {
                            "source": "#main/format_tag",
                            "id": "#main/macs2/format_tag"
                        },
                        {
                            "source": "#main/genome_size",
                            "id": "#main/macs2/genome_size"
                        },
                        {
                            "source": "#main/pvalue",
                            "id": "#main/macs2/pvalue"
                        },
                        {
                            "source": "#main/qvalue",
                            "id": "#main/macs2/qvalue"
                        },
                        {
                            "source": "#main/treatment_bam",
                            "id": "#main/macs2/treatment_bam"
                        }
                    ],
                    "out": [
                        "#main/macs2/peak_bed"
                    ],
                    "id": "#main/macs2"
                },
                {
                    "run": "#multiBamSummary.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/samtools_index/bam_sorted_indexed"
                            ],
                            "id": "#main/multiBamSummary/bamfiles"
                        },
                        {
                            "source": "#main/output",
                            "id": "#main/multiBamSummary/output"
                        }
                    ],
                    "out": [
                        "#main/multiBamSummary/output_npz"
                    ],
                    "id": "#main/multiBamSummary"
                },
                {
                    "run": "#plotCorrelation.cwl",
                    "in": [
                        {
                            "source": "#main/multiBamSummary/output_npz",
                            "id": "#main/plotCorrelation/corData"
                        },
                        {
                            "source": "#main/corMethod",
                            "id": "#main/plotCorrelation/corMethod"
                        },
                        {
                            "source": "#main/whatToPlot",
                            "id": "#main/plotCorrelation/whatToPlot"
                        }
                    ],
                    "out": [
                        "#main/plotCorrelation/output_plotCor"
                    ],
                    "id": "#main/plotCorrelation"
                },
                {
                    "run": "#plotPCA.cwl",
                    "in": [
                        {
                            "source": "#main/multiBamSummary/output_npz",
                            "id": "#main/plotPCA/corData"
                        }
                    ],
                    "out": [
                        "#main/plotPCA/output_plotPCA"
                    ],
                    "id": "#main/plotPCA"
                },
                {
                    "run": "#samtools_index.cwl",
                    "scatter": "#main/samtools_index/bam_sorted",
                    "in": [
                        {
                            "source": [
                                "#main/treatment_bam",
                                "#main/control_bam"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/samtools_index/bam_sorted"
                        }
                    ],
                    "out": [
                        "#main/samtools_index/bam_sorted_indexed"
                    ],
                    "id": "#main/samtools_index"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}