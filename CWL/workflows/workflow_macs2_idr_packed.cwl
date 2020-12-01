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
            "baseCommand": [
                "idr"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "inputBinding": {
                        "position": 5,
                        "prefix": "--idr-threshold"
                    },
                    "id": "#IDR.cwl/idrThreshold"
                },
                {
                    "type": "string",
                    "default": "narrowPeak",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--input-file-type"
                    },
                    "id": "#IDR.cwl/inputFileType"
                },
                {
                    "type": "string",
                    "default": "idrValue.txt",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--output-file"
                    },
                    "id": "#IDR.cwl/outputFile"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--plot"
                    },
                    "id": "#IDR.cwl/plot"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--samples"
                    },
                    "id": "#IDR.cwl/samples"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.outputFile)"
                    },
                    "id": "#IDR.cwl/peaks_idr_scores"
                },
                {
                    "type": "stderr",
                    "id": "#IDR.cwl/stderr_out"
                }
            ],
            "id": "#IDR.cwl"
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
                    "valueFrom": "$(inputs.treatment_bam_1.nameroot + \"_replicate1\")",
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
                    "id": "#MACS2_1.cwl/broad"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--control"
                    },
                    "id": "#MACS2_1.cwl/control_bam_1"
                },
                {
                    "doc": "can be \"ELAND\", \"BED\", \"ELANDMULTI\", \"ELANDEXPORT\", \"ELANDMULTIPET\", \"SAM\", \"BAM\", \"BOWTIE\", \"BAMPE\" or \"BEDPE\"",
                    "type": "string",
                    "default": "AUTO",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "--format"
                    },
                    "id": "#MACS2_1.cwl/format_tag"
                },
                {
                    "doc": "can be \"mm\", \"hs\", \"ce\", \"dm\", or the total number of genomic bp",
                    "type": "string",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--gsize"
                    },
                    "id": "#MACS2_1.cwl/genome_size"
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
                    "id": "#MACS2_1.cwl/pvalue"
                },
                {
                    "type": "float",
                    "default": 0.05,
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--qvalue"
                    },
                    "id": "#MACS2_1.cwl/qvalue"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--treatment"
                    },
                    "id": "#MACS2_1.cwl/treatment_bam_1"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*Peak"
                    },
                    "id": "#MACS2_1.cwl/peak1_bed"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    },
                    "id": "#MACS2_1.cwl/peak1_xls"
                }
            ],
            "id": "#MACS2_1.cwl"
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
                    "valueFrom": "$(inputs.treatment_bam_2.nameroot + \"_replicate2\")",
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
                    "id": "#MACS2_2.cwl/broad"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--control"
                    },
                    "id": "#MACS2_2.cwl/control_bam_2"
                },
                {
                    "doc": "can be \"ELAND\", \"BED\", \"ELANDMULTI\", \"ELANDEXPORT\", \"ELANDMULTIPET\", \"SAM\", \"BAM\", \"BOWTIE\", \"BAMPE\" or \"BEDPE\"",
                    "type": "string",
                    "default": "AUTO",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "--format"
                    },
                    "id": "#MACS2_2.cwl/format_tag"
                },
                {
                    "doc": "can be \"mm\", \"hs\", \"ce\", \"dm\", or the total number of genomic bp",
                    "type": "string",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--gsize"
                    },
                    "id": "#MACS2_2.cwl/genome_size"
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
                    "id": "#MACS2_2.cwl/pvalue"
                },
                {
                    "type": "float",
                    "default": 0.05,
                    "inputBinding": {
                        "position": 6,
                        "prefix": "--qvalue"
                    },
                    "id": "#MACS2_2.cwl/qvalue"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--treatment"
                    },
                    "id": "#MACS2_2.cwl/treatment_bam_2"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*Peak"
                    },
                    "id": "#MACS2_2.cwl/peak2_bed"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_peaks.xls"
                    },
                    "id": "#MACS2_2.cwl/peak2_xls"
                }
            ],
            "id": "#MACS2_2.cwl"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "default": false,
                    "id": "#main/broad"
                },
                {
                    "type": "File",
                    "id": "#main/control_bam_1"
                },
                {
                    "type": "File",
                    "id": "#main/control_bam_2"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/filetype"
                },
                {
                    "type": "string",
                    "id": "#main/genome_size"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/output_basename"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "id": "#main/plot"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#main/pvalue"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "id": "#main/qvalue"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/tagformat"
                },
                {
                    "type": "File",
                    "id": "#main/treatment_bam_1"
                },
                {
                    "type": "File",
                    "id": "#main/treatment_bam_2"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/idr/peaks_idr_scores",
                    "id": "#main/idr_values"
                },
                {
                    "type": "File",
                    "outputSource": "#main/idr/stderr_out",
                    "id": "#main/std_error"
                }
            ],
            "steps": [
                {
                    "run": "#IDR.cwl",
                    "in": [
                        {
                            "source": "#main/filetype",
                            "id": "#main/idr/inputFileType"
                        },
                        {
                            "source": "#main/output_basename",
                            "id": "#main/idr/output_basename"
                        },
                        {
                            "source": [
                                "#main/macs2_1/peak1_bed",
                                "#main/macs2_2/peak2_bed"
                            ],
                            "id": "#main/idr/samples"
                        }
                    ],
                    "out": [
                        "#main/idr/peaks_idr_scores",
                        "#main/idr/stderr_out"
                    ],
                    "id": "#main/idr"
                },
                {
                    "run": "#MACS2_1.cwl",
                    "in": [
                        {
                            "source": "#main/broad",
                            "id": "#main/macs2_1/broad"
                        },
                        {
                            "source": "#main/control_bam_1",
                            "id": "#main/macs2_1/control_bam_1"
                        },
                        {
                            "source": "#main/tagformat",
                            "id": "#main/macs2_1/format_tag"
                        },
                        {
                            "source": "#main/genome_size",
                            "id": "#main/macs2_1/genome_size"
                        },
                        {
                            "source": "#main/pvalue",
                            "id": "#main/macs2_1/pvalue"
                        },
                        {
                            "source": "#main/qvalue",
                            "id": "#main/macs2_1/qvalue"
                        },
                        {
                            "source": "#main/treatment_bam_1",
                            "id": "#main/macs2_1/treatment_bam_1"
                        }
                    ],
                    "out": [
                        "#main/macs2_1/peak1_bed"
                    ],
                    "id": "#main/macs2_1"
                },
                {
                    "run": "#MACS2_2.cwl",
                    "in": [
                        {
                            "source": "#main/broad",
                            "id": "#main/macs2_2/broad"
                        },
                        {
                            "source": "#main/control_bam_2",
                            "id": "#main/macs2_2/control_bam_2"
                        },
                        {
                            "source": "#main/tagformat",
                            "id": "#main/macs2_2/format_tag"
                        },
                        {
                            "source": "#main/genome_size",
                            "id": "#main/macs2_2/genome_size"
                        },
                        {
                            "source": "#main/pvalue",
                            "id": "#main/macs2_2/pvalue"
                        },
                        {
                            "source": "#main/qvalue",
                            "id": "#main/macs2_2/qvalue"
                        },
                        {
                            "source": "#main/treatment_bam_2",
                            "id": "#main/macs2_2/treatment_bam_2"
                        }
                    ],
                    "out": [
                        "#main/macs2_2/peak2_bed"
                    ],
                    "id": "#main/macs2_2"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}