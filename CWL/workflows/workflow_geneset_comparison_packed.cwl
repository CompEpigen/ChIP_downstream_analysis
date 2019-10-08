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
                "bamCoverage"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.bam.nameroot).bigwig",
                    "prefix": "--outFileName",
                    "position": 100
                },
                {
                    "valueFrom": "bigwig",
                    "prefix": "--outFileFormat",
                    "position": 101
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--bam"
                    },
                    "id": "#bamCoverage.cwl/bam"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--binSize"
                    },
                    "id": "#bamCoverage.cwl/binSize"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--extendReads"
                    },
                    "id": "#bamCoverage.cwl/fragmentLength"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--ignoreForNormalization"
                    },
                    "id": "#bamCoverage.cwl/ignoreForNormalization"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*filt.bigwig"
                    },
                    "id": "#bamCoverage.cwl/bigwig"
                }
            ],
            "id": "#bamCoverage.cwl"
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
                "computeMatrix",
                "scale-regions"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--afterRegionStartLength"
                    },
                    "id": "#computeMatrix.cwl/after_region_start_length"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--beforeRegionStartLength"
                    },
                    "id": "#computeMatrix.cwl/before_region_start_length"
                },
                {
                    "type": "string",
                    "default": "computeMatrix.gz",
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--outFileName"
                    },
                    "id": "#computeMatrix.cwl/outFileName"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "referencePoint"
                    },
                    "id": "#computeMatrix.cwl/reference_point"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--regionsFileName"
                    },
                    "id": "#computeMatrix.cwl/regions_bed"
                },
                {
                    "type": "boolean",
                    "default": false,
                    "inputBinding": {
                        "position": 1,
                        "valueFrom": "${\n  if(self){\n    return(\"scale-regions\")\n  }\n  else{\n    return(\"reference-point\")\n  }\n}\n"
                    },
                    "id": "#computeMatrix.cwl/scale_regions_or_use_reference_point"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 10,
                        "prefix": "--scoreFileName"
                    },
                    "id": "#computeMatrix.cwl/scoreFileName"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*Matrix.gz"
                    },
                    "id": "#computeMatrix.cwl/matrix_gzip"
                }
            ],
            "id": "#computeMatrix.cwl"
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
                "plotHeatmap"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--matrixFile"
                    },
                    "id": "#plotHeatmap.cwl/matrixFile"
                },
                {
                    "type": "string",
                    "default": "Heatmap_plot.pdf",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--outFileName"
                    },
                    "id": "#plotHeatmap.cwl/outFileName"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*plot.pdf"
                    },
                    "id": "#plotHeatmap.cwl/plotHeatmap_pdf"
                }
            ],
            "id": "#plotHeatmap.cwl"
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
                "plotProfile"
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--matrixFile"
                    },
                    "id": "#plotProfile.cwl/matrixFile"
                },
                {
                    "type": "string",
                    "default": "Profile_plot.pdf",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--outFileName"
                    },
                    "id": "#plotProfile.cwl/outFileName"
                },
                {
                    "type": "boolean",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--perGroup"
                    },
                    "id": "#plotProfile.cwl/per_group"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*plot.pdf"
                    },
                    "id": "#plotProfile.cwl/plotProfile_pdf"
                }
            ],
            "id": "#plotProfile.cwl"
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
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/bam"
                },
                {
                    "type": "int",
                    "default": 10,
                    "id": "#main/binSize"
                },
                {
                    "type": "int",
                    "default": 200,
                    "id": "#main/fragmentLength"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/ignoreForNormalization"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/regions_bed"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/bamCoverage/bigwig",
                    "id": "#main/bigwig"
                },
                {
                    "type": "File",
                    "outputSource": "#main/computeMatrix/matrix_gzip",
                    "id": "#main/matrix_gzip"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plotHeatmap/plotHeatmap_pdf",
                    "id": "#main/plotHeatmap_pdf"
                },
                {
                    "type": "File",
                    "outputSource": "#main/plotProfile/plotProfile_pdf",
                    "id": "#main/plotProfile_pdf"
                }
            ],
            "steps": [
                {
                    "run": "#bamCoverage.cwl",
                    "scatter": "#main/bamCoverage/bam",
                    "in": [
                        {
                            "source": "#main/samtools_index/bam_sorted_indexed",
                            "id": "#main/bamCoverage/bam"
                        },
                        {
                            "source": "#main/binSize",
                            "id": "#main/bamCoverage/binSize"
                        },
                        {
                            "source": "#main/fragmentLength",
                            "id": "#main/bamCoverage/fragmentLength"
                        }
                    ],
                    "out": [
                        "#main/bamCoverage/bigwig"
                    ],
                    "id": "#main/bamCoverage"
                },
                {
                    "run": "#computeMatrix.cwl",
                    "in": [
                        {
                            "source": "#main/regions_bed",
                            "id": "#main/computeMatrix/regions_bed"
                        },
                        {
                            "source": "#main/bamCoverage/bigwig",
                            "id": "#main/computeMatrix/scoreFileName"
                        }
                    ],
                    "out": [
                        "#main/computeMatrix/matrix_gzip"
                    ],
                    "id": "#main/computeMatrix"
                },
                {
                    "run": "#plotHeatmap.cwl",
                    "in": [
                        {
                            "source": "#main/computeMatrix/matrix_gzip",
                            "id": "#main/plotHeatmap/matrixFile"
                        }
                    ],
                    "out": [
                        "#main/plotHeatmap/plotHeatmap_pdf"
                    ],
                    "id": "#main/plotHeatmap"
                },
                {
                    "run": "#plotProfile.cwl",
                    "in": [
                        {
                            "source": "#main/computeMatrix/matrix_gzip",
                            "id": "#main/plotProfile/matrixFile"
                        }
                    ],
                    "out": [
                        "#main/plotProfile/plotProfile_pdf"
                    ],
                    "id": "#main/plotProfile"
                },
                {
                    "run": "#samtools_index.cwl",
                    "scatter": "#main/samtools_index/bam_sorted",
                    "in": [
                        {
                            "source": "#main/bam",
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