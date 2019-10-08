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
                "computeMatrix",
                "scale-regions"
            ],
            "inputs": [
                {
                    "type": "string",
                    "default": "computeMatrix.gz",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--outFileName"
                    },
                    "id": "#computeMatrix.cwl/outFileName"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--regionsFileName"
                    },
                    "id": "#computeMatrix.cwl/regions_bed"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2,
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
                    "id": "#main/all_samples_in_one_plot"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/bigwigs"
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
                    "run": "#computeMatrix.cwl",
                    "in": [
                        {
                            "source": "#main/regions_bed",
                            "id": "#main/computeMatrix/regions_bed"
                        },
                        {
                            "source": "#main/bigwigs",
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
                            "source": "#main/all_samples_in_one_plot",
                            "id": "#main/plotProfile/all_samples_in_one_plot"
                        },
                        {
                            "source": "#main/computeMatrix/matrix_gzip",
                            "id": "#main/plotProfile/matrixFile"
                        }
                    ],
                    "out": [
                        "#main/plotProfile/plotProfile_pdf"
                    ],
                    "id": "#main/plotProfile"
                }
            ],
            "id": "#main"
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
                    "type": "boolean",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "--plotType",
                        "valueFrom": "overlapped_lines"
                    },
                    "id": "#plotProfile.cwl/all_samples_in_one_plot"
                },
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
        }
    ],
    "cwlVersion": "v1.0"
}