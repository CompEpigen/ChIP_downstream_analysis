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
                    "valueFrom": "${\n  if(self){\n    return(\"scale-regions\")\n  }\n  else{\n    return(\"reference-point\")\n  }\n}\n",
                    "inputBinding": {
                        "position": 1
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
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/after_region_start_length"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "id": "#main/before_region_start_length"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/bigwigs"
                },
                {
                    "type": "boolean",
                    "id": "#main/per_group"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "id": "#main/reference_point"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#main/regions_bed"
                },
                {
                    "type": "boolean",
                    "id": "#main/scale_regions_or_use_reference_point"
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
                            "source": "#main/after_region_start_length",
                            "id": "#main/computeMatrix/after_region_start_length"
                        },
                        {
                            "source": "#main/before_region_start_length",
                            "id": "#main/computeMatrix/before_region_start_length"
                        },
                        {
                            "source": "#main/reference_point",
                            "id": "#main/computeMatrix/reference_point"
                        },
                        {
                            "source": "#main/regions_bed",
                            "id": "#main/computeMatrix/regions_bed"
                        },
                        {
                            "source": "#main/scale_regions_or_use_reference_point",
                            "id": "#main/computeMatrix/scale_regions_or_use_reference_point"
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
                            "source": "#main/computeMatrix/matrix_gzip",
                            "id": "#main/plotProfile/matrixFile"
                        },
                        {
                            "source": "#main/per_group",
                            "id": "#main/plotProfile/per_group"
                        }
                    ],
                    "out": [
                        "#main/plotProfile/plotProfile_pdf"
                    ],
                    "id": "#main/plotProfile"
                }
            ],
            "id": "#main"
        }
    ],
    "cwlVersion": "v1.0"
}