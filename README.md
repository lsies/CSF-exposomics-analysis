# CSF exposomics analysis

Code and numerical inputs for metabolite-driven networks and serum-to-CSF partition modeling.

## Modules

| Directory | Scope | Entry point |
|---|---|---|
| `partition_model` | Sequential descriptor models, LOOCV, compound-deletion analysis and figures | `python run_analysis.py --figures` |
| `network/metabolite_network` | Metabolite mapping, biochemical reactions, PPI, enrichment integration and figures | `python run_analysis.py` |

Run commands within the corresponding module. Use Python 3.12 and install each module's requirements in a separate environment.

## Input scope

The standalone public-expression analysis module is not included. The metabolite-network module retains the archived gene-level expression summaries needed to reproduce its exposure-response overlays and overlap analyses. Their original public-study sources remain documented in the module. This release does not rerun the upstream expression experiments or differential-expression models.

Partition-model inputs are compound-level descriptors and median partition ratios. The same values are provided in Supplementary Data 2, sheet `Model_input`.

The metabolite network starts from exported comparison statistics, feature annotations and database identifiers. It does not reprocess LC-MS raw files or recompute upstream differential-metabolite tests. P values, adjusted P values and VIP statistics are supplied as inputs. The code verifies the comparison-level BH adjustment and rebuilds the downstream network.

Metascape and Cytoscape plugin results are archived external outputs. The Python code imports them and checks associated input lists. New input gene lists require new external analyses. Metascape adjusted P values are read from the original outputs, not recalculated from filtered tables.

Public database responses, source URLs and retrieval metadata accompany the network modules. Third-party data retain their original source terms. The code license does not relicense these data.

## Outputs

Each module uses relative paths and retains inputs, code, results and figure sources. See its README for execution and input descriptions. Clinical participant-level data are not required by these entry points.
