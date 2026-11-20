# CSF exposomics analysis

Code and numerical inputs for metabolite-driven networks and serum-to-CSF partition modeling.

## Modules

| Directory | Scope | Entry point |
|---|---|---|
| `partition_model` | Sequential descriptor models, LOOCV, compound-deletion analysis and figures | `python run_analysis.py --figures` |
| `network/metabolite_network` | Metabolite mapping, biochemical reactions, PPI, enrichment integration and figures | `python run_analysis.py` |

Run commands within the corresponding module. Use Python 3.12 and install each module's requirements in a separate environment.

## Input scope

This release does not include a standalone public gene-expression analysis module. The metabolite-network module retains archived gene-level expression summaries required to reproduce the exposure-response overlays and overlap analyses. The original public studies from which these data were obtained are documented within the module.

**Partition Model** inputs comprise compound-level molecular descriptors and median serum-to-CSF partition ratios. The same data are provided in the `Model_input` sheet of Supplementary Data 2.

[Supplementary Data 1](network/metabolite_network/data/Supplementary%20Data%201.xlsx) provides metabolite group summaries, pairwise comparison statistics, and reaction-to-gene mappings.

**Metabolite Network** inputs comprise exported between-group comparison statistics, feature annotations and database identifiers. P values, adjusted P values and variable importance in projection (VIP) statistics are supplied as inputs. The code checks the Benjamini-Hochberg (BH) adjustment within each comparison and reconstructs the downstream association network.

**Metascape and Cytoscape plugin results** are retained as archived external analysis outputs. The Python code imports these results and checks the corresponding input lists. New input gene lists require rerunning the relevant external analyses. Metascape adjusted P values are read directly from the original outputs rather than recalculated from filtered tables.

The network module also includes **public database responses, original source URLs and retrieval metadata**. Third-party data remain subject to the terms of their original sources.

## Outputs

Each module uses relative paths and retains inputs, code, results and figure sources. See its README for execution and input descriptions. Clinical participant-level data are not required by these entry points.
