# Serum to CSF partition modeling

This module reproduces the sequential descriptor models in Table S11, the final-model coefficients and variance inflation factors in Table S13, and the compound-deletion analysis in Table S14. It also produces the observed versus LOOCV-predicted plot, adjusted LogD effect plot, model diagnostics and sensitivity plots.

## Input

`Data/compound_descriptors.csv` contains the original numerical precision used for model fitting. Each row is a chemical or an isomer-group modeling unit. Required fields are `Compound`, `median`, `LogP`, `HBD`, `LabuteASA`, `HSA_affinity` and `LogD_pH_7.4_estimated`. The corresponding Excel input is supplied in Supplementary Data 2, sheet `Model_input`.

The response is log10 of the median CSF-to-serum concentration ratio. HSA affinity is the absolute docking energy in kcal/mol. Predictors are standardized with population standard deviations. Model 3 adds the square of standardized LogD. Every cross-validation fold estimates standardization from its training data. Compound-deletion analysis repeats LOOCV within each reduced dataset.

## Installation and execution

Tested with Python 3.12. Install dependencies in an isolated environment:

```bash
python -m pip install -r requirements.txt
python run_analysis.py --figures
```

The Excel workbook can be supplied directly:

```bash
python run_analysis.py --input "Supplementary Data 2.xlsx" --figures
```

Paths are resolved relative to the script. `--output` specifies a different result directory. Figures are placed in a sibling `Figures` directory. Each plot script also supports `--input` and `--output` and can be run separately.

## Expected results

| Model | LOOCV Q2 | LOOCV RMSE |
|---|---:|---:|
| Linear LogP | 0.617940 | 0.220929 |
| Linear LogD | 0.693993 | 0.197721 |
| Linear and quadratic LogD | 0.778594 | 0.168183 |

The final model has apparent R2 0.930108 and quadratic coefficient 0.175759. Results are calculated from the supplied input, not populated from these reference values.

## Files

`run_analysis.py` runs model fitting and validation. `code/model_core.py` retains the original fitting and cross-validation functions. `code/plot_loocv.py`, `code/plot_adjusted_effect.py` and `code/plot_diagnostics_sensitivity.py` retain the corresponding original plotting workflows with file-based inputs. `Results` contains numerical outputs and `Figures` contains PDF, PNG and compressed TIFF exports.

This module implements the final models described in Text S4. It does not run external-compound validation or alternative descriptor models outside that specification.
