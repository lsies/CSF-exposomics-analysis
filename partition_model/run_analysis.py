from pathlib import Path
import argparse
import subprocess
import sys
import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))
from model_core import fit_full_statsmodels, run_loocv, fit_final_model, calculate_loocv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=ROOT / 'Data' / 'compound_descriptors.csv')
    parser.add_argument('--output', type=Path, default=ROOT / 'Results')
    parser.add_argument('--figures', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = pd.read_excel(args.input, sheet_name='Model_input') if args.input.suffix.lower() == '.xlsx' else pd.read_csv(args.input)
    required = ['Compound', 'median', 'HBD', 'LabuteASA', 'HSA_affinity', 'LogD_pH_7.4_estimated', 'LogP']
    data = data[required].copy()
    if data['Compound'].duplicated().any() or data.isna().any().any():
        raise ValueError('Input requires unique compounds and complete descriptors.')
    if not np.isfinite(data[required[1:]].to_numpy(dtype=float)).all() or (data['median'] <= 0).any():
        raise ValueError('Descriptors must be finite and partition ratios must be positive.')
    data['HSA_affinity'] = data['HSA_affinity'].abs()
    data['log10_median'] = np.log10(data['median'])
    models = [('Model1', 'LogP', []), ('Model2', 'LogD_pH_7.4_estimated', []), ('Model3', 'LogD_pH_7.4_estimated', ['LogD_pH_7.4_estimated'])]
    summaries = []
    for name, lipophilicity, quadratic in models:
        predictors = ['HBD', 'LabuteASA', 'HSA_affinity', lipophilicity]
        fit, scaler, coefficients, apparent, metrics = fit_full_statsmodels(data, predictors, quadratic)
        predictions, cv = run_loocv(data, predictors, quadratic)
        summaries.append(dict(model=name, R2=fit.rsquared, Adjusted_R2=fit.rsquared_adj, Q2_LOOCV=cv['Q2_or_R2'], RMSE_LOOCV=cv['RMSE'], MAE_LOOCV=cv['MAE'], r_LOOCV=cv['Pearson_r']))
        coefficients.to_csv(output / f'{name}_coefficients.csv', index=False)
        predictions.to_csv(output / f'{name}_loocv_predictions.csv', index=False)
        if name == 'Model3':
            predictions.to_csv(output / 'final_loocv_predictions.csv', index=False)
            apparent.to_csv(output / 'final_apparent_predictions.csv', index=False)
            pd.DataFrame({'descriptor': predictors, 'mean': scaler.mean_, 'scale': scaler.scale_}).to_csv(output / 'scaler_parameters.csv', index=False)
    pd.DataFrame(summaries).to_csv(output / 'Table_S11_model_comparison.csv', index=False)
    fit, scaler, design, y = fit_final_model(data)
    cv = calculate_loocv(data)
    if not np.allclose(cv['predicted'], predictions['LOOCV_pred_log10'], rtol=1e-10, atol=1e-10):
        raise RuntimeError('The two source implementations disagree.')
    coefficients = pd.DataFrame({'term': fit.params.index, 'coefficient': fit.params.values, 'SE': fit.bse.values, 't': fit.tvalues.values, 'p_value': fit.pvalues.values, 'CI_lower': fit.conf_int()[0].values, 'CI_upper': fit.conf_int()[1].values, 'VIF': [np.nan] + [variance_inflation_factor(design.to_numpy(), i) for i in range(1, design.shape[1])]})
    coefficients.to_csv(output / 'Table_S13_final_coefficients.csv', index=False)
    sensitivity = []
    for i, compound in enumerate(data['Compound']):
        reduced = data.drop(index=i).reset_index(drop=True)
        reduced_fit, _, _, _ = fit_final_model(reduced)
        reduced_cv = calculate_loocv(reduced)
        sensitivity.append(dict(Compound_removed=compound, Q2_LOOCV=reduced_cv['Q2'], Delta_Q2=reduced_cv['Q2']-cv['Q2'], RMSE_LOOCV=reduced_cv['RMSE'], Delta_RMSE=reduced_cv['RMSE']-cv['RMSE'], MAE_LOOCV=reduced_cv['MAE'], Pearson_r=reduced_cv['Pearson_r'], Quadratic_beta=reduced_fit.params['z_LogD7.4_squared']))
    pd.DataFrame(sensitivity).to_csv(output / 'Table_S14_compound_sensitivity.csv', index=False)
    if args.figures:
        figures = output.parent / 'Figures'
        for script, source in [('plot_loocv.py', output / 'final_loocv_predictions.csv'), ('plot_adjusted_effect.py', args.input.resolve()), ('plot_diagnostics_sensitivity.py', args.input.resolve())]:
            subprocess.run([sys.executable, str(ROOT / 'code' / script), '--input', str(source), '--output', str(figures)], check=True)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == '__main__':
    main()
