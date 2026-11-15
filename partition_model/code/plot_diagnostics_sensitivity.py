from pathlib import Path
import argparse
import matplotlib
matplotlib.use("Agg")
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parents[1] / "Data/compound_descriptors.csv")
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "Figures")
args = parser.parse_args()
INPUT_PATH = args.input
OUTPUT_ROOT = args.output
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
FONT_FAMILY = 'Arial'
BASE_SIZE = 21
TITLE_SIZE = 24
AXIS_TITLE_SIZE = 24
AXIS_TICK_SIZE = 24
LEGEND_SIZE = 21
ANNOTATION_SIZE = 21
POINT_SIZE = 108
POINT_ALPHA = 0.45
BLUE = '#2878B5'
ORANGE = '#E98B2A'
GREEN = '#3A9D6F'
RED = '#C83E4D'
GRAY = '#666666'
DIAGNOSTIC_WIDTH = 15.0
DIAGNOSTIC_HEIGHT = 11.0
SENSITIVITY_WIDTH = 15.0
SENSITIVITY_HEIGHT = 13.0
OUTPUT_DPI = 600
OUTPUT_DIR = OUTPUT_ROOT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
STANDARDIZATION_DDOF = 0
LABEL_INFLUENTIAL_COMPOUNDS = True
plt.rcdefaults()
plt.rcParams.update({'font.family': FONT_FAMILY, 'font.size': BASE_SIZE, 'axes.titlesize': TITLE_SIZE, 'axes.titleweight': 'normal', 'axes.labelsize': AXIS_TITLE_SIZE, 'axes.labelweight': 'normal', 'xtick.labelsize': AXIS_TICK_SIZE, 'ytick.labelsize': AXIS_TICK_SIZE, 'legend.fontsize': LEGEND_SIZE, 'axes.unicode_minus': True, 'pdf.fonttype': 42, 'ps.fonttype': 42})
sns.set_style('ticks')
df = pd.read_excel(INPUT_PATH, sheet_name='Model_input') if INPUT_PATH.suffix.lower() == '.xlsx' else pd.read_csv(INPUT_PATH)
df.columns = df.columns.astype(str).str.strip().str.replace('\n', '', regex=False)
pass
pass
column_aliases = {'compound': 'Compound', 'Compound_name': 'Compound', 'Median': 'median', 'RCSF_serum_median': 'median', 'R_CSF_serum_median': 'median', 'LogD7.4': 'LogD_pH_7.4_estimated', 'LogD_7.4': 'LogD_pH_7.4_estimated', 'LogD_pH_7_4_estimated': 'LogD_pH_7.4_estimated', 'HBD_count': 'HBD', 'Labute_ASA': 'LabuteASA', 'HSA affinity': 'HSA_affinity', 'HSA_Affinity': 'HSA_affinity'}
df = df.rename(columns=column_aliases)
required_columns = ['Compound', 'median', 'LogD_pH_7.4_estimated', 'HBD', 'LabuteASA', 'HSA_affinity']
missing_columns = [col for col in required_columns if col not in df.columns]
if missing_columns:
    raise ValueError(f'Missing required columns: {missing_columns}\n\nProvide the final model input columns:\nCompound, median, LogD_pH_7.4_estimated, HBD, LabuteASA, HSA_affinity\n\nAvailable columns: {df.columns.tolist()}')
df = df[required_columns].copy()
df['Compound'] = df['Compound'].astype(str).str.strip()
numeric_columns = ['median', 'LogD_pH_7.4_estimated', 'HBD', 'LabuteASA', 'HSA_affinity']
for col in numeric_columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')
if df[required_columns].isna().any().any():
    bad_rows = df.loc[df[required_columns].isna().any(axis=1), required_columns]
    raise ValueError(f'The following rows contain missing or nonnumeric values:\n{bad_rows.to_string(index=False)}')
if len(df) != 12:
    warnings.warn(f'Read {len(df)} modeling units; the final model expects 12.')
if (df['median'] <= 0).any():
    bad_compounds = df.loc[df['median'] <= 0, 'Compound'].tolist()
    raise ValueError(f'The following compounds have nonpositive median values and cannot be log10-transformed:\n{bad_compounds}')
df['HSA_affinity'] = df['HSA_affinity'].abs()
df = df.reset_index(drop=True)
pass
FIRST_ORDER_DESCRIPTORS = ['HBD', 'LabuteASA', 'HSA_affinity', 'LogD_pH_7.4_estimated']

def calculate_scaler(data):
    means = data[FIRST_ORDER_DESCRIPTORS].mean()
    scales = data[FIRST_ORDER_DESCRIPTORS].std(ddof=STANDARDIZATION_DDOF)
    zero_sd = scales.index[np.isclose(scales, 0)].tolist()
    if zero_sd:
        raise ValueError(f'The following descriptors have zero standard deviation and cannot be standardized: {zero_sd}')
    return {'means': means, 'scales': scales}

def create_design_matrix(data, scaler):
    z_hbd = (data['HBD'] - scaler['means']['HBD']) / scaler['scales']['HBD']
    z_asa = (data['LabuteASA'] - scaler['means']['LabuteASA']) / scaler['scales']['LabuteASA']
    z_hsa = (data['HSA_affinity'] - scaler['means']['HSA_affinity']) / scaler['scales']['HSA_affinity']
    z_logd = (data['LogD_pH_7.4_estimated'] - scaler['means']['LogD_pH_7.4_estimated']) / scaler['scales']['LogD_pH_7.4_estimated']
    X = pd.DataFrame({'z_HBD': z_hbd, 'z_LabuteASA': z_asa, 'z_HSA_affinity': z_hsa, 'z_LogD7.4': z_logd, 'z_LogD7.4_squared': z_logd ** 2}, index=data.index)
    X = sm.add_constant(X, has_constant='add')
    return X

def fit_final_model(data):
    scaler = calculate_scaler(data)
    X = create_design_matrix(data, scaler)
    y = np.log10(data['median'].to_numpy(dtype=float))
    model = sm.OLS(y, X).fit()
    return (model, scaler, X, y)

def calculate_loocv(data):
    data = data.reset_index(drop=True)
    observed = np.log10(data['median'].to_numpy(dtype=float))
    predicted = np.full(len(data), np.nan, dtype=float)
    for test_position in range(len(data)):
        train_mask = np.ones(len(data), dtype=bool)
        train_mask[test_position] = False
        train_data = data.loc[train_mask].copy().reset_index(drop=True)
        test_data = data.loc[[test_position]].copy().reset_index(drop=True)
        scaler = calculate_scaler(train_data)
        X_train = create_design_matrix(train_data, scaler)
        y_train = np.log10(train_data['median'].to_numpy(dtype=float))
        model = sm.OLS(y_train, X_train).fit()
        X_test = create_design_matrix(test_data, scaler)
        predicted[test_position] = model.predict(X_test).iloc[0]
    residual = observed - predicted
    tss = np.sum((observed - observed.mean()) ** 2)
    press = np.sum(residual ** 2)
    q2 = 1 - press / tss
    rmse = np.sqrt(np.mean(residual ** 2))
    mae = np.mean(np.abs(residual))
    if np.isclose(np.std(observed), 0) or np.isclose(np.std(predicted), 0):
        pearson_r = np.nan
    else:
        pearson_r = np.corrcoef(observed, predicted)[0, 1]
    return {'observed': observed, 'predicted': predicted, 'residual': residual, 'Q2': q2, 'RMSE': rmse, 'MAE': mae, 'Pearson_r': pearson_r}

def save_figure(fig, basename):
    fig.savefig(OUTPUT_DIR / f'{basename}.png', dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=0.15)
    fig.savefig(OUTPUT_DIR / f'{basename}.tiff', dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=0.15, pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(OUTPUT_DIR / f'{basename}.pdf', bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
full_model, full_scaler, X_full, y_full = fit_final_model(df)
full_loocv = calculate_loocv(df)
full_quadratic_beta = full_model.params['z_LogD7.4_squared']
pass
pass
pass
pass
pass
pass
pass
pass
expected_values = {'R2': 0.930108, 'Adjusted_R2': 0.871865, 'Q2': 0.778594, 'RMSE': 0.168183, 'MAE': 0.144629, 'Quadratic_beta': 0.175759}
actual_values = {'R2': full_model.rsquared, 'Adjusted_R2': full_model.rsquared_adj, 'Q2': full_loocv['Q2'], 'RMSE': full_loocv['RMSE'], 'MAE': full_loocv['MAE'], 'Quadratic_beta': full_quadratic_beta}
mismatched_results = []
for key in expected_values:
    if not np.isclose(actual_values[key], expected_values[key], atol=0.003):
        mismatched_results.append((key, actual_values[key], expected_values[key]))
if mismatched_results:
    pass
    for key, actual, expected in mismatched_results:
        pass
    pass
else:
    pass
influence = full_model.get_influence()
fitted_values = full_model.fittedvalues.to_numpy()
raw_residuals = full_model.resid.to_numpy()
studentized_residuals = influence.resid_studentized_internal
cooks_distance = influence.cooks_distance[0]
leverage = influence.hat_matrix_diag
n = len(df)
p = X_full.shape[1]
cook_threshold = 4 / n
average_leverage = p / n
high_leverage_threshold = 2 * p / n
positions = np.arange(n)
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(DIAGNOSTIC_WIDTH, DIAGNOSTIC_HEIGHT), layout='constrained')
ax = axes[0, 0]
ax.scatter(fitted_values, studentized_residuals, s=POINT_SIZE, color=BLUE, alpha=POINT_ALPHA, edgecolors='none')
ax.axhline(0, color=GRAY, linestyle='--', linewidth=1.4)
ax.axhline(2, color=RED, linestyle=':', linewidth=1.2)
ax.axhline(-2, color=RED, linestyle=':', linewidth=1.2)
if LABEL_INFLUENTIAL_COMPOUNDS:
    for i, compound in enumerate(df['Compound']):
        if abs(studentized_residuals[i]) > 2:
            ax.annotate(compound, (fitted_values[i], studentized_residuals[i]), xytext=(6, 5), textcoords='offset points', fontsize=ANNOTATION_SIZE)
ax.set_xlabel('Fitted $\\log_{10}(R_{\\mathrm{CSF/serum}})$')
ax.set_ylabel('Internally studentized residual')
ax.set_title('A  Residuals versus fitted values', loc='left')
ax = axes[0, 1]
theoretical_quantiles, ordered_values = stats.probplot(studentized_residuals, dist='norm', fit=False)
ax.scatter(theoretical_quantiles, ordered_values, s=POINT_SIZE, color=BLUE, alpha=POINT_ALPHA, edgecolors='none')
q_theoretical = np.quantile(theoretical_quantiles, [0.25, 0.75])
q_observed = np.quantile(ordered_values, [0.25, 0.75])
qq_slope = (q_observed[1] - q_observed[0]) / (q_theoretical[1] - q_theoretical[0])
qq_intercept = q_observed[0] - qq_slope * q_theoretical[0]
line_x = np.array([np.min(theoretical_quantiles), np.max(theoretical_quantiles)])
line_y = qq_intercept + qq_slope * line_x
ax.plot(line_x, line_y, color=RED, linewidth=1.7)
ax.set_xlabel('Theoretical quantiles')
ax.set_ylabel('Ordered studentized residuals')
ax.set_title('B  Normal Q–Q plot', loc='left')
ax = axes[1, 0]
ax.vlines(positions, 0, cooks_distance, color=BLUE, linewidth=1.5, alpha=0.75)
ax.scatter(positions, cooks_distance, s=POINT_SIZE * 0.75, color=BLUE, alpha=POINT_ALPHA, edgecolors='none', zorder=3)
ax.axhline(cook_threshold, color=RED, linestyle='--', linewidth=1.5, label=f'$4/n={cook_threshold:.3f}$')
if LABEL_INFLUENTIAL_COMPOUNDS:
    for i, compound in enumerate(df['Compound']):
        if cooks_distance[i] > cook_threshold:
            ax.annotate(compound, (positions[i], cooks_distance[i]), xytext=(5, 5), textcoords='offset points', fontsize=ANNOTATION_SIZE)
ax.set_xticks(positions)
ax.set_xticklabels(df['Compound'], rotation=48, ha='right')
ax.set_xlabel('Compound')
ax.set_ylabel("Cook's distance")
ax.set_title("C  Cook's distance", loc='left')
ax.legend(frameon=False, loc='best')
ax = axes[1, 1]
ax.vlines(positions, 0, leverage, color=BLUE, linewidth=1.5, alpha=0.75)
ax.scatter(positions, leverage, s=POINT_SIZE * 0.75, color=BLUE, alpha=POINT_ALPHA, edgecolors='none', zorder=3)
ax.axhline(average_leverage, color=GRAY, linestyle=':', linewidth=1.3, label=f'$p/n={average_leverage:.2f}$')
ax.axhline(high_leverage_threshold, color=RED, linestyle='--', linewidth=1.5, label=f'$2p/n={high_leverage_threshold:.2f}$')
if LABEL_INFLUENTIAL_COMPOUNDS:
    for i, compound in enumerate(df['Compound']):
        if leverage[i] > high_leverage_threshold:
            ax.annotate(compound, (positions[i], leverage[i]), xytext=(5, 5), textcoords='offset points', fontsize=ANNOTATION_SIZE)
ax.set_xticks(positions)
ax.set_xticklabels(df['Compound'], rotation=48, ha='right')
ax.set_xlabel('Compound')
ax.set_ylabel('Leverage')
ax.set_title('D  Leverage', loc='left')
ax.legend(frameon=False, loc='best')
for ax in axes.flat:
    ax.grid(axis='y', linestyle=':', alpha=0.18)
    sns.despine(ax=ax)
save_figure(fig, 'Figure_S2_model_diagnostics')
sensitivity_results = []
for removed_position in range(len(df)):
    removed_compound = df.loc[removed_position, 'Compound']
    reduced_df = df.drop(index=removed_position).copy().reset_index(drop=True)
    reduced_model, _, _, _ = fit_final_model(reduced_df)
    reduced_loocv = calculate_loocv(reduced_df)
    sensitivity_results.append({'Compound_removed': removed_compound, 'n_remaining': len(reduced_df), 'Q2_LOOCV': reduced_loocv['Q2'], 'RMSE_LOOCV': reduced_loocv['RMSE'], 'MAE_LOOCV': reduced_loocv['MAE'], 'Pearson_r_LOOCV': reduced_loocv['Pearson_r'], 'Quadratic_beta': reduced_model.params['z_LogD7.4_squared']})
sensitivity_df = pd.DataFrame(sensitivity_results)
x = np.arange(len(sensitivity_df))
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(SENSITIVITY_WIDTH, SENSITIVITY_HEIGHT), sharex=True, layout='constrained')
ax = axes[0]
ax.plot(x, sensitivity_df['Q2_LOOCV'], marker='o', markersize=8, linewidth=1.6, color=BLUE)
ax.axhline(full_loocv['Q2'], color=RED, linestyle='--', linewidth=1.6, label=f"Full model: $Q^2_{{\\mathrm{{LOOCV}}}}={full_loocv['Q2']:.3f}$")
ax.set_ylabel('$Q^2_{\\mathrm{LOOCV}}$')
ax.set_title('A  Cross-validated predictive performance', loc='left')
ax.legend(frameon=False, loc='best')
ax = axes[1]
ax.plot(x, sensitivity_df['RMSE_LOOCV'], marker='o', markersize=8, linewidth=1.6, color=ORANGE)
ax.axhline(full_loocv['RMSE'], color=RED, linestyle='--', linewidth=1.6, label=f"Full model: RMSE$_{{\\mathrm{{LOOCV}}}}={full_loocv['RMSE']:.3f}$")
ax.set_ylabel('RMSE$_{\\mathrm{LOOCV}}$')
ax.set_title('B  Cross-validated prediction error', loc='left')
ax.legend(frameon=False, loc='best')
ax = axes[2]
ax.plot(x, sensitivity_df['Quadratic_beta'], marker='o', markersize=8, linewidth=1.6, color=GREEN)
ax.axhline(full_quadratic_beta, color=RED, linestyle='--', linewidth=1.6, label=f'Full model: $\\beta_{{[z(\\mathrm{{LogD}}_{{7.4}})]^2}}={full_quadratic_beta:.3f}$')
ax.axhline(0, color='black', linestyle=':', linewidth=1.3)
ax.set_ylabel(f'$\\beta_{{[z(\\mathrm{{LogD}}_{{7.4}})]^2}}$')
ax.set_title('C  Quadratic LogD$_{7.4}$ coefficient', loc='left')
ax.legend(frameon=False, loc='best')
ax.set_xticks(x)
ax.set_xticklabels(sensitivity_df['Compound_removed'], rotation=48, ha='right')
ax.set_xlabel('Compound removed')
for ax in axes:
    ax.grid(axis='y', linestyle=':', alpha=0.22)
    sns.despine(ax=ax)
save_figure(fig, 'Figure_S3_compound_sensitivity')
sensitivity_df.to_csv(OUTPUT_DIR / 'Leave_one_compound_out_results.csv', index=False, encoding='utf-8-sig')
diagnostic_df = pd.DataFrame({'Compound': df['Compound'], 'Observed_log10_RCSF_serum': y_full, 'Apparent_fitted_log10': fitted_values, 'Raw_residual': raw_residuals, 'Internally_studentized_residual': studentized_residuals, 'Cook_distance': cooks_distance, 'Leverage': leverage})
diagnostic_df.to_csv(OUTPUT_DIR / 'Final_model_diagnostic_results.csv', index=False, encoding='utf-8-sig')
pass
pass
