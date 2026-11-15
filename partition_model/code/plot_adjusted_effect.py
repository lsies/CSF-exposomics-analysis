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
import statsmodels.api as sm
FONT_FAMILY = 'Arial'
BASE_SIZE = 14
AXIS_TITLE_SIZE = 24
AXIS_TICK_SIZE = 24
ANNOTATION_SIZE = 21
LEGEND_SIZE = 21
FIG_WIDTH = 8.2
FIG_HEIGHT = 7.2
POINT_SIZE = 180
POINT_COLOR = '#2878B5'
POINT_ALPHA = 0.45
LINE_COLOR = '#C83E4D'
LINE_WIDTH = 2.6
CI_COLOR = '#C83E4D'
CI_ALPHA = 0.2
RUG_COLOR = '#333333'
RUG_HEIGHT_FRACTION = 0.025
RUG_LINEWIDTH = 1.2
SHOW_ADJUSTED_POINTS = True
SHOW_RUG = True
SHOW_MODEL_ANNOTATION = True
SHOW_LEGEND = True
ANNOTATION_X = 0.04
ANNOTATION_Y = 0.96
LEGEND_LOCATION = 'upper right'
USE_ORIGINAL_LOGD_SCALE = True
GRID_POINTS = 300
STANDARDIZATION_DDOF = 0
OUTPUT_DPI = 600
OUTPUT_PNG = OUTPUT_ROOT / 'Adjusted_LogD7.4_effect_curve.png'
OUTPUT_PDF = OUTPUT_ROOT / 'Adjusted_LogD7.4_effect_curve.pdf'
OUTPUT_TIFF = OUTPUT_ROOT / 'Adjusted_LogD7.4_effect_curve.tiff'
df = pd.read_excel(INPUT_PATH, sheet_name='Model_input') if INPUT_PATH.suffix.lower() == '.xlsx' else pd.read_csv(INPUT_PATH)
df.columns = df.columns.astype(str).str.strip().str.replace('\n', '', regex=False)
pass
pass
column_aliases = {'compound': 'Compound', 'Compound_name': 'Compound', 'Median': 'median', 'R_CSF_serum_median': 'median', 'RCSF_serum_median': 'median', 'Labute_ASA': 'LabuteASA', 'HSA affinity': 'HSA_affinity', 'HSA_Affinity': 'HSA_affinity', 'LogD7.4': 'LogD_pH_7.4_estimated', 'LogD_7.4': 'LogD_pH_7.4_estimated', 'LogD_pH_7_4_estimated': 'LogD_pH_7.4_estimated', 'HBD_count': 'HBD'}
df = df.rename(columns=column_aliases)
required_columns = ['Compound', 'median', 'LabuteASA', 'HSA_affinity', 'LogD_pH_7.4_estimated', 'HBD']
missing_columns = [col for col in required_columns if col not in df.columns]
if missing_columns:
    raise ValueError(f'Missing required columns:\n{missing_columns}\n\nProvide these six input columns:\nCompound, median, LabuteASA, HSA_affinity, LogD_pH_7.4_estimated, HBD')
plot_df = df[required_columns].copy()
plot_df['Compound'] = plot_df['Compound'].astype('string').str.strip()
numeric_columns = ['median', 'LabuteASA', 'HSA_affinity', 'LogD_pH_7.4_estimated', 'HBD']
for col in numeric_columns:
    plot_df[col] = pd.to_numeric(plot_df[col], errors='coerce')
before_rows = len(plot_df)
plot_df = plot_df.dropna(subset=required_columns).reset_index(drop=True)
removed_rows = before_rows - len(plot_df)
if removed_rows > 0:
    warnings.warn(f'Removed {removed_rows} rows containing missing or nonnumeric values.')
if plot_df.empty:
    raise ValueError('No observations remain for modeling after data cleaning.')
if (plot_df['median'] <= 0).any():
    bad_compounds = plot_df.loc[plot_df['median'] <= 0, 'Compound'].tolist()
    raise ValueError(f'The following compounds have nonpositive median values and cannot be log10-transformed:\n{bad_compounds}')
if len(plot_df) != 12:
    warnings.warn(f'Read {len(plot_df)} compounds; the final model expects 12 compounds. Check the input data.')
pass
plot_df['log10_median'] = np.log10(plot_df['median'])
descriptor_columns = ['LabuteASA', 'HSA_affinity', 'LogD_pH_7.4_estimated', 'HBD']
standardization_parameters = {}
for col in descriptor_columns:
    mean_value = plot_df[col].mean()
    sd_value = plot_df[col].std(ddof=STANDARDIZATION_DDOF)
    if np.isclose(sd_value, 0):
        raise ValueError(f'{col} has zero standard deviation across compounds and cannot be standardized.')
    standardization_parameters[col] = {'mean': mean_value, 'sd': sd_value}
    plot_df[f'z_{col}'] = (plot_df[col] - mean_value) / sd_value
plot_df['z_LogD_squared'] = plot_df['z_LogD_pH_7.4_estimated'] ** 2
model_columns = ['z_LabuteASA', 'z_HSA_affinity', 'z_LogD_pH_7.4_estimated', 'z_HBD', 'z_LogD_squared']
X = plot_df[model_columns].copy()
X = sm.add_constant(X, has_constant='add')
y = plot_df['log10_median']
model = sm.OLS(y, X).fit()
pass
pass
pass
pass
pass
for term, coefficient in model.params.items():
    pass
pass
pass
beta_logd = model.params['z_LogD_pH_7.4_estimated']
beta_logd_squared = model.params['z_LogD_squared']
p_logd = model.pvalues['z_LogD_pH_7.4_estimated']
p_logd_squared = model.pvalues['z_LogD_squared']
ci_logd_squared = model.conf_int().loc['z_LogD_squared']
pass
pass
pass
logd_min = plot_df['LogD_pH_7.4_estimated'].min()
logd_max = plot_df['LogD_pH_7.4_estimated'].max()
logd_grid = np.linspace(logd_min, logd_max, GRID_POINTS)
logd_mean = standardization_parameters['LogD_pH_7.4_estimated']['mean']
logd_sd = standardization_parameters['LogD_pH_7.4_estimated']['sd']
z_logd_grid = (logd_grid - logd_mean) / logd_sd
prediction_X = pd.DataFrame({'const': np.ones(GRID_POINTS), 'z_LabuteASA': np.zeros(GRID_POINTS), 'z_HSA_affinity': np.zeros(GRID_POINTS), 'z_LogD_pH_7.4_estimated': z_logd_grid, 'z_HBD': np.zeros(GRID_POINTS), 'z_LogD_squared': z_logd_grid ** 2})
prediction_X = prediction_X[X.columns]
prediction_results = model.get_prediction(prediction_X).summary_frame(alpha=0.05)
predicted_mean = prediction_results['mean'].to_numpy()
ci_lower = prediction_results['mean_ci_lower'].to_numpy()
ci_upper = prediction_results['mean_ci_upper'].to_numpy()
model_residuals = model.resid.to_numpy()
observed_z_logd = plot_df['z_LogD_pH_7.4_estimated'].to_numpy()
adjusted_points = model.params['const'] + beta_logd * observed_z_logd + beta_logd_squared * observed_z_logd ** 2 + model_residuals
if USE_ORIGINAL_LOGD_SCALE:
    x_curve = logd_grid
    x_points = plot_df['LogD_pH_7.4_estimated'].to_numpy()
    x_label = 'Estimated $\\mathrm{LogD}_{7.4}$'
else:
    x_curve = z_logd_grid
    x_points = observed_z_logd
    x_label = 'Standardized $\\mathrm{LogD}_{7.4}$'
plt.rcdefaults()
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [FONT_FAMILY, 'Liberation Sans', 'DejaVu Sans'], 'font.size': BASE_SIZE, 'axes.unicode_minus': False, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'mathtext.fontset': 'stixsans'})
fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT), constrained_layout=True)
ax.fill_between(x_curve, ci_lower, ci_upper, color=CI_COLOR, alpha=CI_ALPHA, linewidth=0, label='95% CI', zorder=1)
ax.plot(x_curve, predicted_mean, color=LINE_COLOR, linewidth=LINE_WIDTH, label='Adjusted effect', zorder=3)
if SHOW_ADJUSTED_POINTS:
    ax.scatter(x_points, adjusted_points, s=POINT_SIZE, color=POINT_COLOR, edgecolors='none', linewidths=0, alpha=POINT_ALPHA, label='Adjusted observations', zorder=4)
if SHOW_RUG:
    all_y_values = np.concatenate([predicted_mean, ci_lower, ci_upper, adjusted_points])
    y_min = np.nanmin(all_y_values)
    y_max = np.nanmax(all_y_values)
    y_range = y_max - y_min
    if np.isclose(y_range, 0):
        y_range = 1
    rug_bottom = y_min - 0.08 * y_range
    rug_top = rug_bottom + RUG_HEIGHT_FRACTION * y_range
    for x_value in x_points:
        ax.plot([x_value, x_value], [rug_bottom, rug_top], color=RUG_COLOR, linewidth=RUG_LINEWIDTH, alpha=0.85, clip_on=False, zorder=2)
    ax.set_ylim(rug_bottom - 0.025 * y_range, y_max + 0.08 * y_range)
legend = ax.legend(loc='upper right', bbox_to_anchor=(1.02, 0.99), frameon=False, fontsize=LEGEND_SIZE, markerscale=1.2, handlelength=1.2, handletextpad=0.5, labelspacing=0.5, borderaxespad=0.2, scatterpoints=1)
legend.set_zorder(8)
if SHOW_MODEL_ANNOTATION:
    annotation_text = f'$\\beta_{{\\mathrm{{LogD}}^2}}$ = {beta_logd_squared:.3f}\n95% CI: {ci_logd_squared.iloc[0]:.3f}, {ci_logd_squared.iloc[1]:.3f}\n$p$ = {p_logd_squared:.3f}'
    ax.text(ANNOTATION_X, ANNOTATION_Y, annotation_text, transform=ax.transAxes, ha='left', va='top', fontsize=ANNOTATION_SIZE, fontfamily=FONT_FAMILY, linespacing=1.25, bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=1.0, edgecolor='#999999', linewidth=0.8), zorder=10)
ax.set_xlabel(x_label, fontsize=AXIS_TITLE_SIZE, fontweight='normal', labelpad=10)
ax.set_ylabel('Adjusted predicted $\\log_{10}(R_{\\mathrm{CSF/serum}})$', fontsize=AXIS_TITLE_SIZE, fontweight='normal', labelpad=10)
ax.tick_params(axis='both', which='major', labelsize=AXIS_TICK_SIZE, direction='out', length=5, width=1.1, colors='black')
for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
    tick_label.set_fontfamily(FONT_FAMILY)
    tick_label.set_fontsize(AXIS_TICK_SIZE)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_color('black')
ax.grid(False)
fig.savefig(OUTPUT_PNG, dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=0.1, facecolor='white')
fig.savefig(OUTPUT_PDF, bbox_inches='tight', pad_inches=0.1, facecolor='white')
fig.savefig(OUTPUT_TIFF, dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=0.1, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
plt.close('all')
pass
pass
pass
pass
