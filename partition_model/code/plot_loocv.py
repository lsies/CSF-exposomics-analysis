from pathlib import Path
import argparse
import matplotlib
matplotlib.use("Agg")
parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parents[1] / "Results/final_loocv_predictions.csv")
parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "Figures")
args = parser.parse_args()
INPUT_PATH = args.input
OUTPUT_ROOT = args.output
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
FONT_FAMILY = 'Arial'
BASE_SIZE = 14
AXIS_TITLE_SIZE = 21
AXIS_TICK_SIZE = 21
METRICS_TEXT_SIZE = 21
AXIS_TITLE_WEIGHT = 'normal'
AXIS_TICK_WEIGHT = 'normal'
METRICS_TEXT_WEIGHT = 'normal'
POINT_SIZE = 360
POINT_ALPHA = 0.45
POINT_COLOR = '#4C92C3'
REFERENCE_LINE_COLOR = '#666666'
REFERENCE_LINE_WIDTH = 1.5
REFERENCE_LINE_STYLE = '--'
REFERENCE_LINE_ALPHA = 1.0
SHOW_METRICS_BACKGROUND = True
METRICS_X = 0.04
METRICS_Y = 0.96
METRICS_LINE_SPACING = 1.3
METRICS_BACKGROUND_COLOR = 'white'
METRICS_BACKGROUND_ALPHA = 0.9
METRICS_BORDER_COLOR = '#999999'
METRICS_BORDER_WIDTH = 0.8
FIG_WIDTH = 7.2
FIG_HEIGHT = 7.0
LAYOUT_PAD = 0.08
LAYOUT_W_PAD = 0.08
LAYOUT_H_PAD = 0.08
SAVE_PADDING = 0.15
AXIS_PADDING_RATIO = 0.16
AXIS_MINIMUM_PADDING = 0.12
OUTPUT_DPI = 600
output_png = OUTPUT_ROOT / 'Observed_LOOCV_predicted_plot.png'
output_pdf = OUTPUT_ROOT / 'Observed_LOOCV_predicted_plot.pdf'
output_tiff = OUTPUT_ROOT / 'Observed_LOOCV_predicted_plot.tiff'
df = pd.read_excel(INPUT_PATH, sheet_name='Model_input') if INPUT_PATH.suffix.lower() == '.xlsx' else pd.read_csv(INPUT_PATH)
df.columns = df.columns.astype(str).str.strip()
pass
pass
final_model_name = 'Quadratic_LogD_pH_7.4_estimated'
if 'model_name' in df.columns:
    df['model_name'] = df['model_name'].astype(str).str.strip()
    available_models = df['model_name'].dropna().unique().tolist()
    pass
    pass
    df = df.loc[df['model_name'] == final_model_name].copy()
    if df.empty:
        raise ValueError(f'Final model not found: {final_model_name}\nAvailable models: {available_models}')
required_cols = ['Compound', 'Observed_log10_median', 'LOOCV_pred_log10']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    raise ValueError(f'Missing required columns: {missing_cols}\nProvide the following columns: {required_cols}')
plot_df = df[required_cols].copy()
plot_df['Compound'] = plot_df['Compound'].astype('string').str.strip()
plot_df['Observed_log10_median'] = pd.to_numeric(plot_df['Observed_log10_median'], errors='coerce')
plot_df['LOOCV_pred_log10'] = pd.to_numeric(plot_df['LOOCV_pred_log10'], errors='coerce')
plot_df = plot_df.dropna(subset=required_cols).reset_index(drop=True)
if plot_df.empty:
    raise ValueError('No data are available for plotting.')
if len(plot_df) != 12:
    pass
else:
    pass
observed = plot_df['Observed_log10_median'].to_numpy(dtype=float)
predicted = plot_df['LOOCV_pred_log10'].to_numpy(dtype=float)
tss = np.sum((observed - np.mean(observed)) ** 2)
if np.isclose(tss, 0):
    raise ValueError('Observed_log10_median has no variation; LOOCV Q-squared cannot be calculated.')
press = np.sum((observed - predicted) ** 2)
q2_loocv = 1 - press / tss
rmse_loocv = np.sqrt(np.mean((observed - predicted) ** 2))
mae_loocv = np.mean(np.abs(observed - predicted))
if np.isclose(np.std(observed), 0) or np.isclose(np.std(predicted), 0):
    pearson_r = np.nan
else:
    pearson_r = np.corrcoef(observed, predicted)[0, 1]
pass
pass
pass
pass
pass
plt.rcdefaults()
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': [FONT_FAMILY, 'Liberation Sans', 'DejaVu Sans'], 'font.size': BASE_SIZE, 'axes.labelsize': AXIS_TITLE_SIZE, 'xtick.labelsize': AXIS_TICK_SIZE, 'ytick.labelsize': AXIS_TICK_SIZE, 'axes.unicode_minus': False, 'mathtext.fontset': 'dejavusans', 'mathtext.default': 'regular', 'pdf.fonttype': 42, 'ps.fonttype': 42})
fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=150, layout='constrained')
fig.set_constrained_layout_pads(w_pad=LAYOUT_W_PAD, h_pad=LAYOUT_H_PAD, wspace=LAYOUT_PAD, hspace=LAYOUT_PAD)
ax.scatter(observed, predicted, s=POINT_SIZE, c=POINT_COLOR, alpha=POINT_ALPHA, edgecolors='none', linewidths=0, zorder=3)
all_values = np.concatenate([observed, predicted])
data_min = np.min(all_values)
data_max = np.max(all_values)
data_range = data_max - data_min
padding = max(data_range * AXIS_PADDING_RATIO, AXIS_MINIMUM_PADDING)
axis_min = data_min - padding
axis_max = data_max + padding
ax.set_xlim(axis_min, axis_max)
ax.set_ylim(axis_min, axis_max)
ax.plot([axis_min, axis_max], [axis_min, axis_max], linestyle=REFERENCE_LINE_STYLE, linewidth=REFERENCE_LINE_WIDTH, color=REFERENCE_LINE_COLOR, alpha=REFERENCE_LINE_ALPHA, zorder=1)
ax.set_aspect('equal', adjustable='box')
metrics_text = f'LOOCV Q² = {q2_loocv:.3f}\nRMSE = {rmse_loocv:.3f}\nMAE = {mae_loocv:.3f}\nPearson r = {pearson_r:.3f}'
if SHOW_METRICS_BACKGROUND:
    metrics_bbox = dict(boxstyle='round,pad=0.35', facecolor=METRICS_BACKGROUND_COLOR, edgecolor=METRICS_BORDER_COLOR, linewidth=METRICS_BORDER_WIDTH, alpha=METRICS_BACKGROUND_ALPHA)
else:
    metrics_bbox = None
ax.text(METRICS_X, METRICS_Y, metrics_text, transform=ax.transAxes, ha='left', va='top', fontsize=METRICS_TEXT_SIZE, fontfamily=FONT_FAMILY, fontweight=METRICS_TEXT_WEIGHT, color='black', linespacing=METRICS_LINE_SPACING, bbox=metrics_bbox, clip_on=False, zorder=5)
ax.set_xlabel('Observed $\\log_{10}(R_{\\mathrm{CSF/serum}})$', fontsize=AXIS_TITLE_SIZE, fontfamily=FONT_FAMILY, fontweight=AXIS_TITLE_WEIGHT, labelpad=12)
ax.set_ylabel('LOOCV-predicted $\\log_{10}(R_{\\mathrm{CSF/serum}})$', fontsize=AXIS_TITLE_SIZE, fontfamily=FONT_FAMILY, fontweight=AXIS_TITLE_WEIGHT, labelpad=12)
ax.tick_params(axis='both', which='major', labelsize=AXIS_TICK_SIZE, direction='out', length=5, width=1.1, colors='black')
for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
    tick_label.set_fontfamily(FONT_FAMILY)
    tick_label.set_fontsize(AXIS_TICK_SIZE)
    tick_label.set_fontweight(AXIS_TICK_WEIGHT)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_color('black')
ax.grid(False)
fig.canvas.draw()
fig.savefig(output_png, dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=SAVE_PADDING, facecolor='white', edgecolor='none')
fig.savefig(output_pdf, bbox_inches='tight', pad_inches=SAVE_PADDING, facecolor='white', edgecolor='none')
fig.savefig(output_tiff, dpi=OUTPUT_DPI, bbox_inches='tight', pad_inches=SAVE_PADDING, facecolor='white', edgecolor='none', pil_kwargs={'compression': 'tiff_lzw'})
plt.close('all')
pass
pass
pass
pass
