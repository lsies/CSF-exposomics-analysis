import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error

def regression_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = mean_absolute_error(y_true, y_pred)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    q2 = np.nan if ss_tot == 0 else 1 - ss_res / ss_tot
    if np.std(y_true) == 0 or np.std(y_pred) == 0:
        r = np.nan
    else:
        r = np.corrcoef(y_true, y_pred)[0, 1]
    return {'Q2_or_R2': q2, 'RMSE': rmse, 'MAE': mae, 'Pearson_r': r}

def build_design_matrix(df, predictors, quadratic_vars=None, scaler=None):
    quadratic_vars = quadratic_vars or []
    X_raw = df[predictors].copy()
    if scaler is None:
        scaler = StandardScaler()
        X_scaled_array = scaler.fit_transform(X_raw)
    else:
        X_scaled_array = scaler.transform(X_raw)
    X_scaled = pd.DataFrame(X_scaled_array, columns=[f'z_{c}' for c in predictors], index=df.index)
    X_design = X_scaled.copy()
    for var in quadratic_vars:
        if var not in predictors:
            raise ValueError(f'Quadratic variable {var} is absent from predictors. A nonlinear predictor must also retain its linear term.')
        z_col = f'z_{var}'
        quad_col = f'z_{var}^2'
        X_design[quad_col] = X_design[z_col] ** 2
    return (X_design, scaler)

def fit_full_statsmodels(data, predictors, quadratic_vars):
    X_design, scaler = build_design_matrix(data, predictors=predictors, quadratic_vars=quadratic_vars, scaler=None)
    X_const = sm.add_constant(X_design, has_constant='add')
    y = data['log10_median']
    if X_const.shape[1] >= len(data):
        warnings.warn(f'The model has many parameters relative to its observations: n={len(data)}, parameters={X_const.shape[1]}. P values and confidence intervals may be unstable.')
    model = sm.OLS(y, X_const).fit()
    pred_log10 = model.predict(X_const)
    pred_raw = 10 ** pred_log10
    metrics = regression_metrics(y, pred_log10)
    coef_df = pd.DataFrame({'term': model.params.index, 'coef': model.params.values, 'SE': model.bse.values, 't_value': model.tvalues.values, 'p_value': model.pvalues.values, 'CI_2.5%': model.conf_int()[0].values, 'CI_97.5%': model.conf_int()[1].values})
    pred_df = data[['Compound', 'median', 'log10_median']].copy()
    pred_df['Apparent_pred_log10'] = pred_log10
    pred_df['Apparent_pred_median'] = pred_raw
    pred_df['Apparent_residual_log10'] = pred_df['log10_median'] - pred_df['Apparent_pred_log10']
    return (model, scaler, coef_df, pred_df, metrics)

def fit_sklearn_and_predict(train_data, test_data, predictors, quadratic_vars):
    X_train, scaler = build_design_matrix(train_data, predictors=predictors, quadratic_vars=quadratic_vars, scaler=None)
    X_test, _ = build_design_matrix(test_data, predictors=predictors, quadratic_vars=quadratic_vars, scaler=scaler)
    y_train = train_data['log10_median']
    model = LinearRegression()
    model.fit(X_train, y_train)
    pred = model.predict(X_test)
    return pred

def run_loocv(data, predictors, quadratic_vars):
    records = []
    for i in range(len(data)):
        test_data = data.iloc[[i]].copy()
        train_data = data.drop(data.index[i]).copy()
        pred_log10 = fit_sklearn_and_predict(train_data=train_data, test_data=test_data, predictors=predictors, quadratic_vars=quadratic_vars)[0]
        observed_log10 = test_data['log10_median'].iloc[0]
        observed_raw = test_data['median'].iloc[0]
        pred_raw = 10 ** pred_log10
        records.append({'Compound': test_data['Compound'].iloc[0], 'Observed_median': observed_raw, 'Observed_log10_median': observed_log10, 'LOOCV_pred_log10': pred_log10, 'LOOCV_pred_median': pred_raw, 'LOOCV_residual_log10': observed_log10 - pred_log10, 'Observed_over_predicted': observed_raw / pred_raw})
    loocv_df = pd.DataFrame(records)
    metrics = regression_metrics(loocv_df['Observed_log10_median'], loocv_df['LOOCV_pred_log10'])
    return (loocv_df, metrics)

STANDARDIZATION_DDOF = 0
FIRST_ORDER_DESCRIPTORS = ["HBD", "LabuteASA", "HSA_affinity", "LogD_pH_7.4_estimated"]

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
