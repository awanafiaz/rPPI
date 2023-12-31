import numpy as np
from numba import njit
from scipy.stats import norm, binom
from scipy.special import expit
from scipy.optimize import brentq, minimize
from statsmodels.regression.linear_model import OLS, WLS
from statsmodels.stats.weightstats import _zconfint_generic, _zstat_generic
from sklearn.linear_model import LogisticRegression
# from .utils import (
#     safe_expit,
#     safe_log1pexp,
#     compute_cdf,
#     compute_cdf_diff,
#     dataframe_decorator,
#     linfty_dkw,
#     linfty_binom,
#     form_discrete_distribution,
# )

import pandas as pd
import csv

def load_dataset(file_path):
    """
    Load the dataset from a CSV file.
    """
    try:
        dataset = pd.read_csv(file_path)
        return dataset
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None

def compute_cdf(Y, grid, w=None):
    """Computes the empirical CDF of the data.

    Args:
        Y (ndarray): Data.
        grid (ndarray): Grid of values to compute the CDF at.
        w (ndarray, optional): Sample weights.

    Returns:
        tuple: Empirical CDF and its standard deviation at the specified grid points.
    """
    w = np.ones(Y.shape[0]) if w is None else w / w.sum() * Y.shape[0]
    if w is None:
        indicators = (Y[:, None] <= grid[None, :]).astype(float)
    else:
        indicators = ((Y[:, None] <= grid[None, :]) * w[:, None]).astype(float)
        
    with open("indicators_py.csv", mode='w', newline='') as file:
      writer = csv.writer(file)
      writer.writerows(indicators)
        
    return indicators.mean(axis=0), indicators.std(axis=0)


def compute_cdf_diff(Y, Yhat, grid, w=None):
    """Computes the difference between the empirical CDFs of the data and the predictions.

    Args:
        Y (ndarray): Data.
        Yhat (ndarray): Predictions.
        grid (ndarray): Grid of values to compute the CDF at.
        w (ndarray, optional): Sample weights.

    Returns:
        tuple: Difference between the empirical CDFs of the data and the predictions and its standard deviation at the specified grid points.
    """
    w = np.ones(Y.shape[0]) if w is None else w / w.sum() * Y.shape[0]
    indicators_Y = (Y[:, None] <= grid[None, :]).astype(float)
    indicators_Yhat = (Yhat[:, None] <= grid[None, :]).astype(float)
    if w is None:
        return (indicators_Y - indicators_Yhat).mean(axis=0), (
            indicators_Y - indicators_Yhat
        ).std(axis=0)
    else:
        return (w[:, None] * (indicators_Y - indicators_Yhat)).mean(axis=0), (
            w[:, None] * (indicators_Y - indicators_Yhat)
        ).std(axis=0)
      
def rectified_p_value(
    rectifier,
    rectifier_std,
    imputed_mean,
    imputed_std,
    null=0,
    alternative="two-sided",
):
    """Computes a rectified p-value.

    Args:
        rectifier (float or ndarray): Rectifier value.
        rectifier_std (float or ndarray): Rectifier standard deviation.
        imputed_mean (float or ndarray): Imputed mean.
        imputed_std (float or ndarray): Imputed standard deviation.
        null (float, optional): Value of the null hypothesis to be tested. Defaults to `0`.
        alternative (str, optional): Alternative hypothesis, either 'two-sided', 'larger' or 'smaller'.

    Returns:
        float or ndarray: P-value.
    """
    rectified_point_estimate = imputed_mean + rectifier
    
    rectified_std = np.maximum(
        np.sqrt(imputed_std**2 + rectifier_std**2), 1e-16
    )
    
    return _zstat_generic(
        rectified_point_estimate, 0, rectified_std, alternative, null
    )[1]


"""
    QUANTILE ESTIMATION

"""


def _rectified_cdf(Y, Yhat, Yhat_unlabeled, grid, w=None, w_unlabeled=None):
    """Computes the rectified CDF of the data.

    Args:
        Y (ndarray): Gold-standard labels.
        Yhat (ndarray): Predictions corresponding to the gold-standard labels.
        Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
        grid (ndarray): Grid of values to compute the CDF at.
        w (ndarray, optional): Sample weights for the labeled data set.
        w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.

    Returns:
        ndarray: Rectified CDF of the data at the specified grid points.
    """
    w = np.ones(Y.shape[0]) if w is None else w / w.sum() * Y.shape[0]
    w_unlabeled = (
        np.ones(Yhat_unlabeled.shape[0])
        if w_unlabeled is None
        else w_unlabeled / w_unlabeled.sum() * Yhat_unlabeled.shape[0]
    )
    cdf_Yhat_unlabeled, _ = compute_cdf(Yhat_unlabeled, grid, w=w_unlabeled)
    
    cdf_rectifier, _ = compute_cdf_diff(Y, Yhat, grid, w=w)
    
    return cdf_Yhat_unlabeled + cdf_rectifier


def ppi_quantile_pointestimate(
    Y, Yhat, Yhat_unlabeled, q, exact_grid=False, w=None, w_unlabeled=None
):
    """Computes the prediction-powered point estimate of the quantile.

    Args:
        Y (ndarray): Gold-standard labels.
        Yhat (ndarray): Predictions corresponding to the gold-standard labels.
        Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
        q (float): Quantile to estimate.
        exact_grid (bool, optional): Whether to compute the exact solution (True) or an approximate solution based on a linearly spaced grid of 5000 values (False).
        w (ndarray, optional): Sample weights for the labeled data set.
        w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.

    Returns:
        float: Prediction-powered point estimate of the quantile.
    """
    assert len(Y.shape) == 1
    w = np.ones(Y.shape[0]) if w is None else w / w.sum() * Y.shape[0]
    w_unlabeled = (
        np.ones(Yhat_unlabeled.shape[0])
        if w_unlabeled is None
        else w_unlabeled / w_unlabeled.sum() * Yhat_unlabeled.shape[0]
    )
    grid = np.concatenate([Y, Yhat, Yhat_unlabeled], axis=0)
    
    if exact_grid:
        grid = np.sort(grid)
    else:
        grid = np.linspace(grid.min(), grid.max(), 5000)
    rectified_cdf = _rectified_cdf(
        Y, Yhat, Yhat_unlabeled, grid, w=w, w_unlabeled=w_unlabeled
    )
    
    minimizers = np.argmin(np.abs(rectified_cdf - q))
    minimizer = (
        minimizers
        if isinstance(minimizers, (int, np.int64))
        else minimizers[0]
    )
    
    return grid[
        minimizer
    ]  # Find the intersection of the rectified CDF and the quantile


def ppi_quantile_ci(
    Y,
    Yhat,
    Yhat_unlabeled,
    q,
    alpha=0.05,
    exact_grid=False,
    w=None,
    w_unlabeled=None,
):
    """Computes the prediction-powered confidence interval for the quantile.

    Args:
        Y (ndarray): Gold-standard labels.
        Yhat (ndarray): Predictions corresponding to the gold-standard labels.
        Yhat_unlabeled (ndarray): Predictions corresponding to the unlabeled data.
        q (float): Quantile to estimate. Must be in the range (0, 1).
        alpha (float, optional): Error level; the confidence interval will target a coverage of 1 - alpha. Must be in the range (0, 1).
        exact_grid (bool, optional): Whether to use the exact grid of values or a linearly spaced grid of 5000 values.
        w (ndarray, optional): Sample weights for the labeled data set.
        w_unlabeled (ndarray, optional): Sample weights for the unlabeled data set.

    Returns:
        tuple: Lower and upper bounds of the prediction-powered confidence interval for the quantile.
    """
    n = Y.shape[0]
    N = Yhat_unlabeled.shape[0]
    w = np.ones(n) if w is None else w / w.sum() * n
    w_unlabeled = (
        np.ones(N)
        if w_unlabeled is None
        else w_unlabeled / w_unlabeled.sum() * N
    )

    grid = np.concatenate([Y, Yhat, Yhat_unlabeled], axis=0)
    if exact_grid:
        grid = np.sort(grid)
    else:
        grid = np.linspace(grid.min(), grid.max(), 5000)
    cdf_Yhat_unlabeled, cdf_Yhat_unlabeled_std = compute_cdf(
        Yhat_unlabeled, grid, w=w_unlabeled
    )
    cdf_rectifier, cdf_rectifier_std = compute_cdf_diff(Y, Yhat, grid, w=w)
    
    # Calculate rectified p-value for null that the rectified cdf is equal to q
    rec_p_value = rectified_p_value(
        cdf_rectifier,
        cdf_rectifier_std / np.sqrt(n),
        cdf_Yhat_unlabeled,
        cdf_Yhat_unlabeled_std / np.sqrt(N),
        null=q,
        alternative="two-sided",
    )
    
    # Return the min and max values of the grid where p > alpha
    return grid[rec_p_value > alpha][[0, -1]]

def main():
    # Replace 'your_dataset.csv' with the actual path to your dataset file
    file_path = "C:/Users/saler/Google Drive/Fred Hutch/Research/Post-Prediction Inference/rPPI/R/dat.csv"

    # Load the dataset
    dataset = load_dataset(file_path)
    
    if dataset is not None:
      
      labeled_data = dataset[dataset["set"] == "tst"]
    
      X              = np.concatenate((np.ones((labeled_data.shape[0], 1)), labeled_data["X1"].values.reshape(-1, 1)), axis=1)  
      Y              = labeled_data["Y"].values
      Yhat           = labeled_data["Yhat"].values
    
      unlabeled_data = dataset[dataset["set"] == "val"]
    
      X_unlabeled    = np.concatenate((np.ones((unlabeled_data.shape[0], 1)), unlabeled_data["X1"].values.reshape(-1, 1)), axis=1)   
      Yhat_unlabeled = unlabeled_data["Yhat"].values
      
      print(ppi_quantile_pointestimate(Y, Yhat, Yhat_unlabeled, q = 0.5))
      
      print(ppi_quantile_ci(Y, Yhat, Yhat_unlabeled, q = 0.5))

if __name__ == "__main__":
    main()
