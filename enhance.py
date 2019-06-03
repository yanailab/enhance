#!/usr/bin/env python3

# ENHANCE denoising algorithm for single-cell RNA-Seq data
# Algorithm version: 0.1

# Author: Florian Wagner <florian.wagner@nyu.edu>
# Copyright (c) 2019 New York University

# Notes
# =====
# - Python 3 implementation, depends on scikit-learn.
# - The command-line version also depends on click and pandas.

import time
import sys
import os
from math import ceil
from typing import Tuple, Iterable, Optional, Dict, List
import hashlib
import logging

from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from scipy.stats import poisson
import numpy as np

_LOGGER = logging.getLogger(__name__)


def enhance(
        X: np.ndarray,
        target_transcript_count: Optional[int] = 200000,
        max_neighbor_frac: Optional[float] = 0.02,
        pc_var_fold_thresh: Optional[float] = 2.0,
        max_components: Optional[int] = 50,
        k: Optional[int] = None,
        use_double_precision: Optional[bool] = False,
        seed: Optional[int] = 0) -> np.ndarray:
    """Remove technical noise from a scRNA-Seq expression matrix."""
    
    _LOGGER.info('Applying ENHANCE...')
    
    if use_double_precision:
        X = np.array(X, dtype=np.float64, order='C', copy=False)
    else:
        X = np.array(X, dtype=np.float32, order='C', copy=False)
    
    _LOGGER.info('Input matrix hash: %s', get_hash(X))

    t0_total = time.time()

    transcript_count = np.median(X.sum(axis=1))
    _LOGGER.info('The median transcript count of the matrix is %.1f.',
                 transcript_count)
    
    if k is None:
        k = int(ceil(target_transcript_count / transcript_count))

        k_max = int(max_neighbor_frac * X.shape[0])
        if k <= k_max:
            _LOGGER.info(
                'Will perform denoising with k=%d '
                '(value was determined automatically '
                'based on a target transcript count of %d).',
                k, target_transcript_count)
        else:
            _LOGGER.warning(
                'Performing denoising with k=%d, to not exceed %.1f %% of '
                'the total number of cells. However, based on a target '
                'transcript count of %d, we should use k=%d. As a result, '
                'denoising results may be biased towards highly expressed '
                'genes.',
                k_max, 100*max_neighbor_frac,
                target_transcript_count, k)
            k = k_max
    else:
        _LOGGER.info('Will perform denoising with k=%d'
                     '(value was pre-specified).', k)

    # determine number of significant PCs
    _LOGGER.info('Determining the number of significant PCs...')
    sys.stdout.flush()
    
    num_components = determine_num_components(
        X, pc_var_fold_thresh, max_components,
        seed=seed)
    _LOGGER.info('The number of significant PCs is %d.', num_components)
    
    # aggregate cells
    _LOGGER.info('Aggregating cells...'); sys.stdout.flush()
    X_agg, cell_sizes = knn_aggregate(X, k, num_components, seed=seed)
    
    # denoise using PCA
    _LOGGER.info('Removing noise using PCA...'); sys.stdout.flush()
    D, scores, components, mean = denoise_pca(
        X_agg, num_components, cell_sizes,
        seed=seed)
    
    t1_total = time.time()
    _LOGGER.info('ENHANCE took %.1f s.', t1_total-t0_total)

    denoised_hash = get_hash(D)
    _LOGGER.info('Denoised matrix hash: %s', denoised_hash)

    return D, scores, components, mean, cell_sizes


def get_hash(X: np.ndarray) -> str:
    return str(hashlib.md5(X.data.tobytes()).hexdigest())


def determine_num_components(
        X: np.ndarray,
        var_fold_thresh: float = 2.0,
        max_components: int = 50,
        seed: int = 0) -> int:
    """Determine the number of significant principal components."""
    
    transcript_count = np.median(X.sum(axis=1))

    # apply PCA to real matrix
    _, real_pca_model = apply_pca(X, max_components, transcript_count, seed)
    
    # simulate pure noise matrix
    np.random.seed(seed)
    mean = normalize(X, transcript_count).mean(axis=0)
    X_noise = np.empty(X.shape, dtype=X.dtype)
    for i in range(X.shape[0]):
        X_noise[i, :] = poisson.rvs(mean)

    # apply PCA to pure noise matrix
    _, random_pca_model = apply_pca(X_noise, max_components, transcript_count, seed)
    var_thresh = var_fold_thresh * random_pca_model.explained_variance_[0]

    # determine number of components
    num_components = np.sum(real_pca_model.explained_variance_ >= var_thresh)
    
    return num_components


def normalize(
        X: np.ndarray,
        transcript_count: float = None) -> np.ndarray:
    """Perform median-normalization."""
    
    num_transcripts = X.sum(axis=1)

    if transcript_count is None:
        transcript_count = np.median(num_transcripts)

    N = ((transcript_count / num_transcripts) * X.T).T
    return N


def ft_transform(X: np.ndarray) -> np.ndarray:
    """Apply the Freeman-Tukey transformation."""
    
    # work around a bug where np.sqrt() says input is invalid for arrays
    # of type np.float32 that contain zeros
    invalid_errstate = 'warn'
    if np.issubdtype(X.dtype, np.float32):
        if np.amin(X) >= 0:
            invalid_errstate = 'ignore'
    with np.errstate(invalid=invalid_errstate):
        T = np.sqrt(X) + np.sqrt(X + 1)
    
    return T


def apply_pca(
        X, num_components: int = 50,
        transcript_count = None,
        seed: int = 0) -> Tuple[np.ndarray, PCA]:
    """Apply principal component analysis."""
    
    pca_model = PCA(
        n_components=num_components,
        svd_solver='randomized',
        random_state=seed)
    
    X_trans = ft_transform(normalize(X, transcript_count))
    scores = pca_model.fit_transform(X_trans)

    return scores, pca_model


def knn_aggregate(
        X: np.ndarray, k: int, num_components: int,
        seed: int = 0) -> np.ndarray:
    """Aggregate measurements from nearest neighbors."""

    transcript_count = np.median(X.sum(axis=1))
    
    scores, _ = apply_pca(X, num_components, transcript_count, seed=seed)
    X_agg, _ = aggregate_neighbors(X, scores, k)
    
    _, pca_model = apply_pca(X_agg, num_components, transcript_count, seed=seed)
    input_matrix = ft_transform(normalize(X, transcript_count))
    scores = pca_model.transform(input_matrix)
    X_agg, cell_sizes = aggregate_neighbors(X, scores, k)
    
    return X_agg, cell_sizes


def aggregate_neighbors(
        X: np.ndarray, scores: np.ndarray, k: int) \
        -> Tuple[np.ndarray, np.ndarray]:
    """Sub-routine for nearest neighbor aggregation."""
    
    num_transcripts = X.sum(axis=1)
    dtype = X.dtype

    # make sure score matrix is C-contiguous
    scores = np.array(scores, dtype=dtype, order='C', copy=False)
    
    # work around a bug where np.sqrt() says input is invalid for arrays
    # of type np.float32 that contain zeros
    invalid_errstate = 'warn'
    if np.issubdtype(scores.dtype, np.float32):
        invalid_errstate = 'ignore'
    with np.errstate(invalid=invalid_errstate):
        D = pairwise_distances(scores, n_jobs=1, metric='euclidean')    

    S = np.argsort(D, axis=1, kind='mergesort')
    X_agg = np.empty(X.shape, dtype=dtype)
    cell_sizes = np.empty(X.shape[0], dtype=dtype)
    for i in range(X.shape[0]):
        ind = S[i, :k]
        X_agg[i, :] = np.sum(X[ind, :], axis=0, dtype=dtype)
        cell_sizes[i] = np.median(num_transcripts[ind])
    
    return X_agg, cell_sizes


def restore_matrix(
        scores: np.ndarray, components: np.ndarray,
        mean: np.ndarray, cell_sizes: np.ndarray) -> np.ndarray:
    """Restore the expression matrix from PCA results and cell sizes."""

    # transform from PC space to original space
    D = scores.dot(components)

    # add gene means
    D = D + mean

    # invert the Freeman-Tukey transform
    D[D < 1] = 1
    D = np.power(D, 2)
    D = np.power(D-1, 2) / (4*D)
    
    D = ((cell_sizes / D.sum(axis=1)) * D.T).T
    return D


def denoise_pca(
        X: np.ndarray, num_components: int,
        cell_sizes: np.ndarray,
        seed: int = 0) -> np.ndarray:
    """Denoise data using PCA."""
    
    scores, pca_model = apply_pca(X, num_components, seed=seed)
    components = pca_model.components_
    mean = pca_model.mean_
    
    D = restore_matrix(scores, components, mean, cell_sizes)
    
    return D, scores, components, mean


def write_factorized(
        file_path: str,
        scores: np.ndarray, components: np.ndarray,
        mean: np.ndarray, cell_sizes: np.ndarray,
        cells: Iterable[str], genes: Iterable[str],
        compressed: bool = True) -> None:
    """Write ENHANCE results in factorized form."""

    file_path = os.path.expanduser(file_path)

    data = {}
    data['scores'] = np.array(scores, copy=False)
    data['components'] =  np.array(components, copy=False)
    data['mean'] = np.array(mean, copy=False)
    data['cell_sizes'] = np.array(cell_sizes, copy=False)
    data['cells'] = np.array(list(cells))
    data['genes'] = np.array(list(genes))

    if compressed:
        np.savez_compressed(file_path, **data)
    else:
        np.savez(file_path, **data)

    
def read_factorized(file_path: str) \
        -> Tuple[np.ndarray, List[str], List[str], Dict[str, np.ndarray]]:
    """Read ENHANCE output in factorized form."""
    file_path = os.path.expanduser(file_path)
    data = np.load(file_path)
    scores = data['scores']
    components = data['components']
    mean = data['mean']
    cell_sizes = data['cell_sizes']
    cells = data['cells'].tolist()
    genes = data['genes'].tolist()

    D = restore_matrix(scores, components, mean, cell_sizes)
    return D, cells, genes, data


if __name__ == '__main__':
    import click
    import pandas as pd

    ### set up logger
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    logger.handlers = []

    # create the formatter
    log_fmt = '[%(asctime)s] %(levelname)s: %(message)s'
    log_datefmt = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(log_fmt, log_datefmt)

    # create and attach a StreamHandler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    @click.command()
    @click.option('-f', '--fpath', help='The input UMI-count matrix.')
    @click.option('-o', '--saveto', default=None, help='The output matrix.')
    @click.option('--transcript-count', type=int, default=200000,
                  help='The target median transcript count for determining the'
                       'number of neighbors to use for aggregation.'
                       '(Ignored if "--num-neighbors" is specified.)')
    @click.option('--max-neighbor-frac', type=float, default=0.02,
                  help='The maximum number of neighbors to use for '
                       'aggregation, relative to the total number of cells in '
                       'the dataset. '
                       '(Ignored if "--num-neighbors" is specified.)')
    @click.option('--pc-var-fold-thresh', type=float, default=2.0,
                  help='The fold difference in variance required for relevant '
                       'PCs, relative to the variance of the first PC of '
                       'a simulated dataset containing only noise.')
    @click.option('--max-components', type=int, default=50,
                  help='The maximum number of principal components to use.')
    @click.option('--num-neighbors', type=int, default=None, show_default=True,
                  help='The number of neighbors to use for aggregation.')
    @click.option('--sep', default='\t', show_default=False,
                  help='Separator used in input file. The output file will '
                       'use this separator as well.  [default: \\t]')
    @click.option('--use-double-precision', is_flag=True,
                  help='Whether to use double-precision floating point format.'
                       ' (This doubles the amount of memory required.)')
    @click.option('--seed', default=0, show_default=True,
                  help='Seed for pseudo-random number generator.')
    @click.option('--test', is_flag=True,
                  help='Test if results for test data are correct.')
    def main(
            fpath, saveto,
            transcript_count, max_neighbor_frac, pc_var_fold_thresh,
            max_components, num_neighbors,
            sep, use_double_precision, seed, test):

        if use_double_precision:
            dtype = np.float64
        else:
            dtype = np.float32

        if saveto is None:
            _LOGGER.warning('No output file specified! Will not store results.')

        fpath_expanded = os.path.expanduser(fpath)
        file_size = os.path.getsize(fpath_expanded) / 1e6
        _LOGGER.info('Reading the expression matrix (%.1f MB) from "%s"...',
                     file_size, fpath)
        matrix = pd.read_csv(fpath_expanded, index_col=0, sep=sep)
        p, n = matrix.shape
        _LOGGER.info('The expression matrix contains %d genes and %d cells.',
                     p, n)

        result = enhance(
            matrix.values.T,
            target_transcript_count=transcript_count,
            max_neighbor_frac=max_neighbor_frac,
            pc_var_fold_thresh=pc_var_fold_thresh,
            max_components=max_components,
            k=num_neighbors,
            use_double_precision=use_double_precision,
            seed=seed)

        D = result[0]

        if saveto is not None:
            _LOGGER.info('Writing the denoised expression matrix to "%s"...',
                        saveto)
            sys.stdout.flush()
            matrix = pd.DataFrame(
                D.T, index=matrix.index, columns=matrix.columns)
            saveto_expanded = os.path.expanduser(saveto)
            matrix.to_csv(saveto_expanded, sep=sep)
            file_size = os.path.getsize(saveto_expanded) / 1e6
            _LOGGER.info('File size: %.1f MB.', file_size)

        if test:
            denoised_hash = get_hash(D)
            if denoised_hash == '3b77e847e27875a159c628759d1b60fe':
                _LOGGER.info('Output is identical!')
            else:
                _LOGGER.warning('Output is not identical!')

    main()
