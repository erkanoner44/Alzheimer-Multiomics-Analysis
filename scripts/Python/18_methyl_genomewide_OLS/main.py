#!/usr/bin/env python3
"""GSE134379 MTG: genome-wide age/sex/plate-adjusted exploratory beta OLS.

NOT raw-array QC; detection p columns duplicated; no cell-composition adjustment.
"""
import argparse
import csv
import gzip
import json
import math
import re
import traceback
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests

FILES = {'AD': 'GSE134379_processedSamples_mtgAD.txt.gz', 'ND': 'GSE134379_processedSamples_mtgND.txt.gz'}
HEADER_RE = re.compile(r'^Beta value for\s+MTG\s+(.+?)\s*$', re.I)
CG_RE = re.compile(r'^cg\d+$', re.I)


def write_audit(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')


def header_info(path, group, md):
    if not path.is_file():
        raise FileNotFoundError(str(path))
    with gzip.open(path, 'rt', encoding='utf-8-sig', errors='replace') as f:
        h = f.readline().rstrip('\r\n').split('\t')
    if not h or len(h) % 2:
        raise ValueError(f'{group}: beta/detection header length invalid')
    ids = []
    for j in range(0, len(h), 2):
        match = HEADER_RE.fullmatch(h[j].strip())
        if match is None or h[j + 1].strip().lower() != 'detection p value':
            raise ValueError(f'{group}: unexpected header pair {h[j:j+2]}')
        ids.append(match.group(1).strip())
    if len(ids) != len(set(ids)):
        raise ValueError(f'{group}: duplicate subject IDs')
    target = set(md.loc[md.diagnosis == group, 'subject_id'])
    if set(ids) != target:
        raise ValueError(f'{group}: GEO metadata subject mismatch: {list(set(ids)^target)[:5]}')
    return ids, len(h) + 1


def design_metadata(path):
    if not path.is_file():
        raise FileNotFoundError(f'Module17 verified metadata not found: {path}')
    df = pd.read_csv(path, dtype={'subject_id': str, 'plate': str}, encoding='utf-8-sig')
    required = {'subject_id', 'diagnosis', 'age_numeric', 'gender', 'plate', 'brain_region'}
    if not required.issubset(df.columns):
        raise ValueError(f'Missing metadata columns: {sorted(required-set(df.columns))}')
    df = df.loc[df.brain_region == 'MTG'].copy()
    if df.subject_id.isna().any() or df.subject_id.duplicated().any():
        raise ValueError('Missing/duplicate MTG subject ID')
    if Counter(df.diagnosis) != Counter({'AD': 225, 'ND': 179}):
        raise ValueError(f'Incorrect MTG group sizes: {Counter(df.diagnosis)}')
    if not df.gender.isin(['M', 'F']).all():
        raise ValueError('Gender must be M/F')
    if df[['age_numeric', 'plate']].isna().any().any():
        raise ValueError('Age / plate is missing')
    df['age_numeric'] = pd.to_numeric(df.age_numeric, errors='raise')
    if not np.isfinite(df.age_numeric).all():
        raise ValueError('Age is not finite')
    return df


def design(md, ids_ad, ids_nd):
    ordered = md.set_index('subject_id').loc[ids_ad + ids_nd]
    if not ordered.diagnosis.iloc[:225].eq('AD').all() or not ordered.diagnosis.iloc[225:].eq('ND').all():
        raise ValueError('Sample order and diagnosis disagree')
    cov = pd.DataFrame({'intercept': np.ones(len(ordered)),
                        'diagnosis_AD': ordered.diagnosis.eq('AD').astype(float).to_numpy(),
                        'age_per_10_years': ((ordered.age_numeric - ordered.age_numeric.mean())/10).to_numpy(),
                        'male': ordered.gender.eq('M').astype(float).to_numpy()})
    dummy = pd.get_dummies(ordered.plate.astype(str), prefix='plate', drop_first=True, dtype=float).reset_index(drop=True)
    cov = pd.concat([cov, dummy], axis=1)
    X = cov.to_numpy(dtype=np.float64, copy=True)
    rank = int(np.linalg.matrix_rank(X))
    condition = float(np.linalg.cond(X))
    if rank != X.shape[1] or not np.isfinite(condition) or condition > 1e8:
        raise ValueError(f'Rank/conditioning invalid: rank={rank}, cond={condition}')
    return X, {'n_columns': X.shape[1], 'rank': rank, 'condition_number': condition,
               'columns': list(cov.columns), 'plate_by_group': pd.crosstab(md.plate, md.diagnosis).to_dict()}


def parse_beta(values):
    # Explicitly return writable float64 array, independent of Pandas read-only buffers.
    ar = np.empty(len(values), dtype=np.float64)
    for i, v in enumerate(values):
        try:
            ar[i] = float(v)
        except (ValueError, TypeError):
            ar[i] = np.nan
    invalid = (~np.isfinite(ar)) | (ar < 0) | (ar > 1)
    num_invalid = int(invalid.sum())
    ar[invalid] = np.nan
    return ar, num_invalid


def tests_block(Y, X, min_AD, min_ND):
    """Vectorized OLS for fully observed rows; stable per-row QR fallback for missingness."""
    k, nsamples = Y.shape
    p = X.shape[1]
    counts_AD = np.isfinite(Y[:, :225]).sum(axis=1)
    counts_ND = np.isfinite(Y[:, 225:]).sum(axis=1)
    estimates = np.full(k, np.nan)
    pvals = np.full(k, np.nan)
    low = np.full(k, np.nan)
    high = np.full(k, np.nan)
    means_AD = np.divide(np.nansum(Y[:, :225], axis=1), counts_AD,
                         out=np.full(k, np.nan), where=counts_AD > 0)
    means_ND = np.divide(np.nansum(Y[:, 225:], axis=1), counts_ND,
                         out=np.full(k, np.nan), where=counts_ND > 0)
    states = np.full(k, 'COVERAGE_FAIL', dtype=object)
    eligible = (counts_AD >= min_AD) & (counts_ND >= min_ND)
    # Constant beta probes cannot support a meaningful inferential variance.
    spread = np.nanmax(Y, axis=1) - np.nanmin(Y, axis=1)
    constant = np.isfinite(spread) & (spread < 1e-12)
    eligible = eligible & ~constant
    states[constant] = 'CONSTANT_BETA'
    complete = eligible & np.isfinite(Y).all(axis=1)
    if complete.any():
        yc = Y[complete]
        xtx_inv = np.linalg.inv(X.T @ X)
        b = (xtx_inv @ X.T @ yc.T).T
        residual = yc - b @ X.T
        df = nsamples - p
        mse = np.einsum('ij,ij->i', residual, residual) / df
        se2 = mse * xtx_inv[1, 1]
        okay = np.isfinite(se2) & (se2 > 0)
        ix = np.flatnonzero(complete)
        estimates[ix] = b[:, 1]
        states[ix] = np.where(okay, 'OK', 'ZERO_VARIANCE')
        if okay.any():
            valid_ix = ix[okay]
            stderr = np.sqrt(se2[okay])
            t = b[okay, 1] / stderr
            pvals[valid_ix] = 2 * student_t.sf(np.abs(t), df)
            crit = student_t.ppf(0.975, df)
            low[valid_ix] = b[okay, 1] - crit * stderr
            high[valid_ix] = b[okay, 1] + crit * stderr
    # Fit only rows with missing values individually, using QR-based least squares.
    for i in np.flatnonzero(eligible & ~complete):
        valid = np.isfinite(Y[i])
        xv, yv = X[valid], Y[i, valid]
        df = len(yv) - p
        if df <= 0 or np.linalg.matrix_rank(xv) != p:
            states[i] = 'DESIGN_RANK_FAIL'
            continue
        b, _, _, _ = np.linalg.lstsq(xv, yv, rcond=None)
        estimates[i] = b[1]
        residual = yv - xv @ b
        mse = float(residual @ residual) / df
        # QR-factor covariance avoids unstable unguarded pseudo-inversion.
        _, R = np.linalg.qr(xv, mode='reduced')
        rinv = np.linalg.solve(R, np.eye(p))
        se2 = mse * float((rinv @ rinv.T)[1, 1])
        if not np.isfinite(se2) or se2 <= 0:
            states[i] = 'ZERO_VARIANCE'
            continue
        se = math.sqrt(se2)
        pvals[i] = 2 * student_t.sf(abs(b[1] / se), df)
        crit = student_t.ppf(0.975, df)
        low[i], high[i] = b[1] - crit*se, b[1] + crit*se
        states[i] = 'OK'
    return pd.DataFrame({'n_AD': counts_AD, 'n_ND': counts_ND,
                         'mean_beta_AD': means_AD, 'mean_beta_ND': means_ND,
                         'delta_beta_unadjusted': means_AD-means_ND,
                         'adjusted_AD_beta': estimates, 'adjusted_p': pvals,
                         'CI95_low': low, 'CI95_high': high,
                         'status': states})


def analyze(ad_path, nd_path, wa, wn, X, out, chunk, coverage, audit):
    interim = out / 'MTG_genomewide_OLS_intermediate.csv.gz'
    full_result = out / 'MTG_genomewide_age_sex_plate_adjusted.csv.gz'
    significant = out / 'MTG_genomewide_BH_FDR005_absDeltaBeta002.csv'
    n_all, n_noncg, n_cg, invalid_betas = 0, 0, 0, 0
    ids_seen = set()
    written = False
    print('Genome-wide CpG scan started; no network download required.', flush=True)
    with gzip.open(ad_path, 'rt', encoding='utf-8-sig', errors='replace') as fa, \
         gzip.open(nd_path, 'rt', encoding='utf-8-sig', errors='replace') as fn, \
         gzip.open(interim, 'wt', encoding='utf-8', newline='') as fout:
        next(fa); next(fn)
        buffer_ids, buffer_y = [], []
        for lineno, (line_a, line_n) in enumerate(zip(fa, fn), start=2):
            a = line_a.rstrip('\r\n').split('\t')
            n = line_n.rstrip('\r\n').split('\t')
            if len(a) != wa or len(n) != wn:
                raise ValueError(f'Line {lineno}: invalid width {len(a)}, {len(n)}')
            if a[0] != n[0]:
                raise ValueError(f'Line {lineno}: AD/ND probe order mismatch {a[0]}, {n[0]}')
            n_all += 1
            probe = a[0]
            if not CG_RE.fullmatch(probe):
                n_noncg += 1
                continue
            if probe in ids_seen:
                raise ValueError(f'Duplicate CpG ID: {probe}')
            ids_seen.add(probe)
            vals, rejected = parse_beta(a[1::2] + n[1::2])
            invalid_betas += rejected
            buffer_ids.append(probe)
            buffer_y.append(vals)
            n_cg += 1
            if len(buffer_ids) >= chunk:
                batch = tests_block(np.vstack(buffer_y), X, math.ceil(225*coverage), math.ceil(179*coverage))
                batch.insert(0, 'CpG', buffer_ids)
                batch.to_csv(fout, index=False, header=not written)
                written = True
                buffer_ids, buffer_y = [], []
                if n_cg % (chunk*10) == 0:
                    print('CpG models processed:', n_cg, flush=True)
        if fa.readline() or fn.readline():
            raise ValueError('AD and ND files have unequal row counts')
        if buffer_ids:
            batch = tests_block(np.vstack(buffer_y), X, math.ceil(225*coverage), math.ceil(179*coverage))
            batch.insert(0, 'CpG', buffer_ids)
            batch.to_csv(fout, index=False, header=not written)
    if not n_cg:
        raise ValueError('No cg probes found')
    df = pd.read_csv(interim)
    if len(df) != n_cg:
        raise RuntimeError('Results row count mismatch')
    valid = df.adjusted_p.notna()
    if not valid.any():
        raise RuntimeError('No valid adjusted p values')
    df['BH_FDR_genomewide'] = np.nan
    df.loc[valid, 'BH_FDR_genomewide'] = multipletests(df.loc[valid, 'adjusted_p'].to_numpy(), method='fdr_bh')[1]
    df.to_csv(full_result, index=False, compression='gzip')
    hits = df.loc[(df.BH_FDR_genomewide < 0.05) & (df.adjusted_AD_beta.abs() >= 0.02)].copy()
    hits.sort_values('BH_FDR_genomewide').to_csv(significant, index=False)
    audit.update({'total_rows_in_matrix': n_all, 'non_cg_rows_excluded': n_noncg,
                  'cg_rows': n_cg, 'invalid_beta_cells_set_missing': invalid_betas,
                  'p_values_calculated': int(valid.sum()),
                  'coverage_or_model_failures': int((~valid).sum()),
                  'n_BH_FDR_genomewide_lt005': int((df.BH_FDR_genomewide < 0.05).sum()),
                  'n_BH_FDR_genomewide_lt005_and_abs_adjusted_beta_ge002': len(hits),
                  'result_all': full_result.name, 'result_hits': significant.name})
    interim.unlink(missing_ok=True)
    print('Finished:', n_cg, 'cg probes;', int(valid.sum()), 'tested;', len(hits), 'screened hits.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True, help='Verified Module17 MTG_covariates_verified.csv')
    parser.add_argument('--out', type=Path, default=Path('Module18_results'))
    parser.add_argument('--inspect', action='store_true', help='Metadata and matrix headers only')
    parser.add_argument('--analyze', action='store_true', help='Run genome-wide adjusted OLS and BH')
    parser.add_argument('--chunk', type=int, default=500)
    parser.add_argument('--min-coverage', type=float, default=0.8)
    args = parser.parse_args()
    if args.inspect == args.analyze:
        parser.error('Choose exactly one of --inspect or --analyze')
    if args.chunk <= 0 or not (0 < args.min_coverage <= 1):
        parser.error('Invalid --chunk or --min-coverage')
    args.out.mkdir(parents=True, exist_ok=True)
    path_audit = args.out / 'module18_audit.json'
    audit = {'status': 'STARTED', 'started_utc': datetime.now(timezone.utc).isoformat(),
             'dataset': 'GSE134379', 'region': 'MTG',
             'model': 'OLS: methylation beta ~ diagnosis_AD + age_per_10_years + male + plate fixed effects',
             'FDR_scope': 'BH across ALL successfully tested cg probes (genome-wide among analysed cg probes)',
             'limitations': ['Processed beta matrices lack validated detection p values; no detection-p QC.',
                             'No cell-type fraction adjustment; brain cell mixture confounding persists.',
                             'Plate adjusted; bead-chip within plate is not additionally modeled.',
                             'Postmortem independent cohort; not individual-level multiomics correlation.',
                             'No raw-array normalization, cross-reactive-probe or SNP-overlap probe filtering.',
                             'OLS on beta (not M values); associations are not causal.']}
    write_audit(path_audit, audit)
    try:
        md = design_metadata(args.metadata)
        pa, pn = (args.data_dir / FILES[g] for g in ('AD', 'ND'))
        ids_a, wa = header_info(pa, 'AD', md)
        ids_n, wn = header_info(pn, 'ND', md)
        if set(ids_a) & set(ids_n):
            raise ValueError('AD and ND overlap')
        X, design_audit = design(md, ids_a, ids_n)
        audit.update({'n_subjects': len(md), 'n_AD': len(ids_a), 'n_ND': len(ids_n),
                      'n_beta_columns': len(ids_a) + len(ids_n), 'design': design_audit,
                      'minimum_coverage_per_group': args.min_coverage})
        if args.inspect:
            audit['status'] = 'INSPECT_PASS_ANALYSIS_NOT_RUN'
            print('PASS: 404 samples and full-rank design. Analysis NOT run.', flush=True)
        else:
            analyze(pa, pn, wa, wn, X, args.out, args.chunk, args.min_coverage, audit)
            audit['status'] = 'COMPLETED_GENOMEWIDE_EXPLORATORY'
    except Exception as exc:
        audit['status'] = 'FAILED'
        audit['error'] = str(exc)
        audit['traceback'] = traceback.format_exc(limit=6)
        raise
    finally:
        audit['ended_utc'] = datetime.now(timezone.utc).isoformat()
        write_audit(path_audit, audit)
        print('Audit:', path_audit, 'Status:', audit['status'], flush=True)

if __name__ == '__main__':
    main()
