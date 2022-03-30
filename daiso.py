
# Degradation Aware ISOform abundance estimation -------------------------------

# Libraries --------------------------------------------------------------------

import os
import sys
import time
import argparse
import numpy as np
import pysam as ps
import pandas as pd
import scipy.sparse as sp
import scipy.special as spf
import matplotlib.pyplot as plt

from tqdm import tqdm
from math import isclose
from scipy.stats import dirichlet, linregress, beta
import scipy.stats as st
from sklearn.preprocessing import normalize
from statsmodels.distributions.empirical_distribution import ECDF

VERSION = '0.0.2'

# Helper functions -------------------------------------------------------------

def get_AS(tags):
    '''
    Args:
    - tags (list of tuple): pysam read.tags
    Returns:
    - alignment score (int)
    '''
    return [i[1] for i in tags if i[0] == 'AS'][0]

def check_unit_interval(value):
    value = float(value)
    if value < 0 or value > 1:
        raise argparse.ArgumentTypeError('Degradation rate must be in the open interval (0,1)')
    return value

def flex_read_len_prob(args, curr_read_len, curr_iso_len, ecdf):
    '''
    Args:
    - read_length (float) in kb
    - isoform_length (float) in kb
    - ecdf (ECDF)
    - delta (float) in kb for approximating probabilities
    Return:
    - probability of observing read given isoform
    '''
    max_len = max(ecdf.x)
    delta = args.delta / 1e3
    full_len_tol = float(args.full_len_tol) / 1e3
    # Case 0: read length greater than isoform length - 0 prob
    if curr_read_len > curr_iso_len:
        return 0
    # Case 1: read length == isoform length
    full_len_prob = max(1 - ecdf(curr_iso_len), 1 - ecdf(max_len - delta))
    if isclose(curr_read_len, curr_iso_len, abs_tol = full_len_tol):
        return full_len_prob
    # Case 2: read length < isoform length
    if curr_read_len < max_len:
        lim_low = ecdf(max(0,curr_read_len - delta))
        lim_upp = ecdf(curr_read_len + delta)
        degraded_prob_unnorm = lim_upp - lim_low
    else:
        degraded_prob_unnorm = 1 - ecdf(max_len - 1e-3)
    return degraded_prob_unnorm / ecdf(curr_iso_len)

def exact_read_len_prob(args, curr_read_len, curr_iso_len, avg_deg_rate):
    '''
    Args:
    - read_id (int): read id
    - iso_id (int): isoform id
    - read_len (dict): read length map read_id:read_length
    - iso_len (dict): isoform length map iso_id:iso_length
    - iso_d_rate (np array (float): degradation rates for isoforms
    - max_len (float): length in kb
    '''
    max_len = 1 / avg_deg_rate
    full_len_tol = float(args.full_len_tol) / 1e3
    if isclose(curr_read_len, curr_iso_len, abs_tol = full_len_tol): # full length read
        return 1 - avg_deg_rate * curr_read_len
    elif curr_read_len < min(max_len, curr_iso_len): # degraded read
        #return avg_deg_rate/ (curr_iso_len * 1e3)
        return avg_deg_rate / 1e3
    return 0

def parser_init():
    '''
    Description: Initialize a parser
    Returns: argparse parser
    '''
    parser = argparse.ArgumentParser(
        prog='daiso',
        usage='%(prog)s -a alignment.bam [options]',
        description = 'Degradation-aware isoform quantification',
        add_help=False
        )

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    argset1 = parser.add_argument_group('Parameters for alignment filtering')
    argset2 = parser.add_argument_group('Parameters for degradation estimation')
    argset3 = parser.add_argument_group('Parameters for inference')
    general = parser.add_argument_group('General help')

    required.add_argument(
        '-a', '--alignment',
        help='BAM file (sorted and indexed) containing reads aligned to the transcriptome'
    )

    general.add_argument(
        '-h', '--help', action='help', default=argparse.SUPPRESS,
        help='Show this help message and exit'
    )
    general.add_argument(
        '-v', '--version', action='store_true',
        help='Print version'
    )

    optional.add_argument(
        '-o', '--output', default='out',
        help='Output prefix. The count and degradation plot files will be written out to <output>_counts.tsv and <output>_survival.png respectively (default: out)'
    )
    optional.add_argument(
        '--no_quant', action='store_true',
        help='Estimate degradation only. No quantification is run. (default: False)'
    )
    optional.add_argument(
        '--seed', default=42, type=int,
        help='Set seed (default: 42)'
    )

    argset1.add_argument(
        '--filter_distance', default=50, type=float,
        help='Filter alignments where end position is further than this distance away from the annotated 3 prime end (default: 50 bp)'
    )
    argset1.add_argument(
        '--filter_score', default=0.95, type=float,
        help='Filter alignments where the alignment score is lower than this fraction of the best alignment score for this read (default: 0.95)'
    )

    argset2.add_argument(
        '--deg_rate', type=check_unit_interval,
        help='Degradation rate in (0,1). If this argument is supplied, degradation rate will not be estimated from the data and deg_const will be True'
    )
    argset2.add_argument(
        '--deg_const', action='store_true',
        help='Uses exact constant degradation model for read length-isoform agreement (default:False)'
    )
    argset2.add_argument(
        '--bin_size', default=500, type=float,
        help='Bin isoform lengths into bins of this size (default: 500 bp)'
    )
    argset2.add_argument(
        '--min_read_count', default=5, type=int,
        help='Only estimate degradation with isoforms with at least this many reads (default: 5)'
    )
    argset2.add_argument(
        '--min_iso_count', default=1, type=int,
        help='Only estimate degradation with length bins supported by this many isoforms (default: 1)'
    )
    argset2.add_argument(
        '--full_len_tol', default=50, type=float,
        help='Consider a read full length if it is within this many bp of the annotated isoform length (default: 50 bp)'
    )
    argset2.add_argument(
        '--delta', default=50, type=float,
        help='Left and right shift for estimating exact probabilities for ecdf (default: 50 bp)'
    )
    argset2.add_argument(
        '--return_survival', action='store_false',
        help='Return the survival function (1-ecdf). Writes values to <output>_survival.tsv (default: True)'
    )

    argset3.add_argument(
        '--inference', default='EM', choices=['EM','VB'],
        help='Inference algorithm - expectation maximization or variational Bayesian inference (default: EM)'
    )
    argset3.add_argument(
        '--max_iter', default=200, type=int,
        help='Maximum number of iterations (default: 200)'
    )
    argset3.add_argument(
        '--return_loglik', action='store_false',
        help='Return the log likelihood over iterations for EM. Writes values to <output>_loglik.tsv (default: True)'
    )
    argset3.add_argument(
        '--prior', default='gamma_hyper', choices=['symmetric', 'gamma_hyper'],
        help='Choice of prior - symmetric prior or a gamma hyperprior (default: gamma_hyper)'
    )
    argset3.add_argument(
        '--alpha_zero', default=1, type=float,
        help='Value of the concentration parameter for the symmetric prior (default: 1)'
    )
    argset3.add_argument(
        '--gamma_rate', default=5, type=float,
        help='Rate parameter for the gamma hyperprior (default: 5)'
    )
    argset3.add_argument(
        '--gamma_scale', default=5, type=float,
        help='Scale parameter for the gamma hyperprior (default: 5)'
    )
    argset3.add_argument(
        '--return_cred', action='store_true',
        help='Return credible interals for parameters (default: False)'
    )
    argset3.add_argument(
        '--cred_int', default=0.95, type=float,
        help='Width of credible interval (default: 0.95)'
    )

    return parser

def process_reads(args):
    # Extract from SAM:
    # - aligned read lengths <scriptl> from SAM file
    # - observed alignment <alpha> between read and transcript isoform
    # - isoform lengths
    alnm_file = args.alignment
    samfile = ps.AlignmentFile(alnm_file, 'rb')
    # Map from reference name to reference length
    ref_len_map = dict(zip(samfile.references, samfile.lengths))
    # Map from read query name from pysam (str) : max AlignmentScore (int)
    read_score = dict()

    # Step 1 - Record maximum alignment score ----------------------------------
    for read in tqdm(samfile.fetch()):
        read_nm = read.query_name
        score = get_AS(read.tags)
        if read_nm not in read_score:
            read_score[read_nm] = score
        else:
            if score > read_score[read_nm]:
                read_score[read_nm] = score

    # Map from reference id from pysam (int) : isoform id (int)
    iso_map = dict()
    # Map from isoform id (int) : isoform name (str)
    iso_nm = dict()
    # Map from isoform id (int) : isoform len (int)
    iso_len = dict()
    # Map from read query name from pysam (str) : read_id (int)
    read_map = dict()
    # Map from read_id (int) : read_len (int)
    read_len = dict()
    # alpha matrix - map from read_id (int) : compatible isoforms with score fraction (list(tuple(int, float)))
    alpha_dict = dict()

    # Step 2 - Parse and filter reads and construct alignment ------------------
    # Store read lengths
    read_counter = 0
    read_id = 0
    iso_id = 0
    for read in tqdm(samfile.fetch()):
        read_nm = read.query_name
        # Check if read passes filters
        curr_score = get_AS(read.tags)
        fail_score = curr_score < args.filter_score * read_score[read_nm]
        fail_strand = read.is_reverse
        fail_map = read.is_unmapped
        fail_pos = ref_len_map[read.reference_name] - read.reference_end > args.filter_distance
        if (fail_score or fail_strand or fail_map or fail_pos):
            continue
        # Encode read query name with integers read_id
        if read_nm not in read_map:
            read_map[read_nm] = read_id
            read_len[read_id] = read.query_alignment_length
            read_id += 1
        # Encode reference id with integers iso_id
        iso_ref_id = read.reference_id
        if iso_ref_id not in iso_map:
            iso_map[iso_ref_id] = iso_id
            iso_len[iso_id] = ref_len_map[read.reference_name]
            iso_nm[iso_id] = read.reference_name
            iso_id += 1
        # Record observed alignment alpha
        if read_map[read_nm] not in alpha_dict:
            alpha_dict[read_map[read_nm]] = [
                (iso_map[iso_ref_id], curr_score / read_score[read_nm])
            ]
        else:
            alpha_dict[read_map[read_nm]].append(
                (iso_map[iso_ref_id], curr_score / read_score[read_nm])
            )
        # Increment counter
        read_counter += 1

    # Get number of reads and number of mapped isoforms
    num_read = len(read_map)
    num_iso = len(iso_map)
    print('Number of reads: {}'.format(num_read))
    print('Number of isoforms: {}'.format(num_iso))

    # Construct alpha
    alpha_mat = sp.dok_matrix((num_read, num_iso), dtype = np.float32)
    for read_id, alnms in tqdm(alpha_dict.items()):
        for alnm in alnms:
            alpha_mat[read_id, alnm[0]] = alnm[1]

    alpha_dok = alpha_mat
    alpha_mat = alpha_mat.tocsc()

    return {
        'alpha_mat':alpha_mat,
        'alpha_dok':alpha_dok,
        'iso_len':iso_len,
        'read_len':read_len,
        'num_read':num_read,
        'num_iso':num_iso,
        'iso_nm':iso_nm
    }

def estimate_degradation(args, alpha_mat, iso_len, read_len):
    # Step 3 - Estimate global deg. rate -------------------------------------------
    # Extract train data - find uniquely-mapped-to isoforms
    alpha_mat_t = alpha_mat.transpose()
    uniq_iso_mat = alpha_mat_t @ alpha_mat
    # Get counts for iso
    uniq_iso_counts = np.array(uniq_iso_mat.sum(axis=0)).flatten()
    # Get isos that appear only once
    uniq_iso_nz = uniq_iso_mat.nonzero()[0]
    uniq_iso_id,uniq_iso_rep = np.unique(uniq_iso_nz, return_counts= True)
    # Apply count filter
    uniq_iso_keep = np.logical_and(uniq_iso_rep == 1, uniq_iso_counts > args.min_read_count)
    uniq_iso = uniq_iso_id[uniq_iso_keep]
    uniq_iso_lens = np.array([iso_len[i] for i in uniq_iso])
    # Now, we sample isoforms evenly across lengths
    uniq_iso_lens = np.round(uniq_iso_lens / args.bin_size) * args.bin_size
    iso_sample = pd.DataFrame(zip(uniq_iso, uniq_iso_lens), columns = ['iso_id', 'iso_len'])
    # Determine sample size
    iso_group_size = iso_sample.groupby('iso_len', group_keys=False).size()
    iso_group_size = iso_group_size[iso_group_size > args.min_iso_count]
    print('Isoform count by length bin:\n {}'.format(iso_group_size))
    iso_len_max = max(iso_group_size.index)
    iso_sample = iso_sample[iso_sample['iso_len'] <= iso_len_max]
    length_sample_size = np.round(np.median(iso_group_size.to_numpy())).astype(int)
    print('Sampling {} isoforms per length bin'.format(length_sample_size))
    iso_sample = iso_sample.groupby('iso_len', group_keys=False).apply(
        lambda x:x.sample(n=length_sample_size, random_state=args.seed) if x.shape[0] >= length_sample_size else x
    )
    uniq_iso = iso_sample.iso_id.to_numpy()
    alpha_uniq = alpha_mat[:,uniq_iso] # Subset to unique isoforms
    # Train only on degraded reads
    train_data = []
    train_data_iso = []
    # Get the index of reads and shifted isoforms
    read_idx_arr, iso_idx_arr = alpha_uniq.nonzero()
    for i in tqdm(range(len(read_idx_arr))):
        # Full length read, ignore
        read_idx = read_idx_arr[i]
        iso_idx = uniq_iso[iso_idx_arr[i]]
        if isclose(read_len[read_idx], iso_len[iso_idx], abs_tol = args.full_len_tol / 1e3):
            continue
        train_data.append(read_len[read_idx])
        train_data_iso.append(iso_idx)

    print('Estimating degradation from {} reads from {} isoforms'.format(
        len(train_data), len(np.unique(train_data_iso))
    ))

    # Step 4 - Density estimation --------------------------------------------------
    ecdf = ECDF([x/1e3 for x in train_data])
    # Report the average estimated degradation
    ex = ecdf.x
    ex[ex < 0] = 0
    ey = ecdf.y
    d_rate_avg_est = round(linregress(ex,ey).slope,3)

    print('Average estimated degradation rate: {} ncov / kb'.format(d_rate_avg_est))
    plt.figure(figsize = (10,5))
    plt.plot(ecdf.x,1-ecdf.y)
    plt.title('Survival function on training reads with d = {}'.format(d_rate_avg_est))
    plt.savefig('%s_survival.png'%args.output)
    plt.close()

    survival = pd.DataFrame(
        {'x': ex,
         'ncov': 1-ecdf.y
        }
    )
    survival.to_csv('%s_survival.tsv'%args.output, sep = '\t', index = False)
    return {
        'ecdf': ecdf,
        'avg_deg_rate': d_rate_avg_est
    }

def length_model(args, alpha_dok, iso_len, read_len, ecdf, avg_deg_rate):
    # Construct rho
    num_read = len(read_len)
    num_iso = len(iso_len)
    rho_mat = sp.dok_matrix((num_read, num_iso), dtype = np.float32)
    for read_id, iso_id in tqdm(alpha_dok.keys()):
        curr_read_len = read_len[read_id] / 1e3
        curr_iso_len = iso_len[iso_id] / 1e3
        if args.deg_const:
            rho_mat[read_id, iso_id] = exact_read_len_prob(
                args, curr_read_len, curr_iso_len, avg_deg_rate
            )
        else:
            rho_mat[read_id, iso_id] = flex_read_len_prob(
                args, curr_read_len, curr_iso_len, ecdf
            )
    # Sparse column matrix
    rho_mat = rho_mat.tocsc()
    return rho_mat

def run_em(args, rho_mat, alpha_mat, num_read, num_iso, iso_nm):
    theta_vec = sp.csc_matrix(np.ones(num_iso) / num_iso)
    # Step 5 - Run EM --------------------------------------------------------------
    em_iter = 0
    converged = False
    log_lik = []
    while not converged:
        print('EM iteration {}'.format(em_iter))
        # E step: ------------------------------------------------------------------
        # row-wise multiply alpha_mat with theta
        gamma = alpha_mat.multiply(theta_vec)
        # element-wise multiply with rho
        gamma = gamma.multiply(rho_mat)
        gamma = normalize(gamma, norm = 'l1', axis = 1)
        # Re-estimation ------------------------------------------------------------
        # M step: ------------------------------------------------------------------
        theta_update = gamma.sum(axis = 0) / num_read
        theta_diff = np.linalg.norm(theta_update - theta_vec)
        print('Theta diff: {}'.format(theta_diff))
        theta_vec = theta_update
        # Calculate log-likelihood: ------------------------------------------------
        ll = alpha_mat.multiply(theta_vec)
        ll = ll.multiply(rho_mat).sum(axis = 1)
        ll = np.where(ll > 0., ll, 1e-5)
        ll = np.sum(np.log(ll))
        print('Log likelihood: {}\n'.format(ll))
        log_lik.append(ll)
        # Convergence on iterations: -----------------------------------------------
        em_iter += 1
        if em_iter >= args.max_iter:
            converged = True

    # Step 6 - Collate and write out results ---------------------------------------
    iso_counts = (theta_vec * num_read).tolist()[0]
    iso_out = []
    for i in range(num_iso):
        iso_out.append(iso_nm[i])

    count_table = pd.DataFrame(
        {'transcript_id': iso_out,
         'estimated_count': iso_counts
        }
    )
    count_table.to_csv('%s_counts.tsv'%args.output, sep = '\t', index = False)

    with open('%s_loglik.tsv'%args.output, 'w') as f:
        for x in log_lik:
            f.write("%s\n" % x)

def run_vb(args, rho_mat, alpha_mat, num_read, num_iso, iso_nm):

    if args.prior == 'symmetric':
        # Draw theta from symmetric Dirichlet with conc. theta_alpha
        alpha_zero = np.ones(num_iso) * args.alpha_zero
    elif args.prior == 'gamma_hyper':
        # Draw alphas from gamma prior with mean 1
        alpha_zero = st.gamma.rvs(
            a = args.gamma_rate, scale = 1 / args.gamma_scale,
            size = num_iso, random_state = args.seed
        )
    alpha_vec = alpha_zero
    theta_vec = sp.csc_matrix(
        dirichlet.rvs(alpha_zero, size = 1, random_state = args.seed)[0]
    )

    # Run VBEM
    vb_iter = 0
    converged = False
    log_lik = []
    while not converged:
        print('VB iteration {}'.format(vb_iter))
        # VBE step: ----------------------------------------------------------------
        gamma = alpha_mat.multiply(rho_mat)
        # Possible numerical underflow
        phi = sp.csc_matrix(
            np.exp(spf.digamma(alpha_vec) - spf.digamma(alpha_vec.sum()))
        )
        gamma = gamma.multiply(phi)
        # Posterior responsibilities
        gamma = normalize(gamma, norm = 'l1', axis = 1)
        # VBM step: ----------------------------------------------------------------
        # Update the parameters of Dirichlet
        alpha_vec = alpha_zero + gamma.sum(axis = 0)
        # Theta: -------------------------------------------------------------------
        theta_exp = alpha_vec / (alpha_zero.sum() + num_read)
        theta_update = theta_exp * num_read
        theta_diff = np.linalg.norm(theta_update - theta_vec)
        print('Theta diff: {}'.format(theta_diff))
        theta_vec = theta_update
        # Convergence: -------------------------------------------------------------
        vb_iter += 1
        if vb_iter >= args.max_iter:
            converged = True

    iso_counts = theta_vec.tolist()[0]
    iso_out = []
    for i in range(num_iso):
        iso_out.append(iso_nm[i])

    count_table = pd.DataFrame(
        {'transcript_id': iso_out,
         'estimated_count': iso_counts
        }
    )

    # 95% credible intervals for each theta
    if args.return_cred:
        cred_int = []
        alpha_vec = alpha_vec.tolist()[0]
        for i in range(len(alpha_vec)):
            cred_int.append(beta.ppf([(1-args.cred_int)/2, (1+args.cred_int)/2],
                alpha_vec[i], sum(alpha_vec) - alpha_vec[i]))

        cred_int = np.array(cred_int)
        theta_low = cred_int[:,0]
        theta_high = cred_int[:,1]
        iso_low = (theta_low * num_read).tolist()
        iso_high = (theta_high * num_read).tolist()
        count_table['counts_low'] = iso_low
        count_table['counts_upp'] = iso_high

    count_table.to_csv('%s_counts.tsv'%args.output, sep = '\t', index = False)

def main():
    parser = parser_init()
    args = parser.parse_args()

    if args.version:
        print('v%s' % VERSION)
        sys.exit()

    if args.alignment is None:
        print('Please provide a path to an alignment file via --alignment')
        sys.exit()

    print(args)
    read_stats = process_reads(args)

    alpha_mat = read_stats['alpha_mat']
    alpha_dok = read_stats['alpha_dok']
    iso_len = read_stats['iso_len']
    read_len = read_stats['read_len']
    num_read = read_stats['num_read']
    num_iso = read_stats['num_iso']
    iso_nm = read_stats['iso_nm']

    if args.deg_rate is not None:
        # If the degradation is specified, don't estimate, use exact.
        print('Degradation provided. Using exact model.')
        args.deg_const = True
        ecdf = None
        avg_deg_rate = args.deg_rate
    else:
        # If the degradation is not specified and deg_const is False, use the emp.
        # If the degradation is not specified and deg_const is True, use the avg deg rate and exact.
        print('Degradation not provided. Estimating degradation.')
        deg_stats = estimate_degradation(args, alpha_mat, iso_len, read_len)
        ecdf = deg_stats['ecdf']
        avg_deg_rate = deg_stats['avg_deg_rate']

    if args.no_quant:
        # Just estimate degradation rate
        sys.exit()

    rho_mat = length_model(args, alpha_dok, iso_len, read_len, ecdf, avg_deg_rate)

    if args.inference == 'EM':
        run_em(args, rho_mat, alpha_mat, num_read, num_iso, iso_nm)

    if args.inference == 'VB':
        run_vb(args, rho_mat, alpha_mat, num_read, num_iso, iso_nm)

if __name__ == '__main__':
    main()
