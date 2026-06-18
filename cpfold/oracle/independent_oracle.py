#!/usr/bin/env python3
"""
INDEPENDENT ORACLE for ChromoPainter corrected chunk counts.

Purpose: validate `fs cp` (dense AND -fold) against a reference that shares NO code
with finestructure4 - it implements the Li-Stephens copying HMM from the model
definition and computes the per-pop chunk counts by EXACT brute-force enumeration
over all K^N donor paths. This catches a shared bug in the emission / transition /
chunk-count code (which a fold-vs-dense comparison cannot, since both share cp_emis).

Tiny case (cpfold/oracle/data.*): 1 recipient (TGT) vs 4 donor haplotypes in 2 pops,
6 loci. It exercises BOTH special emission codes with independent code:
  - the missing allele (9) at a recipient locus (TGT[2]) and a donor locus (B2[3]);
  - the gap allele (8) at a recipient locus (TGT[4]) with one matching donor (A1[4]).

Cross-checked by a SECOND, formula-free Monte-Carlo estimator in this same file
(importance sampling from the prior transition model, weighting by the emission
product, and counting chunk starts from the SAMPLED jump indicator rather than the
brute force's analytical jump/(no-jump+jump) split). Both agree with fs:
  brute-force   popA=1.4538194056  popB=0.2377183618   (exact; the harness gates on this)
  monte-carlo   popA~1.45432       popB~0.23953        (4e6 samples; agrees to ~2e-3)
  fs dense/fold popA, popB         (printed, 6 dp)
Regenerate the committed expected.txt from this script:  python3 independent_oracle.py
"""
import numpy as np
from itertools import product

# ---- DATA (cpfold/oracle/data.phase) ----
TGT = [0, 0, 9, 1, 8, 0]          # recipient; 9 = missing, 8 = gap
donors = {
    'A1': [0, 0, 0, 1, 8, 0],     # popA; 8 = gap allele (matches the recipient gap at locus 4)
    'A2': [0, 1, 0, 1, 0, 1],     # popA
    'B1': [1, 1, 1, 0, 1, 0],     # popB
    'B2': [1, 0, 1, 9, 1, 0],     # popB; 9 = missing donor allele
}
donor_names = ['A1', 'A2', 'B1', 'B2']
donor_mat = np.array([donors[n] for n in donor_names])   # (K=4, N=6)
pos = [100, 200, 400, 700, 1100, 1600]
lam = [1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5]

# ---- PARAMETERS (fs cp ... -j -i 0 -n 100 -M 0.01) ----
K, N = 4, 6
rho = 100.0
mut = 0.01
SMALL_NUM = 1e-20                      # emission floor for the gap allele (matches fs)
copy_prob = np.full(K, 1.0 / K)        # jump-target prob = 1/K
copy_probSTART = np.full(K, 1.0 / K)   # locus-0 prior   = 1/K
popA_idx, popB_idx = [0, 1], [2, 3]

# ---- EMISSION e(r,d): identical to the model, re-derived from scratch ----
def emit(r, d):
    if r == 9:    return 1.0                                   # recipient missing => uninformative
    if r == 8:    return (1.0 - SMALL_NUM) if r == d else SMALL_NUM   # gap: matches only a gap donor
    if r == d:    return 1.0 - mut
    return mut

E = np.array([[emit(TGT[l], donor_mat[i][l]) for i in range(K)] for l in range(N)])

# ---- TRANSITION: T[l] over interval (l, l+1] ----
T = np.array([1.0 - np.exp(-(pos[l + 1] - pos[l]) * rho * lam[l]) for l in range(N - 1)])

def trans_prob(j, i, l):   # P(i at l+1 | j at l)
    return (1.0 - T[l]) * (1.0 if i == j else 0.0) + T[l] * copy_prob[i]

# ---- METHOD 1: EXACT BRUTE FORCE over all K^N paths ----
def brute_force():
    Z = 0.0
    chunk_accum = np.zeros(K)
    for path in product(range(K), repeat=N):
        w = copy_probSTART[path[0]] * E[0][path[0]]
        for l in range(1, N):
            w *= trans_prob(path[l - 1], path[l], l - 1) * E[l][path[l]]
        if w == 0.0:
            continue
        Z += w
        starts = np.zeros(K)
        starts[path[0]] += 1.0                         # locus 0 is always a chunk start
        for l in range(1, N):
            j, i = path[l - 1], path[l]
            if i != j:
                starts[i] += 1.0                       # donor change => new chunk
            else:                                      # same donor: split no-jump vs jump-repick
                jump = T[l - 1] * copy_prob[i]
                starts[i] += jump / ((1.0 - T[l - 1]) + jump)
        chunk_accum += w * starts
    return chunk_accum / Z, np.log(Z)

# ---- METHOD 2: formula-free Monte Carlo (independent of the analytical split) ----
# Importance sampling: propose paths from the prior (copy_probSTART then the
# (1-T) stay / T jump-to-copy_prob transition), weight by the emission product, and
# count a chunk start at locus l IFF the proposal sampled a jump entering l. The jump
# indicator comes from the proposal, NOT from the brute force's jump/(no-jump+jump)
# expectation, so a bug in that analytical formula would show up as MC != brute-force.
def monte_carlo(n_samples, seed=0):
    rng = np.random.default_rng(seed)
    path = np.empty((n_samples, N), dtype=np.int64)
    jumped = np.zeros((n_samples, N), dtype=bool)   # jumped[:,l] = recombination entering locus l
    cps, cpj = np.cumsum(copy_probSTART), np.cumsum(copy_prob)
    path[:, 0] = np.searchsorted(cps, rng.random(n_samples))
    w = E[0][path[:, 0]].astype(float).copy()
    for l in range(1, N):
        stay = rng.random(n_samples) < (1.0 - T[l - 1])
        path[:, l] = path[:, l - 1]                 # no jump => copy the previous donor
        jm = ~stay
        if jm.any():
            path[jm, l] = np.searchsorted(cpj, rng.random(int(jm.sum())))
        jumped[:, l] = jm
        w *= E[l][path[:, l]]
    starts = np.zeros((n_samples, K))
    np.add.at(starts, (np.arange(n_samples), path[:, 0]), 1.0)   # locus 0 is always a start
    for l in range(1, N):
        jl = jumped[:, l]
        np.add.at(starts, (np.nonzero(jl)[0], path[jl, l]), 1.0)
    return (w[:, None] * starts).sum(0) / w.sum()

if __name__ == '__main__':
    cc, loglik = brute_force()
    mc = monte_carlo(4_000_000, seed=0)
    popA = sum(cc[i] for i in popA_idx)
    popB = sum(cc[i] for i in popB_idx)
    mcA = sum(mc[i] for i in popA_idx)
    mcB = sum(mc[i] for i in popB_idx)
    print(f"forward_loglik = {loglik:.15f}")
    print("per-donor (brute-force exact):")
    for i, n in enumerate(donor_names):
        print(f"  {n}: {cc[i]:.10f}   (mc {mc[i]:.6f})")
    print(f"popA (A1+A2) = {popA:.10f}   (mc {mcA:.6f})")
    print(f"popB (B1+B2) = {popB:.10f}   (mc {mcB:.6f})")
