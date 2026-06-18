/* ChromoPainterFold.h - exact block-fold chunk-count (see ChromoPainterFold.c) */
#ifndef CHROMOPAINTERFOLD_H
#define CHROMOPAINTERFOLD_H

/* Numeric constants used by the block fold (ChromoPainterFold.c). cp_emis() is the
   fold's emission; it must match the dense forward-backward emission in
   ChromoPainterSampler.c (the validation harness checks fold == dense). */
#ifndef SMALL_NUM
#define SMALL_NUM 1e-20   /* emission floor for the gap allele (code 8) */
#endif
#ifndef MIN_NE
#define MIN_NE 1e-8       /* lower floor on the N_e (-in) E-M estimate */
#endif

/* emission e(recipient allele r, donor allele d, donor mutation m):
   r==9 -> 1.0 (recipient missing ignored); r==8 -> gap (SMALL_NUM); else (1-m)/m. */
static inline double cp_emis(int r, int d, double m){
    if(r==9) return 1.0;
    if(r==8) return (r==d)?(1-SMALL_NUM):SMALL_NUM;
    return (r==d)?(1-m):m;
}

/* Compute per-population corrected chunk counts, expected differences, expected chunk
   lengths, the per-pop start term, N_e, the regional bootstrap, the forward log-
   likelihood (*out_loglik, = the dense forward's return) and, optionally, -s samples,
   via the O(N*Umean) block fold - identical to the dense forward-backward-chunkcount.
   Grouping is by (local substring, parameter class), so it is exact for uniform
   copy_prob + global mutation and for fixed per-donor copy_prob / mutation.

   out_ccpop[p] is the per-pop chunk count INCLUDING the start term; out_start[p]
   (optional) is the per-pop start term alone, so out_ccpop-out_start is the per-pop
   posterior EXCLUDING start (the dense copy_prob_new total that drives -ip).

   Regional bootstrap (optional, NULL/0 to skip): with out_regfinal/out_regsq (npop)
   and out_numreg, the fold banks the cumulative per-locus total into a region each time
   it crosses region_size (rounding 1e-7) and accumulates the per-pop sums and sums-of-
   squares, dropping the final partial region as the dense does. FP-equivalent (~1e-9),
   NOT bit-identical: a boundary can land one locus off the dense, so the raw per-region
   files may differ - the conserved quantities are the aggregate sum and chromocombine's c.

   Samples (optional): samplesTOT>0 with out_samples!=NULL draws samplesTOT copying paths
   into out_samples[s*nloci+l] (donor index); distributionally identical to the dense. */
void cpfold_perpop(signed char *newh, signed char **existing_h, int nhaps, int nloci,
                   double *TransProb, double *MutProb_vec, double *copy_prob, double *copy_probSTART,
                   double *pos, double *lambda, double delta, double rhobar,
                   int *pop_vec_in, int ndonorpops, int Ustar,
                   double *out_ccpop, double *out_start, double *out_ndiff, double *out_nlen, double *out_Ne,
                   double *out_loglik, double *out_etp, double *out_ecp,
                   double region_size, double *out_regfinal, double *out_regsq, int *out_numreg,
                   int samplesTOT, int *out_samples);

/* Defensive teardown of the engine's static donor buffer (no-op if none). Idempotent. */
void cpfold_cleanup(void);

#endif
