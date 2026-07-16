/* ChromoPainterFold.c - exact block-fold of ChromoPainter's chunk-count.
 *
 * Produces the per-population corrected_chunk_count identically to the dense
 * forward-backward-chunkcount, in O(N*Umean) instead of O(N*K). Also returns the
 * expected differences, expected chunk lengths, the per-pop start term, the N_e
 * (-in) and global-mutation (-iM) EM quantities, the regional bootstrap, the
 * forward log-likelihood, and (optionally) -s copying-path samples, so -fold drives
 * the full -i N -ip -im -in -iM EM loop and -s. Exact for uniform and for fixed
 * per-donor copy_prob (-p) / mutation (-m). Reuses the caller's exact TransProb.
 * cpfold_perpop() is the entry point, called from sampler() under -fold.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "ChromoPainterFold.h"

/* fail-fast allocators: the engine's arrays are O(N*K) (gigabytes at chr scale),
   so abort with a clear message rather than deref a NULL on OOM. */
static void *cf_xmalloc(size_t n){ void *p=malloc(n); if(!p){ fprintf(stderr,"ChromoPainterFold: OOM (malloc %zu)\n",n); exit(1);} return p; }
static void *cf_xcalloc(size_t n,size_t s){ void *p=calloc(n,s); if(!p){ fprintf(stderr,"ChromoPainterFold: OOM (calloc %zu x %zu)\n",n,s); exit(1);} return p; }
static void *cf_xrealloc(void *q,size_t n){ void *p=realloc(q,n); if(!p){ free(q); fprintf(stderr,"ChromoPainterFold: OOM (realloc %zu)\n",n); exit(1);} return p; }

/* engine globals (set per call by cpfold_perpop; single-threaded engine) */
static int K, N, npop;
/* per-donor copy_prob / copy_probSTART / mutation: one global value in the common
   case, per-donor under -p / -m. The engine groups state by (substring,param-class)
   so all three are constant within each group. */
static double *cf_cp, *cf_cps, *cf_mut;
static int *cf_pclass, cf_nclass;  /* [K] dense id of distinct (cp,cps,mut) tuples; seeds the grouping */
static uint8_t *donors;     /* [N*K] locus-major, built from existing_h per call */
static int *pop_vec;        /* [K] donor population (indexes the per-pop outputs) */
static double *T;           /* [N-1] transition (the caller's TransProb) */
static double *pos_g, *lam_g, delta_g, rho_g;  /* for the N_e (-in) EM update */

/* emission is cp_emis() from ChromoPainterFold.h - shared with the dense FB. */

typedef struct { int s,e,U; } Block;
typedef struct {
    int nb, Umax;
    Block *blk;
    int *gidB;
    int *sizeg,*repg;
    int *cellB, *celloff, *cellGb, *cellG2, maxcell;
    int *gperm, *goff;   /* group->donor CSR per block (for -fold sampling): block b
                            group g's donors are gperm[b*K + goff[b*(Umax+1)+g ..+1]). */
} Groups;

/* panel-fixed (g_b,g_{b+1}) contingency for the L2 boundary join-fold */
static void build_contingency(Groups *G){
    int nb=G->nb, Umax=G->Umax; int *gidB=G->gidB;
    G->celloff=cf_xmalloc(sizeof(int)*(size_t)(nb+1));
    G->cellB=cf_xmalloc(sizeof(int)*(size_t)nb*K);
    int *cmap=cf_xmalloc(sizeof(int)*(size_t)Umax*Umax);
    int *cstamp=cf_xcalloc((size_t)Umax*Umax,sizeof(int)); int gen=0;
    size_t capcell=(size_t)nb*8+16, total=0; int maxcell=0;
    int *cGb=cf_xmalloc(sizeof(int)*capcell), *cG2=cf_xmalloc(sizeof(int)*capcell);
    G->celloff[0]=0;
    for(int b=0;b<nb-1;b++){
        gen++; int nc=0; int *gA=gidB+(size_t)b*K, *gB=gidB+(size_t)(b+1)*K; int *cb=G->cellB+(size_t)b*K;
        for(int i=0;i<K;i++){ int key=gA[i]*Umax+gB[i];
            if(cstamp[key]!=gen){ cstamp[key]=gen; cmap[key]=nc;
                if(total+nc>=capcell){ capcell*=2; cGb=cf_xrealloc(cGb,sizeof(int)*capcell); cG2=cf_xrealloc(cG2,sizeof(int)*capcell); }
                cGb[total+nc]=gA[i]; cG2[total+nc]=gB[i]; nc++; }
            cb[i]=cmap[key]; }
        total+=nc; G->celloff[b+1]=total; if(nc>maxcell)maxcell=nc;
    }
    G->celloff[nb]=total;
    G->cellGb=cGb; G->cellG2=cG2; G->maxcell=maxcell<1?1:maxcell;
    free(cmap); free(cstamp);
}
static void free_groups(Groups *G){ free(G->blk);free(G->gidB);free(G->sizeg);free(G->repg);
    free(G->cellB);free(G->celloff);free(G->cellGb);free(G->cellG2);
    free(G->gperm);free(G->goff); }

/* PBWT adaptive variable-length blocks: grow each block by incrementally
   splitting the substring grouping column-by-column; cut when U reaches Ustar.
   16-way allele code (al&15) keeps all cp allele symbols 0-5,8,9 distinct. */
static Groups build_groups_adaptive(int Ustar){
    Groups G; G.nb=0;
    int cap_blk=64; G.blk=cf_xmalloc(sizeof(Block)*cap_blk); int *gidB=NULL;
    int *gid=cf_xmalloc(sizeof(int)*K), *ng=cf_xmalloc(sizeof(int)*K);
    /* map indexed by ek=gid[i]*16+code; gid[i] < K always (<= K groups), and the
       transient group count between block-start cuts can exceed Ustar on
       multi-allelic data, so size by K (not Ustar) to be overflow-proof. */
    int mapsz=K*16+16;
    int *mp=cf_xmalloc(sizeof(int)*mapsz), *mstamp=cf_xcalloc(mapsz,sizeof(int)); int gen=0, Umax=0;
    /* seed each block by param-class so every group is parameter-homogeneous
       (copy_prob/mutation group-constant); nclass=1 (uniform) gives substring-only. */
    int a=0; for(int i=0;i<K;i++) gid[i]=cf_pclass[i]; int U=cf_nclass; int b=0;
    while(b<N){
        gen++; int newU=0; const uint8_t *col=donors+(size_t)b*K;
        for(int i=0;i<K;i++){ int code=col[i]&15; int ek=gid[i]*16+code;
            if(mstamp[ek]!=gen){ mstamp[ek]=gen; mp[ek]=newU++; } ng[i]=mp[ek]; }
        if(newU>Ustar && (b-a)>=2 && b<N-1){
            if(G.nb>=cap_blk){ cap_blk*=2; G.blk=cf_xrealloc(G.blk,sizeof(Block)*cap_blk); }
            G.blk[G.nb].s=a; G.blk[G.nb].e=b; G.blk[G.nb].U=U;
            gidB=cf_xrealloc(gidB,(size_t)(G.nb+1)*K*sizeof(int)); memcpy(gidB+(size_t)G.nb*K,gid,K*sizeof(int));
            if(U>Umax)Umax=U; G.nb++;
            a=b; for(int i=0;i<K;i++) gid[i]=cf_pclass[i]; U=cf_nclass;
        } else { int *t=gid;gid=ng;ng=t; U=newU; b++; }
    }
    if(G.nb>=cap_blk){ cap_blk*=2; G.blk=cf_xrealloc(G.blk,sizeof(Block)*cap_blk); }
    G.blk[G.nb].s=a; G.blk[G.nb].e=N; G.blk[G.nb].U=U;
    gidB=cf_xrealloc(gidB,(size_t)(G.nb+1)*K*sizeof(int)); memcpy(gidB+(size_t)G.nb*K,gid,K*sizeof(int));
    if(U>Umax)Umax=U; G.nb++;
    free(gid);free(ng);free(mp);free(mstamp);
    G.gidB=gidB; G.Umax=Umax;
    G.sizeg=cf_xmalloc(sizeof(int)*(size_t)G.nb*Umax); G.repg=cf_xmalloc(sizeof(int)*(size_t)G.nb*Umax);
    for(int bb=0;bb<G.nb;bb++){ int Ub=G.blk[bb].U; int *gb=gidB+(size_t)bb*K;
        int *sz=G.sizeg+(size_t)bb*Umax,*rp=G.repg+(size_t)bb*Umax;
        for(int g=0;g<Ub;g++){sz[g]=0;rp[g]=-1;}
        for(int i=0;i<K;i++){int g=gb[i];sz[g]++; if(rp[g]<0)rp[g]=i;}
    }
    /* group->donor CSR (gperm/goff), panel-fixed, for the -fold path sampler. */
    G.goff=cf_xmalloc(sizeof(int)*(size_t)G.nb*(Umax+1)); G.gperm=cf_xmalloc(sizeof(int)*(size_t)G.nb*K);
    { int *cur=cf_xmalloc(sizeof(int)*(Umax+1));
      for(int bb=0;bb<G.nb;bb++){ int Ub=G.blk[bb].U; int *gb=gidB+(size_t)bb*K;
        int *sz=G.sizeg+(size_t)bb*Umax, *go=G.goff+(size_t)bb*(Umax+1), *gp=G.gperm+(size_t)bb*K;
        go[0]=0; for(int g=0;g<Ub;g++){ go[g+1]=go[g]+sz[g]; cur[g]=go[g]; }
        for(int i=0;i<K;i++){ int g=gb[i]; gp[cur[g]++]=i; } }
      free(cur); }
    build_contingency(&G);
    return G;
}

/* FOLDED forward+backward+chunkcount -> per-pop chunk counts (ccpop[npop]).
   Direct port of cpfold.c fold_cc; rh = the recipient allele row. */
/* etp_out[N-1] (per-locus transition prob, for -d) and ecp_out[N*npop] (per-locus
   per-pop copy posterior, for -b) are optional: filled only when non-NULL. */
/* Sample one donor ~ rescaled forward a[l][.] WITHOUT materializing the length-K
   vector: draw a group g prop to (size_g*GF + PF*Sentry_g), then a donor within g
   prop to (GF + PF*a_entry(i)) over that group's members (gperm/goff CSR). gw is a
   >=Umax scratch buffer. Uses the global rand() (seeded by -S). */
static int fold_sample_donor(int l, Groups *Gr, const double *GF, const double *PF,
                             const double *AENTRY, const double *SentS, const size_t *off,
                             const int *locblk, double *gw){
    int Umax=Gr->Umax, b=locblk[l], sb=Gr->blk[b].s, U=Gr->blk[b].U;
    const double *GFr=GF+off[b]+(size_t)(l-sb)*U, *PFr=PF+off[b]+(size_t)(l-sb)*U;
    const double *ent=AENTRY+(size_t)b*K, *Sx=SentS+(size_t)b*Umax;
    const int *sz=Gr->sizeg+(size_t)b*Umax, *go=Gr->goff+(size_t)b*(Umax+1), *gp=Gr->gperm+(size_t)b*K;
    double gtot=0.0;
    for(int g=0;g<U;g++){ double m=sz[g]*GFr[g]+PFr[g]*Sx[g]; if(m<0.0)m=0.0; gw[g]=m; gtot+=m; }
    if(gtot<=0.0) return gp[go[0]];   /* deep underflow: every group mass floored to 0 */
    double u=((double)rand()/RAND_MAX)*gtot, c=0.0; int g=0;
    for(g=0;g<U;g++){ c+=gw[g]; if(u<=c) break; } if(g>=U) g=U-1;
    int lo=go[g], hi=go[g+1]; double dtot=0.0;
    for(int t=lo;t<hi;t++){ double m=GFr[g]+PFr[g]*ent[gp[t]]; if(m>0.0)dtot+=m; }
    double u2=((double)rand()/RAND_MAX)*dtot, c2=0.0; int sel=gp[hi-1];
    for(int t=lo;t<hi;t++){ double m=GFr[g]+PFr[g]*ent[gp[t]]; if(m>0.0)c2+=m; if(u2<=c2){ sel=gp[t]; break; } }
    return sel;
}

/* samplesTOT>0 + sout!=NULL: also draw samplesTOT copying paths (hierarchical FFBS
   over the affine forward, no full Alphamat) into sout[s*N+l] = donor index. */
static void fold_cc(const uint8_t *rh, double *ccpop, double *startpop, double *ndiff, double *nlen, double *Ne_out, double *loglik_out, double *etp_out, double *ecp_out, double region_size, double *out_regfinal, double *out_regsq, int *out_numreg, int samplesTOT, int *sout, Groups *Gr){
    #define EM(L,II) cp_emis(rh[(L)], donors[(size_t)(L)*K+(II)], cf_mut[(II)])
    int nb=Gr->nb, Umax=Gr->Umax; Block *blk=Gr->blk; int *gidB=Gr->gidB;
    int *sizeg=Gr->sizeg, *repg=Gr->repg;
    /* per-(block,group) copy_prob / copy_probSTART (group-constant: each group is
       parameter-homogeneous; = the uniform scalar in the common case). */
    double *CPG=cf_xmalloc((size_t)nb*Umax*sizeof(double)), *CPSG=cf_xmalloc((size_t)nb*Umax*sizeof(double));
    size_t *off=cf_xmalloc(sizeof(size_t)*nb), tot=0;
    for(int b=0;b<nb;b++){ off[b]=tot; tot += (size_t)(blk[b].e-blk[b].s)*blk[b].U; }
    double *GF=cf_xmalloc(tot*sizeof(double)), *PF=cf_xmalloc(tot*sizeof(double)), *Eg=cf_xmalloc(tot*sizeof(double));
    double *AENTRY=cf_xmalloc(sizeof(double)*(size_t)nb*K), *WB=cf_xmalloc(sizeof(double)*(size_t)nb*K);
    double *CEXIT=cf_xmalloc(sizeof(double)*(size_t)nb*K);
    double *As=cf_xmalloc(sizeof(double)*N), *Bs=cf_xmalloc(sizeof(double)*N);
    double *SentS=cf_xmalloc((size_t)nb*Umax*sizeof(double)), *SwS=cf_xmalloc((size_t)Umax*sizeof(double));   /* SwS: per-block scratch */
    double *G=cf_xmalloc(Umax*sizeof(double)), *P=cf_xmalloc(Umax*sizeof(double));
    double *GBp=cf_xmalloc(Umax*sizeof(double)), *PBp=cf_xmalloc(Umax*sizeof(double));
    double *GBc=cf_xmalloc(Umax*sizeof(double)), *PBc=cf_xmalloc(Umax*sizeof(double));
    double *szP=cf_xmalloc((size_t)Umax*npop*sizeof(double)), *SaP=cf_xmalloc((size_t)Umax*npop*sizeof(double));
    double *SwP=cf_xmalloc((size_t)Umax*npop*sizeof(double)), *SawP=cf_xmalloc((size_t)Umax*npop*sizeof(double));
    double *Mc=cf_xmalloc((size_t)Gr->maxcell*npop*sizeof(double)), *Mce=cf_xmalloc((size_t)Gr->maxcell*npop*sizeof(double));
    double *sumA_arr = (samplesTOT>0 && sout) ? cf_xmalloc(sizeof(double)*N) : NULL;  /* rescaled forward sum S~[l], for sampling */
    /* regional bootstrap accumulators (the dense's regional_chunk_count_sum +
       total_regional_chunk_count): per-pop chunk count banked into a region of
       region_size cumulative count. regreg resets to 0 at each bank. */
    double *regreg = cf_xcalloc(npop, sizeof(double)); double totalreg = 0.0;

    /* emission precompute: each block writes its own disjoint Eg / CPG / CPSG region, so
       it is the one embarrassingly-parallel phase (the forward/backward recurrences below
       are sequential in both the locus and the block dimension). Threaded under -fopenmp;
       a no-op single-threaded, where the engine is meant to be scaled by running recipients
       concurrently rather than by intra-recipient threads. */
#pragma omp parallel for schedule(dynamic)
    for(int b=0;b<nb;b++){ int sb=blk[b].s,eb=blk[b].e,U=blk[b].U; int *rp=repg+(size_t)b*Umax;
        double *base=Eg+off[b]; double *cpg=CPG+(size_t)b*Umax, *cpsg=CPSG+(size_t)b*Umax;
        for(int g=0;g<U;g++){ cpg[g]=cf_cp[rp[g]]; cpsg[g]=cf_cps[rp[g]]; }
        for(int l=sb;l<eb;l++){ double *row=base+(size_t)(l-sb)*U; for(int g=0;g<U;g++) row[g]=EM(l, rp[g]); } }

    double *aprev=cf_xcalloc(K,sizeof(double));
    for(int b=0;b<nb;b++){
        int sb=blk[b].s, eb=blk[b].e, U=blk[b].U; int *gb=gidB+(size_t)b*K; int *sz=sizeg+(size_t)b*Umax;
        double *ent=AENTRY+(size_t)b*K, *Sx=SentS+(size_t)b*Umax;
        double *GFb=GF+off[b], *PFb=PF+off[b], *Egb=Eg+off[b];
        double *cpg=CPG+(size_t)b*Umax, *cpsg=CPSG+(size_t)b*Umax;
        for(int g=0;g<U;g++) Sx[g]=0.0;
        for(int i=0;i<K;i++){ double a=(b==0)?0.0:aprev[i]; ent[i]=a; if(b>0) Sx[gb[i]]+=a; }
        for(int l=sb;l<eb;l++){
            double rlm1=(l>=2)?exp(As[l-2]-As[l-1]):((l==1)?exp(-As[0]):0.0);
            double *GFr=GFb+(size_t)(l-sb)*U, *PFr=PFb+(size_t)(l-sb)*U, *Er=Egb+(size_t)(l-sb)*U;
            double sumA=0.0;
            for(int g=0;g<U;g++){
                double eg=Er[g], gv,pv;
                if(l==sb && b==0){ gv=cpsg[g]*eg; pv=0.0; }     /* copy_probSTART at locus 0 */
                else if(l==sb){ gv=eg*cpg[g]; pv=eg*(1-T[l-1])*rlm1; }
                else { double fac=eg*(1-T[l-1])*rlm1; gv=eg*cpg[g]+fac*G[g]; pv=fac*P[g]; }
                G[g]=gv; P[g]=pv; GFr[g]=gv; PFr[g]=pv; sumA += sz[g]*gv + pv*Sx[g];
            }
            if(sumA_arr) sumA_arr[l]=sumA;
            double s=(l<N-1)?sumA*T[l]:sumA; As[l]=(l>=1?As[l-1]:0.0)+log(s);
        }
        double *GFl=GFb+(size_t)(eb-1-sb)*U, *PFl=PFb+(size_t)(eb-1-sb)*U;
        for(int i=0;i<K;i++){ int g=gb[i]; aprev[i]=GFl[g]+PFl[g]*ent[i]; }
    }

    /* ---- -fold path sampling: hierarchical FFBS over the affine forward ----
       For each path, sample the donor at N-1 ~ a[N-1][.], then walk backward: stay on
       the current donor j with prob (1-T)*a[l][j] / ((1-T)*a[l][j] + T*copy_prob[j]*S~[l])
       (no recombination), else resample i ~ a[l][.] (recombination). All quantities use
       the rescaled forward; the per-locus rescale cancels in both ratios. */
    if(samplesTOT>0 && sout){
        int *locblk=cf_xmalloc(sizeof(int)*N);
        for(int b=0;b<nb;b++) for(int l=blk[b].s;l<blk[b].e;l++) locblk[l]=b;
        double *gw=cf_xmalloc(sizeof(double)*Umax);
        for(int s=0;s<samplesTOT;s++){
            int j=fold_sample_donor(N-1, Gr, GF, PF, AENTRY, SentS, off, locblk, gw);
            sout[(size_t)s*N+(N-1)]=j;
            for(int l=N-2;l>=0;l--){
                int b=locblk[l], sb=blk[b].s, U=blk[b].U, gj=gidB[(size_t)b*K+j];
                double a_lj=GF[off[b]+(size_t)(l-sb)*U+gj] + PF[off[b]+(size_t)(l-sb)*U+gj]*AENTRY[(size_t)b*K+j];
                double w_norec=(1.0-T[l])*a_lj, w_rec=T[l]*cf_cp[j]*sumA_arr[l];
                double u=(double)rand()/RAND_MAX; int i;
                if(u*(w_norec+w_rec) < w_norec) i=j;   /* no recombination: copying continues */
                else i=fold_sample_donor(l, Gr, GF, PF, AENTRY, SentS, off, locblk, gw);
                sout[(size_t)s*N+l]=i; j=i;
            }
        }
        free(locblk); free(gw);
    }

    double Asf=As[N-1];
    for(int p=0;p<npop;p++) ccpop[p]=0.0;
    if(startpop) for(int p=0;p<npop;p++) startpop[p]=0.0;
    if(out_regfinal){ for(int p=0;p<npop;p++){ out_regfinal[p]=0.0; out_regsq[p]=0.0; } }
    if(out_numreg) *out_numreg=0;
    /* EM quantities (-in N_e, -iM global mutation): N_e from per-locus etp (the
       chunkcount total) rho-weighted by the dense gd_l; per-pop expected_differences
       from e_a_l_bc = a_l*c_l*KC folded over the moments. */
    double tot_prob_Ne=0.0, tot_gd=0.0;
    for(int p=0;p<npop;p++){ ndiff[p]=0.0; nlen[p]=0.0; }
    { double S_end=exp(As[N-2]-Asf);   /* last locus N-1: c=1, a[N-1]=aprev (post-forward) */
      for(int i=0;i<K;i++){ double am=aprev[i]*S_end;   /* = e_a_l_bc at l=N-1 (copy posterior) */
        if(rh[N-1]!=donors[(size_t)(N-1)*K+i]) ndiff[pop_vec[i]]+=am;
        if(ecp_out) ecp_out[(size_t)(N-1)*npop+pop_vec[i]]+=am; } }
    double *cexit=cf_xmalloc(sizeof(double)*K); for(int i=0;i<K;i++) cexit[i]=1.0;
    { double sb0=0.0; for(int i=0;i<K;i++) sb0+=T[N-2]*cf_cp[i]*EM(N-1,i); Bs[N-1]=log(sb0); }
    for(int b=nb-1;b>=0;b--){
        int sb=blk[b].s, eb=blk[b].e, U=blk[b].U; int *gb=gidB+(size_t)b*K;
        int *sz=sizeg+(size_t)b*Umax;
        int top=(eb-1<N-2)?eb-1:N-2;
        double *w=WB+(size_t)b*K, *Sw=SwS, *ent=AENTRY+(size_t)b*K;
        double *GFb=GF+off[b], *PFb=PF+off[b], *Egb=Eg+off[b]; double *cpg=CPG+(size_t)b*Umax;
        for(int g=0;g<U;g++) Sw[g]=0.0;
        for(int gp=0;gp<U*npop;gp++){ szP[gp]=0.0; SaP[gp]=0.0; SwP[gp]=0.0; SawP[gp]=0.0; }
        for(int i=0;i<K;i++){ double wi=EM(top+1,i)*cexit[i]; w[i]=wi; int g=gb[i], p=pop_vec[i]; double ei=ent[i];
            Sw[g]+=wi; int idx=g*npop+p; szP[idx]+=1.0; SaP[idx]+=ei; SwP[idx]+=wi; SawP[idx]+=ei*wi; }
        for(int l=top;l>=sb;l--){
            double rb=(l+2<=N-1)?exp(Bs[l+2]-Bs[l+1]):exp(-Bs[N-1]);
            double *Er1=(l+1<eb)?Egb+(size_t)(l+1-sb)*U:0;
            if(l==top){ for(int g=0;g<U;g++){ GBc[g]=1.0; PBc[g]=(1-T[l])*rb; } }
            else { for(int g=0;g<U;g++){ double fac=(1-T[l])*Er1[g]*rb; GBc[g]=1.0+fac*GBp[g]; PBc[g]=fac*PBp[g]; } }
            if(l>0){ double *Er=Egb+(size_t)(l-sb)*U; double sB=0.0;
                for(int g=0;g<U;g++) sB+=cpg[g]*Er[g]*(sz[g]*GBc[g]+PBc[g]*Sw[g]); Bs[l]=Bs[l+1]+log(T[l-1]*sB); }
            double BsR=(l+2<=N-1)?Bs[l+2]:0.0, Asm1=(l>=1)?As[l-1]:0.0;
            double KF1=exp(As[l]+BsR-Asf), KF=exp(Asm1+BsR-Asf), om=(1-T[l]);
            double *GFr=GFb+(size_t)(l-sb)*U, *PFr=PFb+(size_t)(l-sb)*U;
            double etpl=0.0;     /* per-locus chunkcount total = expected_transition_prob[l] (for N_e) */
            /* expected_chunk_length integrand = G_l*0.5*(e_a_lp1_bp + e_a_l_bc) (the
               tp_from_i_to_i terms cancel); G_l = 100*(pos[l+1]-pos[l])*delta*lambda[l] (cM). */
            double Glh = (lam_g[l]>=0) ? 0.5*100.0*(pos_g[l+1]-pos_g[l])*delta_g*lam_g[l] : 0.0;
            if(l < eb-1){
                double *GFr1=GFb+(size_t)(l+1-sb)*U, *PFr1=PFb+(size_t)(l+1-sb)*U;
                for(int g=0;g<U;g++){
                    double eg=Er1[g];
                    double XG=GFr1[g]*KF1 - GFr[g]*KF*eg*om, XP=PFr1[g]*KF1 - PFr[g]*KF*eg*om;
                    double GBl1,PBl1; if(l+1==N-1){GBl1=1.0;PBl1=0.0;} else {GBl1=GBp[g];PBl1=PBp[g];}
                    double cA=GBl1*XG, cB=GBl1*XP, cC=PBl1*XG, cD=PBl1*XP; int base=g*npop;
                    double aA=KF1*GFr1[g]*GBl1, aB=KF1*GFr1[g]*PBl1, aC=KF1*PFr1[g]*GBl1, aD=KF1*PFr1[g]*PBl1; /* e_a_lp1_bp */
                    for(int p=0;p<npop;p++){ double s_=szP[base+p]; if(s_==0.0) continue;
                        double inc=cA*s_ + cB*SaP[base+p] + cC*SwP[base+p] + cD*SawP[base+p];
                        ccpop[p]+=inc; etpl+=inc; regreg[p]+=inc;
                        nlen[p]+=Glh*(aA*s_ + aB*SwP[base+p] + aC*SaP[base+p] + aD*SawP[base+p]); }
                }
            } else if(b+1<nb){
                int co=Gr->celloff[b], nc=Gr->celloff[b+1]-co; int *cidB=Gr->cellB+(size_t)b*K;
                double *cx2=CEXIT+(size_t)(b+1)*K;
                double *GF2=GF+off[b+1], *PF2=PF+off[b+1], *Eg2=Eg+off[b+1];
                for(int x=0;x<nc*npop;x++){ Mc[x]=0.0; Mce[x]=0.0; }
                for(int i=0;i<K;i++){ int c=cidB[i], p=pop_vec[i]; double cv=cx2[i];
                    Mc[c*npop+p]+=cv; Mce[c*npop+p]+=cv*ent[i]; }
                int *cGb=Gr->cellGb+co, *cG2=Gr->cellG2+co;
                for(int c=0;c<nc;c++){ int g=cGb[c], g2=cG2[c];
                    double C0=KF1*GF2[g2], C1=KF1*PF2[g2]-KF*om*Eg2[g2];
                    double cf0=C0+C1*GFr[g], cf1=C1*PFr[g]; int base=c*npop;
                    /* e_a_lp1_bp at boundary: a_{eb}=GF2+PF2*a_l, a_l=GFr+PFr*ent, so
                       a_{eb}*cx2 summed = (GF2+PF2*GFr)*Mc + PF2*PFr*Mce, times KF1. */
                    double bA=KF1*(GF2[g2]+PF2[g2]*GFr[g]), bB=KF1*PF2[g2]*PFr[g];
                    for(int p=0;p<npop;p++){ double inc=cf0*Mc[base+p]+cf1*Mce[base+p]; ccpop[p]+=inc; etpl+=inc; regreg[p]+=inc;
                        nlen[p]+=Glh*(bA*Mc[base+p]+bB*Mce[base+p]); } }
            }
            /* N_e: dense gd_l rho-weight (pos/delta/lambda). gd>0 always (duplicate SNP
               positions are rejected/jittered at read time), so the gd>0 test is defensive,
               skipping only the gd==0 limit (0/0) the dense never reaches either. */
            if(lam_g[l]>=0){ double gd=(pos_g[l+1]-pos_g[l])*delta_g*lam_g[l];
                if(gd>0.0){ tot_gd+=gd; tot_prob_Ne+=(rho_g*gd/(1.0-exp(-rho_g*gd)))*etpl; } }
            if(etp_out) etp_out[l]=etpl;   /* per-locus transition prob (for -d), l in 0..N-2 */
            /* regional bootstrap bank: replicate the dense's per-locus region cut
               (ChromoPainterSampler.c backwardAlgorithm, the region-bank loop). etpl is this
               locus's total chunk count (= dense total_prob); the walk is the same
               reverse N-2..0 order, so the region partition matches the dense up to
               the ~1e-9 FP gap in the per-locus total. rounding_val=1e-7. The final
               partial region (totalreg < region_size at the end) is NOT banked, as in
               the dense. */
            if(out_regfinal){ totalreg += etpl;
              if(totalreg + 1e-7 >= region_size){
                for(int p=0;p<npop;p++){ out_regfinal[p]+=regreg[p]; out_regsq[p]+=regreg[p]*regreg[p]; regreg[p]=0.0; }
                totalreg=0.0; if(out_numreg)(*out_numreg)++;
              } }
            /* per-pop expected_differences (e_a_l_bc=a_l*c_l*KC, mismatched donors) +
               the e_a_l_bc half of expected_chunk_length (all donors). CURRENT GBc/PBc.
               Summed over groups per pop, e_a_l_bc is the -b copy posterior at this locus. */
            { double KC=KF/rb;   /* == exp(Asm1+Bs[l+1]-Asf); one fewer exp per locus */
              int *rp=repg+(size_t)b*Umax;
              for(int g=0;g<U;g++){
                double gA=GFr[g],pA=PFr[g],gB=GBc[g],pB=PBc[g]; int base=g*npop;
                int mm = (rh[l]!=donors[(size_t)l*K+rp[g]]);
                for(int p=0;p<npop;p++){ double sz_=szP[base+p]; if(sz_==0.0) continue;
                    double albc=KC*(gA*gB*sz_ + gA*pB*SwP[base+p] + pA*gB*SaP[base+p] + pA*pB*SawP[base+p]);
                    nlen[p]+=Glh*albc; if(mm) ndiff[p]+=albc;
                    if(ecp_out) ecp_out[(size_t)l*npop+p]+=albc; } } }
            { double *t; t=GBp;GBp=GBc;GBc=t; t=PBp;PBp=PBc;PBc=t; }
        }
        double *cxb=CEXIT+(size_t)b*K;
        for(int i=0;i<K;i++){ int g=gb[i]; double cv=GBp[g]+PBp[g]*w[i]; cexit[i]=cv; cxb[i]=cv; }
    }
    { double k0=exp(Bs[1]-Asf); int *gb0=gidB; double *ent0=AENTRY, *cx0=CEXIT; double *GF0=GF+off[0], *PF0=PF+off[0];
      /* start term (= dense copy_prob_newSTART, the locus-0 posterior): folded into
         ccpop to give the full corrected chunk count, and reported separately in
         startpop so the caller can split off the no-start copy_prob_new (-ip). */
      for(int i=0;i<K;i++){ int g=gb0[i]; double a0=GF0[g]+PF0[g]*ent0[i]; double st=a0*cx0[i]*k0;
        ccpop[pop_vec[i]] += st; if(startpop) startpop[pop_vec[i]] += st; } }

    *Ne_out = (tot_gd>0.0) ? tot_prob_Ne/tot_gd : 0.0;   /* N_e EM estimate (-in); dense floor */
    if(*Ne_out < MIN_NE) *Ne_out = MIN_NE;
    if(loglik_out) *loglik_out = Asf;   /* forward log-likelihood (= dense Alphasum) */

    free(off);free(GF);free(PF);free(Eg);free(AENTRY);free(WB);free(CEXIT);free(As);free(Bs);
    free(SentS);free(SwS);free(G);free(P);free(GBp);free(PBp);free(GBc);free(PBc);
    free(szP);free(SaP);free(SwP);free(SawP);free(aprev);free(cexit);free(Mc);free(Mce);free(CPG);free(CPSG);
    free(sumA_arr);free(regreg);
    #undef EM
}

/* Entry point. Builds the locus-major donor buffer from existing_h (hap-major), the
   recipient row from newh, the param-class seed and the (panel-fixed) grouping, then
   runs the fold. out_ccpop[ndonorpops]. */
void cpfold_perpop(signed char *newh, signed char **existing_h, int nhaps, int nloci,
                   double *TransProb, double *MutProb_vec, double *copy_prob, double *copy_probSTART,
                   double *pos, double *lambda, double delta, double rhobar,
                   int *pop_vec_in, int ndonorpops, int Ustar,
                   double *out_ccpop, double *out_start, double *out_ndiff, double *out_nlen, double *out_Ne,
                   double *out_loglik, double *out_etp, double *out_ecp,
                   double region_size, double *out_regfinal, double *out_regsq, int *out_numreg,
                   int samplesTOT, int *out_samples){
    K=nhaps; N=nloci; npop=ndonorpops;
    cf_cp=copy_prob; cf_cps=copy_probSTART; cf_mut=MutProb_vec;
    pop_vec=pop_vec_in; T=TransProb;
    pos_g=pos; lam_g=lambda; delta_g=delta; rho_g=rhobar;
    /* dense id of the distinct (cp,cps,mut) tuples; the grouping seeds on it so every
       group is parameter-homogeneous. nclass=1 (uniform) is the fast common path. */
    cf_pclass=cf_xmalloc((size_t)K*sizeof(int)); cf_nclass=0;
    { int *rep=cf_xmalloc((size_t)K*sizeof(int));
      for(int i=0;i<K;i++){ int c=-1;
        for(int j=0;j<cf_nclass;j++){ int r=rep[j];
          if(copy_prob[i]==copy_prob[r]&&copy_probSTART[i]==copy_probSTART[r]&&MutProb_vec[i]==MutProb_vec[r]){ c=j; break; } }
        if(c<0){ c=cf_nclass; rep[cf_nclass++]=i; } cf_pclass[i]=c; }
      free(rep); }
    /* locus-major donor buffer from existing_h (hap-major), rebuilt every call. */
    donors=cf_xmalloc((size_t)N*K);
    /* hap-major existing_h -> locus-major donors. Cache-blocked: parallelize over
       disjoint locus tiles (no false sharing on the shared donors buffer) and tile the
       donor dimension so writes to donors[l*K + ib..ie] land contiguously (one cache
       line filled at a time) instead of strided by K. The naive strided-write transpose
       was ~half the chr-scale fold wall-clock. Output is byte-identical. */
    /* Tile sizes: TI=64 is one cache line, so the strided donor writes fill whole
       lines; TL=512 keeps the read block (TI*TL = 32 KiB) inside L1. Wall-clock is
       flat across TL in 256..4096 on tested hardware (a modern L2 holds the block at
       any of these), so this is a robustness choice sized to the smallest common L1,
       not a machine-tuned constant. */
    { const int TL=512, TI=64;
#pragma omp parallel for schedule(static)
      for(int lb=0; lb<N; lb+=TL){ int le=(lb+TL<N)?lb+TL:N;
        for(int ib=0; ib<K; ib+=TI){ int ie=(ib+TI<K)?ib+TI:K;
          for(int l=lb; l<le; l++){ uint8_t *out=donors+(size_t)l*K;
            for(int i=ib; i<ie; i++) out[i]=(uint8_t)existing_h[i][l]; } } } }
    uint8_t *rh=cf_xmalloc(N); for(int l=0;l<N;l++) rh[l]=(uint8_t)newh[l];

    Groups Gr=build_groups_adaptive(Ustar);
    fold_cc(rh, out_ccpop, out_start, out_ndiff, out_nlen, out_Ne, out_loglik, out_etp, out_ecp,
            region_size, out_regfinal, out_regsq, out_numreg, samplesTOT, out_samples, &Gr);

    free_groups(&Gr); free(donors); donors=NULL; free(cf_pclass); free(rh);
}

/* Defensive teardown of the engine's static donor buffer. cpfold_perpop frees it per
   call, so this is normally a no-op; idempotent and safe to call at teardown. */
void cpfold_cleanup(void){ free(donors); donors=NULL; }
