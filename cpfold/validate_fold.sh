#!/usr/bin/env bash
# validate_fold.sh - exactness gate for the exact block fold (-fold) vs the
# dense ChromoPainter forward-backward already in the tree.
#
# For each (dataset, mode) it runs two ChromoPainter configurations and asserts
# the standard output files agree to printed precision (relative tolerance 1e-6):
#   1. dense  (default forward-backward)      -- the reference
#   2. -fold  (the exact O(N*Umean) block fold)
# plus an independent from-scratch oracle (oracle/independent_oracle.py) that
# shares no code with finestructure4, to catch a bug common to dense and fold.
#
# The committed oracle/ case (oracle/*) is self-contained and always runs. The win/
# fold-vs-dense matrix needs developer-generated panels in DATADIR (win*.phase/.recom);
# when absent, those sections are skipped and the script reports how many checks ran -
# it never prints success having validated nothing.
#
# Usage:  ./validate_fold.sh [FS_BINARY] [DATADIR]
#   FS_BINARY defaults to ../fs ; DATADIR defaults to ./win
set -u
FS="${1:-../fs}"
D="${2:-./win}"
ID="$D/cp.idfile"; POP="$D/cp.poplist"
NE=400000; MUT=0.0006338578
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
fail=0; checks=0   # checks = comparisons actually run; 0 => nothing validated (a failure)

# Standard per-pop output files the fold reproduces. The dense painter computes the
# forward-backward in FLOAT (the linear-space rescale: float Alphamat + double logScale)
# while the fold accumulates its chunk counts in DOUBLE, so the two agree to printed
# precision but NOT byte-for-byte: on these data the chunk counts differ only in the last
# printed digit (max rel ~1e-9). The comparison is therefore relative-tolerance (cmp_num,
# 1e-6), not byte-identity; the fold's ABSOLUTE exactness is anchored by the independent
# double-precision oracle below (fold == oracle to <1e-5). EMprobs (cmp_emprobs), the
# regional bootstrap (cmp_region) and the -b/-d per-locus outputs (gzcmp) check at 1e-6 too.
EXTS="chunkcounts chunklengths mutationprobs prop"

run() { # env_prefix phase recom extra_args outprefix
  env $1 OMP_NUM_THREADS=1 "$FS" cp -g "$2" -r "$3" -t "$ID" -f "$POP" 0 0 \
      -s 0 $4 -n "$NE" -M "$MUT" -o "$5" >/dev/null 2>&1
}

# cmp_num: numeric-file comparator - same column count per row, every numeric field
# within rel tol 1e-6 (dense float vs fold double => last-digit differences), non-numeric
# fields (recipient IDs) exact, NaN/Inf rejected, same line count. Returns 0 = match.
cmp_num() { # ref cand
  awk '
    function abs(x){return x<0?-x:x}
    FNR==NR{ nf[FNR]=NF; for(i=1;i<=NF;i++) v[FNR,i]=$i; refn=FNR; next }
    { candn=FNR;
      if(NF!=nf[FNR]) bad=1;
      for(i=1;i<=NF;i++){
        if($i ~ /[nN][aA][nN]|[iI][nN][fF]/ || v[FNR,i] ~ /[nN][aA][nN]|[iI][nN][fF]/) bad=1;
        else if($i ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && v[FNR,i] ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/){
          d=abs(($i+0)-(v[FNR,i]+0)); rel=d/(abs(v[FNR,i]+0)+1e-300);
          # float dense vs double fold: equal if within 3e-6 absolute OR 1e-6 relative.
          # (sub-1 posteriors printed at 6 decimals can flip 1-2 units in the last place;
          # large chunk counts differ ~1e-9 relative. Fails only if BOTH bounds exceeded.)
          if(d>3e-6 && rel>1e-6) bad=1;
        } else if(v[FNR,i]!=$i) bad=1; } }
    END{ if(refn!=candn) bad=1; exit bad?1:0 }' "$1" "$2"
}
cmp_set() { # ref prefix label
  local ok=1 n=0
  for e in $EXTS; do
    [ -f "$1.$e.out" ] || continue          # dense did not produce this ext: not applicable
    n=$((n+1))
    if [ ! -f "$2.$e.out" ]; then ok=0; echo "    MISSING (fold): $e"; continue; fi
    if ! cmp_num "$1.$e.out" "$2.$e.out"; then ok=0; echo "    DIFFER (>1e-6): $e"; fi
  done
  checks=$((checks+1))
  if [ $n -eq 0 ]; then echo "  FAIL  $3 (no dense output - run failed?)"; fail=1; return; fi
  if [ $ok -eq 1 ]; then echo "  PASS  $3"; else echo "  FAIL  $3"; fail=1; fi
}

# EMprobs comparator: enforce the SCHEMA (same column count per row - the fold
# regression was emitting 3 cols instead of 5) and that every field matches the
# reference, with a relative tolerance of 1e-6 on numeric fields (the N_e column
# carries FP-reorder noise in its last %.10lf digit across dense/fold).
cmp_emprobs() { # ref prefix label
  local r="$1.EMprobs.out" f="$2.EMprobs.out"; checks=$((checks+1))
  [ -f "$r" ] && [ -f "$f" ] || { echo "  FAIL  $3 (EMprobs missing)"; fail=1; return; }
  if awk '
      function abs(x){return x<0?-x:x}
      FNR==NR{ nf[FNR]=NF; for(i=1;i<=NF;i++) v[FNR,i]=$i; refn=FNR; next }
      { candn=FNR;
        if(NF!=nf[FNR]){ bad=1 }                                  # schema: same column count
        for(i=1;i<=NF;i++){
          if($i ~ /[nN][aA][nN]|[iI][nN][fF]/ || v[FNR,i] ~ /[nN][aA][nN]|[iI][nN][fF]/){ bad=1 }  # reject NaN/Inf
          else if($i ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && v[FNR,i] ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/){
            d=abs(($i+0)-(v[FNR,i]+0)); rel=(abs(v[FNR,i]+0)>1?d/abs(v[FNR,i]+0):d);
            if(rel>1e-6) bad=1;
          } else if(v[FNR,i]!=$i){ bad=1 } } }
      END{ if(refn!=candn) bad=1; exit bad?1:0 }' "$r" "$f"; then   # same line count (no dropped rows)
    echo "  PASS  $3 (EMprobs schema + values, N_e tol 1e-6)"
  else echo "  FAIL  $3 (EMprobs schema/value/linecount mismatch)"; fail=1; fi
}

for ds in win200 win800 win win20k; do
  ph="$D/$ds.phase"; rc="$D/$ds.recom"
  [ -f "$ph" ] && [ -f "$rc" ] || continue
  echo "=== dataset $ds ($(sed -n 2p "$ph") SNPs) ==="
  # The fold drives the full per-pop E-M loop exactly: -ip (copy proportions, from
  # the per-pop posterior chunk count + start term) and -im (per-pop mutation, from
  # the per-pop expected differences), alongside -in (N_e) and -iM (global mutation),
  # so .prop and every other file tracks the dense exactly. (-im uses -i 6, not -i 10:
  # on the tiny down-sampled win20k the per-pop mutation over-iterates past a valid
  # rate by -i 10 in BOTH backends identically - a property of the data, not the
  # fold - which the tolerance check still passes but is not a meaningful config.)
  for mode in "-i 0" "-i 6 -in -iM" "-i 10 -ip" "-i 10 -ip -in -iM" "-i 6 -im" "-i 6 -ip -im -in"; do
    run ""         "$ph" "$rc" "$mode"       "$TMP/dense"
    run ""         "$ph" "$rc" "$mode -fold" "$TMP/fold"
    cmp_set     "$TMP/dense" "$TMP/fold" "fold == dense   [$mode]"
    cmp_emprobs "$TMP/dense" "$TMP/fold" "fold == dense   [$mode]"
  done
done

# The win-based sections below need the generated win panel; gate on it so a fresh
# clone (oracle/ committed, win/ absent) validates via the oracle alone, without
# running fs on a missing panel.
if [ -f "$D/win.phase" ]; then
# --- per-pop fold exactness (PR2b): FIXED per-pop copy probs (-p) and/or FIXED
#     per-pop mutation rates (-m). The fold groups donors by (substring, pop),
#     so dense == fold at the per-pop output granularity. -m excludes -M. ---
runpp() { # poplist extra_args outprefix
  env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" \
      -t "$ID" -f "$1" 0 0 -s 0 -i 0 -n "$NE" $2 -o "$3" >/dev/null 2>&1
}
PR="$D/cp.poplist.prior"; PM="$D/cp.poplist.mut"; PB="$D/cp.poplist.pm"
echo "=== per-pop fold == dense (win, -i 0) ==="
if [ -f "$PR" ]; then
  runpp "$PR" "-M $MUT -p"        "$TMP/d_p";  runpp "$PR" "-M $MUT -p -fold"        "$TMP/f_p"
  cmp_set "$TMP/d_p"  "$TMP/f_p"  "fold == dense   [-p]"
fi
if [ -f "$PM" ]; then
  runpp "$PM" "-m $MUT"           "$TMP/d_m";  runpp "$PM" "-m $MUT -fold"           "$TMP/f_m"
  cmp_set "$TMP/d_m"  "$TMP/f_m"  "fold == dense   [-m]"
fi
if [ -f "$PB" ]; then
  runpp "$PB" "-p -m $MUT"        "$TMP/d_pm"; runpp "$PB" "-p -m $MUT -fold"        "$TMP/f_pm"
  cmp_set "$TMP/d_pm" "$TMP/f_pm" "fold == dense   [-p -m]"
fi

# --- all-vs-all (-a): -fold is exact here too (each individual is its own donor
#     pop, so the within-pop redistribution is a no-op). ---
echo "=== all-vs-all fold == dense (win, -a 0 0 -i 0) ==="
env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" -a 0 0 \
    -s 0 -i 0 -n "$NE" -M "$MUT"       -o "$TMP/d_a" >/dev/null 2>&1
env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" -a 0 0 \
    -s 0 -i 0 -n "$NE" -M "$MUT" -fold -o "$TMP/f_a" >/dev/null 2>&1
cmp_set "$TMP/d_a" "$TMP/f_a" "fold == dense   [-a 0 0]"

# --- regional bootstrap: -fold now PRODUCES the regional files (it replicates the
#     dense's per-locus region banking). They must be PRESENT and match the dense.
#     On this tiny win data the region boundaries do not straddle the FP gap, so the
#     files are byte-identical; cmp_region uses a 1e-6 numeric tolerance anyway, since
#     at chr scale a boundary can shift one locus (the real gate there is chromocombine
#     `c`, robust to such shifts - not exercised on this tiny data). ---
echo "=== -fold regional bootstrap == dense + per-locus outputs ==="
# cmp_region: header identical; per row same recipient + num.regions; the per-pop
# region totals (and sums-of-squares) within rel 1e-6. Same line/column count.
cmp_region() { # ref prefix label
  local ok=1; checks=$((checks+1))
  for e in regionchunkcounts regionsquaredchunkcounts; do
    if [ ! -f "$2.$e.out" ]; then echo "    MISSING (fold): $e"; ok=0; continue; fi
    if [ ! -f "$1.$e.out" ]; then echo "    MISSING (dense): $e"; ok=0; continue; fi
    awk '
      function abs(x){return x<0?-x:x}
      FNR==NR{ nf[FNR]=NF; for(i=1;i<=NF;i++) v[FNR,i]=$i; refn=FNR; next }
      { candn=FNR;
        if(NF!=nf[FNR]) bad=1;                                   # same column count
        for(i=1;i<=NF;i++){
          if($i ~ /[nN][aA][nN]|[iI][nN][fF]/){ bad=1 }
          else if($i ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && v[FNR,i] ~ /^-?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/){
            d=abs(($i+0)-(v[FNR,i]+0)); rel=(abs(v[FNR,i]+0)>1?d/abs(v[FNR,i]+0):d);
            if(rel>1e-6) bad=1;
          } else if(v[FNR,i]!=$i){ bad=1 } } }
      END{ if(refn!=candn) bad=1; exit bad?1:0 }' "$1.$e.out" "$2.$e.out" \
      || { echo "    DIFFER: $e"; ok=0; }
  done
  if [ $ok -eq 1 ]; then echo "  PASS  $3"; else echo "  FAIL  $3"; fail=1; fi
}
env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" 0 0 \
    -s 0 -i 0 -k 5 -n "$NE" -M "$MUT"       -o "$TMP/reg_d" >/dev/null 2>&1
env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" 0 0 \
    -s 0 -i 0 -k 5 -n "$NE" -M "$MUT" -fold -o "$TMP/reg_f" >/dev/null 2>&1
cmp_region "$TMP/reg_d" "$TMP/reg_f" "fold regional == dense   [-k 5]"
# -b (.copyprobsperlocus, per-locus per-pop copy posterior) and -d (.transitionprobs,
# per-locus transition prob) ARE produced by the fold, matching the dense to 1e-6.
gzcmp() { # ext label  (compares <prefix>.ext.gz under $TMP/bd_d vs $TMP/bd_f, rel 1e-6)
  checks=$((checks+1))
  if [ ! -s "$TMP/bd_d.$1" ] || [ ! -s "$TMP/bd_f.$1" ]; then echo "  FAIL  $2 (missing .gz: $1)"; fail=1; return; fi
  gunzip -c "$TMP/bd_d.$1" > "$TMP/gz_d" 2>/dev/null; gunzip -c "$TMP/bd_f.$1" > "$TMP/gz_f" 2>/dev/null
  if [ ! -s "$TMP/gz_d" ] || [ ! -s "$TMP/gz_f" ]; then echo "  FAIL  $2 (empty after gunzip: $1)"; fail=1; return; fi
  if cmp_num "$TMP/gz_d" "$TMP/gz_f"; then echo "  PASS  $2"; else echo "  FAIL  $2"; fail=1; fi
}
for mode in "-i 0" "-i 6 -in -iM"; do
  env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" 0 0 \
      -s 0 $mode -n "$NE" -M "$MUT" -b -d       -o "$TMP/bd_d" >/dev/null 2>&1
  env OMP_NUM_THREADS=1 "$FS" cp -g "$D/win.phase" -r "$D/win.recom" -t "$ID" -f "$POP" 0 0 \
      -s 0 $mode -n "$NE" -M "$MUT" -b -d -fold -o "$TMP/bd_f" >/dev/null 2>&1
  gzcmp copyprobsperlocus.out.gz "fold -b == dense   [$mode]"
  gzcmp transitionprobs.out.gz   "fold -d == dense   [$mode]"
done
fi   # end win-panel-dependent sections

# --- committed example data (examples/example1): real fold == dense on the maintainer's
#     own ChromoPainter example (all-vs-all on a few recipients), so a fresh clone / CI
#     exercises the painter on real data with NO extra data shipped. ---
EX="$(dirname "$0")/../examples/example1"
if [ -f "$EX/example_cp.phase" ]; then
  echo "=== examples/example1: real fold == dense (all-vs-all, -a 1 4) ==="
  for mode in "-i 0" "-i 2 -in -iM"; do
    env OMP_NUM_THREADS=1 "$FS" cp -g "$EX/example_cp.phase" -r "$EX/example_cp.recombfile" \
        -t "$EX/example_cp.ids" -a 1 4 -s 0 $mode       -o "$TMP/ex_d" >/dev/null 2>&1
    env OMP_NUM_THREADS=1 "$FS" cp -g "$EX/example_cp.phase" -r "$EX/example_cp.recombfile" \
        -t "$EX/example_cp.ids" -a 1 4 -s 0 $mode -fold -o "$TMP/ex_f" >/dev/null 2>&1
    cmp_set     "$TMP/ex_d" "$TMP/ex_f" "fold == dense   [example1 -a 1 4, $mode]"
    cmp_emprobs "$TMP/ex_d" "$TMP/ex_f" "fold == dense   [example1 -a 1 4, $mode]"
  done
fi

# --- independent oracle: fs (dense AND fold) vs a from-scratch reference that
#     shares NO code with finestructure4. The reference (oracle/independent_oracle.py)
#     computes the per-pop chunk counts by exact brute-force enumeration of all 4^6
#     donor paths, cross-checked by a formula-free Monte Carlo. This catches a shared
#     bug in the emission / transition / chunk-count code that a fold-vs-dense diff
#     cannot (both share cp_emis). The case includes a missing allele (9) and a gap allele (8). ---
OD="$(dirname "$D")/oracle"
if [ -f "$OD/data.phase" ]; then
  echo "=== independent oracle (brute-force + MC) vs fs dense and fold ==="
  EXP_A=1.4538194056; EXP_B=0.2377183618   # from oracle/expected.txt (regenerate: python3 oracle/independent_oracle.py)
  orun() { # extra outprefix
    env OMP_NUM_THREADS=1 $1 "$FS" cp -g "$OD/data.phase" -r "$OD/data.recom" \
        -t "$OD/id.txt" -f "$OD/poplist.txt" 0 0 -j -s 0 -i 0 -n 100 -M 0.01 $2 -o "$3" >/dev/null 2>&1
  }
  ocheck() { # chunkcounts-file label  (self-contained: sets fail=1 on any miss)
    checks=$((checks+1))
    if [ ! -f "$1" ]; then echo "  FAIL  $2 (no output)"; fail=1; return; fi
    awk -v ea="$EXP_A" -v eb="$EXP_B" -v lab="$2" '
      function abs(x){return x<0?-x:x}
      END{
        if(a=="" ){ print "  FAIL  "lab" (no output)"; exit 1 }
        da=abs(a-ea); db=abs(b-eb);
        if(da<1e-5 && db<1e-5) printf "  PASS  %s (popA %.6f popB %.6f, |d|<1e-5 vs exact)\n",lab,a,b;
        else { printf "  FAIL  %s popA=%.6f(d=%.1e) popB=%.6f(d=%.1e)\n",lab,a,da,b,db; exit 1 }
      }
      $1=="TGT"{a=$2;b=$3}' "$1" || fail=1
  }
  orun ""        ""      "$TMP/odense"; ocheck "$TMP/odense.chunkcounts.out" "dense == oracle"
  orun ""        "-fold" "$TMP/ofold";  ocheck "$TMP/ofold.chunkcounts.out"  "fold  == oracle"

  # -fold PATH SAMPLING (-s): the hierarchical sampler draws copying paths from the
  # affine forward (no Alphamat). Samples are distributionally - not byte - identical
  # to the dense, so validate statistically: the empirical per-locus popA copy
  # frequency from many samples must converge to the exact analytic posterior (the
  # -fold -b output). Tested single-block (-fold) and multi-block (-foldU 2). popA
  # donor output values for this case are {2,3}; loci at pos 100..1600.
  echo "=== -fold path sampling draws from the correct posterior (oracle) ==="
  orun "" "-fold -b" "$TMP/osb"
  gunzip -c "$TMP/osb.copyprobsperlocus.out.gz" > "$TMP/oanl.txt" 2>/dev/null
  sampcheck() { # foldargs label  (self-contained: sets fail=1; guards both inputs)
    checks=$((checks+1))
    if [ ! -s "$TMP/oanl.txt" ]; then echo "  FAIL  $2 (analytic -fold -b posterior missing)"; fail=1; return; fi
    env OMP_NUM_THREADS=1 "$FS" cp -g "$OD/data.phase" -r "$OD/data.recom" -t "$OD/id.txt" \
        -f "$OD/poplist.txt" 0 0 -j -s 100000 -S 9 -i 0 -n 100 -M 0.01 $1 -o "$TMP/osamp" >/dev/null 2>&1
    gunzip -c "$TMP/osamp.samples.out.gz" > "$TMP/osamp.txt" 2>/dev/null
    if [ ! -s "$TMP/osamp.txt" ]; then echo "  FAIL  $2 (no samples)"; fail=1; return; fi
    awk -v lab="$2" '
      function abs(x){return x<0?-x:x}
      FNR==NR{ if($1+0>0 && $0 !~ /[A-Za-z]/) anl[$1]=$2; next }   # analytic: pos -> popA posterior
      /^[0-9]+ / && NF==7 { n++; pp[1]=100;pp[2]=200;pp[3]=400;pp[4]=700;pp[5]=1100;pp[6]=1600;
        for(l=1;l<=6;l++){ d=$(l+1); if(d==2||d==3) cnt[l]++ } }
      END{ if(n==0){print "  FAIL  "lab" (no samples)"; exit 1}
        mx=0; for(l=1;l<=6;l++){ e=cnt[l]/n; dd=abs(e-anl[pp[l]]); if(dd>mx)mx=dd }
        if(mx<0.03) printf "  PASS  %s (max|emp-analytic popA|=%.4f, n=%d)\n",lab,mx,n;
        else { printf "  FAIL  %s (max diff %.4f)\n",lab,mx; exit 1 } }' "$TMP/oanl.txt" "$TMP/osamp.txt" || fail=1
  }
  sampcheck "-fold"    "fold sampling ~ posterior   [1 block]"
  sampcheck "-foldU 2" "fold sampling ~ posterior   [multi-block]"
fi

echo
if [ "$checks" -eq 0 ]; then
  echo "NO CHECKS RAN - no test data found (win/ panels are developer-generated; the committed"
  echo "oracle/ case should always run - is FS_BINARY=$FS valid?). Nothing was validated."
  exit 1
fi
if [ $fail -eq 0 ]; then echo "ALL EXACT ($checks checks; fold == dense to 1e-6 at printed precision; fold == oracle exactly)"; else echo "SOME CONFIGS DIFFER - see FAIL lines above"; fi
exit $fail
