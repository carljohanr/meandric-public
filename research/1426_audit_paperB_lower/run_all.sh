#!/usr/bin/env bash
# 1426 (JB): full reproduction of the audit (about 3 minutes, peak 2.5 GB, 4 threads). Run from this folder.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs work
E=../1421_export_lower_certificate/export
PY=${PYTHON:-python3}
VPY=$PY  # mpmath is listed in requirements.txt
if [ ! -f "$E/parent_r6.bin" ]; then
  echo "This audit needs a locally generated full certificate export; see docs/REPRODUCING.md." >&2
  exit 2
fi
# Original source-certificate hashes are recorded in SPEC.md. Locally generated
# candidates can differ and must pass the mathematical checks below.
$PY law.py $E/law.json > logs/law.log                       # C1
$PY prep.py                                                 # primes (kept if work/primes.json exists), law tables
$PY parent.py $E/parent_r6.bin > logs/parent.log            # C2, C3 exact; parent densities for C
$PY walks.py $E > logs/walks.log                            # C4; edge terms; closed-form mass
clang -O3 -march=native -w -o work/flowcheck flowcheck.c -lpthread
clang -O2 -w -o work/primtest primtest.c -lpthread
$VPY primtest.py > logs/primtest.log                        # unit test of log, beta mantissas, beta mod p
./work/flowcheck 6 -1 4 1 > logs/flow_r6_parent.log
./work/flowcheck 8 0 4 1 > logs/flow_r8_F0.log
./work/flowcheck 9 1 4 1 > logs/flow_r9_F1.log
./work/flowcheck 10 2 4 1 > logs/flow_r10_F2.log
./work/flowcheck 10 3 4 1 > logs/flow_r10_F3.log
./work/flowcheck 10 4 4 1 > logs/flow_r10_F4.log
$PY bound.py logs/flow_r8_F0.log 0 6106888878721/500000000000 > logs/bound_F0.log
$PY bound.py logs/flow_r9_F1.log 1 3053649017533/250000000000 > logs/bound_F1.log
$PY bound.py logs/flow_r10_F2.log 2 12216487956539/1000000000000 > logs/bound_F2.log
$PY bound.py logs/flow_r10_F3.log 3 3054585940667/250000000000 > logs/bound_F3.log
$PY bound.py logs/flow_r10_F4.log 4 6109550428547/500000000000 > logs/bound_F4.log
# second, independent pair of random primes: JB_WORK=$PWD/work2, rerun prep/parent/walks, then flowcheck with work2 files.
# damage tests: see logs/damage.log for the exact commands (walks.py --flip/--flip-last/--force-include/--scale-coef,
# parent.py --perturb, JB_ZETA=1/7, a damaged copy of law.json).
