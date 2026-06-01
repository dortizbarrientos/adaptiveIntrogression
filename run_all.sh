#!/usr/bin/env bash
# =============================================================================
#  run_all.sh  --  reproduce the study on the REAL Populus data.
#
#  This is the entry point for the paper's results. It runs every pipeline layer
#  and reads Yupeng's real data in 04_futureClimate/data_real/. For the synthetic
#  test harness (how we verify the code), see run_mock.sh instead.
#
#  USAGE
#    bash run_all.sh                       # all stages (SLiM skips if absent)
#    bash run_all.sh analytic math dfe     # just the Python theory layers
#    bash run_all.sh fitness               # just the climate-fitness figure
#    SEEDS="1 2 3" bash run_all.sh sim aggregate
# =============================================================================
set -uo pipefail
cd "$(dirname "$0")"
export ROOT="$(pwd)"
source "$ROOT/pipeline_lib.sh"

STAGES="${*:-env analytic math dfe sim aggregate panel calibration fitness}"
run_stages "$STAGES" "04_futureClimate/data_real"

hr "DONE (real data)"
echo "figures -> $FIGDIR/ ;  grid output -> $GRIDDIR/ ;  data -> $DATADIR/"
