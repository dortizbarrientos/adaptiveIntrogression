#!/usr/bin/env bash
# =============================================================================
# setup_python_env.sh
#
# One command to (re)build a reproducible Python environment for the
# adaptiveIntrogression analysis pipeline, then prove every package imports in
# the SAME interpreter you'll actually run. Solves the "I have tskit but the
# script can't see it" problem by pinning ONE venv and one activation step.
#
# USAGE
#   ./setup_python_env.sh            # create/refresh ./.venv and install deps
#   ./setup_python_env.sh --scan     # ALSO scan your .py files for real imports
#   ./setup_python_env.sh --recreate # delete and rebuild the venv from scratch
#
# After it finishes, activate the env in your shell with:
#   source .venv/bin/activate
# and run your scripts as usual (python3 aggregate_pl_sensitivity.py ...).
#
# DESIGN NOTES (so the script is legible, not magic):
#   * `set -euo pipefail` makes the script fail fast and loud: -e exit on error,
#     -u error on unset variables, -o pipefail catch errors mid-pipe.
#   * The venv is IDEMPOTENT: re-running reuses it unless --recreate is passed.
#   * The smoke test imports each package by its IMPORT name (not pip name) and
#     prints a per-package OK/MISSING table -- so a partial install gives a
#     precise diagnosis instead of a mystery.
#   * On success it writes requirements.lock.txt (exact versions) -- reuse that
#     for byte-for-byte reproducible environments on other machines.
# =============================================================================

set -euo pipefail

# ---- configuration ----------------------------------------------------------
VENV_DIR=".venv"
REQ_FILE="requirements.txt"
LOCK_FILE="requirements.lock.txt"
PYTHON_BIN="${PYTHON_BIN:-python3}"   # override with: PYTHON_BIN=python3.12 ./setup_python_env.sh

SCAN=0
RECREATE=0
for arg in "$@"; do
  case "$arg" in
    --scan)     SCAN=1 ;;
    --recreate) RECREATE=1 ;;
    *) echo "Unknown option: $arg"; echo "Use --scan and/or --recreate."; exit 2 ;;
  esac
done

echo "=============================================================="
echo " Python environment setup  (adaptiveIntrogression)"
echo "=============================================================="

# ---- 0. sanity: interpreter exists and is a sane version --------------------
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: '$PYTHON_BIN' not found on PATH." >&2
  echo "       Set PYTHON_BIN to your interpreter, e.g.:" >&2
  echo "       PYTHON_BIN=/opt/homebrew/bin/python3.12 ./setup_python_env.sh" >&2
  exit 1
fi
PYVER="$("$PYTHON_BIN" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
echo "  interpreter : $PYTHON_BIN  (Python $PYVER)"
echo "  venv dir    : $VENV_DIR"

# ---- 1. (optional) scan real imports before building anything ---------------
if [[ "$SCAN" -eq 1 ]]; then
  echo
  echo "--- import scan: third-party modules found in *.py here -----------------"
  # Grep import lines, strip to top-level module, drop stdlib + local names.
  # This is a HEURISTIC aid, not a resolver: it shows what your code imports so
  # you can reconcile requirements.txt against reality.
  STDLIB="$("$PYTHON_BIN" - <<'PY'
import sys
mods = getattr(sys, "stdlib_module_names", None)
print(" ".join(sorted(mods)) if mods else "")
PY
)"
  # CRITICAL: prune .venv and other dot-dirs, or we scan thousands of INSTALLED
  # library files instead of just the project's own scripts. The leading module
  # name must start with a letter or underscore (not a digit), and we keep only
  # the first dotted component (e.g. `from os.path` -> `os`).
  found="$(find . -name '*.py' -not -path './.venv/*' -not -path './.*/*' 2>/dev/null \
            | xargs grep -hoE '^[[:space:]]*(import|from)[[:space:]]+[a-zA-Z_][a-zA-Z0-9_]*' 2>/dev/null \
            | sed -E 's/^[[:space:]]*(import|from)[[:space:]]+//' \
            | sort -u || true)"
  if [[ -z "$found" ]]; then
    echo "  (no .py files found in $(pwd))"
  else
    while IFS= read -r mod; do
      [[ -z "$mod" ]] && continue
      # skip stdlib
      if [[ " $STDLIB " == *" $mod "* ]]; then continue; fi
      # skip obvious local modules (a matching file exists here)
      if [[ -f "./${mod}.py" ]]; then continue; fi
      echo "    third-party? $mod"
    done <<< "$found"
    echo "  -> reconcile the above against $REQ_FILE (pip name may differ from import name)."
  fi
fi

# ---- 2. create or refresh the venv ------------------------------------------
echo
if [[ "$RECREATE" -eq 1 && -d "$VENV_DIR" ]]; then
  echo "  --recreate: removing existing $VENV_DIR"
  rm -rf "$VENV_DIR"
fi
if [[ ! -d "$VENV_DIR" ]]; then
  echo "  creating venv ..."
  "$PYTHON_BIN" -m venv "$VENV_DIR"
else
  echo "  reusing existing venv (pass --recreate to rebuild)"
fi

# Activate for the rest of THIS script. (User still activates in their shell.)
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
echo "  active python: $(command -v python)"

# ---- 3. install -------------------------------------------------------------
echo
echo "  upgrading pip ..."
python -m pip install --upgrade pip --quiet

if [[ ! -f "$REQ_FILE" ]]; then
  echo "ERROR: $REQ_FILE not found in $(pwd)." >&2
  exit 1
fi
echo "  installing from $REQ_FILE ..."
# We do NOT use --quiet here: if a build fails (e.g. no wheel for a new Python),
# you want to see which package and why. We don't abort the script on a single
# failure -- the smoke test below reports exactly what made it in.
set +e
python -m pip install -r "$REQ_FILE"
PIP_RC=$?
set -e
[[ "$PIP_RC" -ne 0 ]] && echo "  NOTE: pip returned $PIP_RC -- see smoke test for what installed."

# ---- 4. smoke test: import each expected package in THIS interpreter --------
echo
echo "--- import smoke test ---------------------------------------------------"
# Map: "pip_name:import_name". Edit if you change requirements.txt.
python - <<'PY'
checks = [
    ("numpy", "numpy"),
    ("pandas", "pandas"),
    ("matplotlib", "matplotlib"),
    ("tskit", "tskit"),
    ("pyslim", "pyslim"),
    ("msprime", "msprime"),
]
import importlib
width = max(len(p) for p, _ in checks)
all_core_ok = True
core = {"numpy", "pandas", "matplotlib"}
for pip_name, imp in checks:
    try:
        m = importlib.import_module(imp)
        ver = getattr(m, "__version__", "?")
        print(f"    {pip_name:<{width}}  OK      {ver}")
    except Exception as e:
        flag = "MISSING (core!)" if pip_name in core else "missing (optional)"
        print(f"    {pip_name:<{width}}  {flag}")
        if pip_name in core:
            all_core_ok = False
import sys
sys.exit(0 if all_core_ok else 3)
PY
SMOKE_RC=$?

# ---- 5. lock exact versions on success --------------------------------------
echo
if [[ "$SMOKE_RC" -eq 0 ]]; then
  python -m pip freeze > "$LOCK_FILE"
  echo "  wrote $LOCK_FILE (exact versions for reproducible reuse)."
  echo
  echo "  SUCCESS. Activate the env in your shell before running scripts:"
  echo "      source $VENV_DIR/bin/activate"
else
  echo "  Core packages missing -- env is NOT ready. See the table above." >&2
  echo "  Common cause on bleeding-edge Python: no prebuilt wheel yet." >&2
  echo "  Try a slightly older interpreter, e.g.:" >&2
  echo "      PYTHON_BIN=python3.12 ./setup_python_env.sh --recreate" >&2
  exit 3
fi
