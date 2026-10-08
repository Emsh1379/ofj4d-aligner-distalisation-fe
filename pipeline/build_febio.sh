#!/usr/bin/env bash
# Build FEBio (https://github.com/febiosoftware/FEBio) from source with Intel MKL
# PARDISO, for an external cross-validation of this model against an established
# biomechanics solver.  Needs git, cmake >= 3.10, g++ with C++17, python3 with pip.
#
# usage: bash pipeline/build_febio.sh [WORKDIR]        (default: $HOME/febio)
#   NO_MKL=1   build without MKL (FEBio's slower built-in solvers only)
#   JOBS=N     parallel compile jobs (default: all cores)
#
# Result: WORKDIR/febio.sh, a wrapper that sets the library path and runs FEBio,
#   e.g.  WORKDIR/febio.sh -i model.feb
set -euo pipefail
WORK=$(mkdir -p "${1:-$HOME/febio}" && cd "${1:-$HOME/febio}" && pwd)
JOBS=${JOBS:-$(nproc)}
cd "$WORK"
[ -d FEBio ] || git clone --depth 1 https://github.com/febiosoftware/FEBio.git

CMAKE_ARGS=(-DCMAKE_BUILD_TYPE=Release)
MKL="$WORK/mkl"
if [ "${NO_MKL:-0}" != 1 ]; then
  if [ ! -f "$MKL/include/mkl.h" ] || [ ! -f "$MKL/lib/libmkl_core.a" ]; then
    # FEBio links MKL statically (it only searches lib*.a on Linux).  pip --target
    # drops wheel data files, which is where the libraries live, so unpack by hand.
    # MKL < 2025: oneMKL 2025+ no longer declares the legacy sparse routines FEBio's
    # ILU0 preconditioner calls (mkl_dcsrtrsv); with 2026.1 the build stops in NumCore.
    mkdir -p "$MKL/wheels" "$MKL/include" "$MKL/lib"
    python3 -m pip download --no-deps -d "$MKL/wheels" "mkl-include<2025" "mkl-static<2025" "intel-openmp<2025"
    for w in "$MKL"/wheels/*.whl; do python3 -m zipfile -e "$w" "$MKL/unpacked"; done
    for d in "$MKL"/unpacked/*.data/data/include; do cp -r "$d"/. "$MKL/include/"; done
    for d in "$MKL"/unpacked/*.data/data/lib; do cp -r "$d"/. "$MKL/lib/"; done
  fi
  CMAKE_ARGS+=(-DUSE_MKL=ON -DMKL_INC="$MKL/include" -DMKL_LIB_DIR="$MKL/lib"
               -DMKL_OMP_LIB="$MKL/lib/libiomp5.so")
else
  CMAKE_ARGS+=(-DUSE_MKL=OFF)
fi

cmake -S FEBio -B FEBio/build "${CMAKE_ARGS[@]}"
cmake --build FEBio/build -j "$JOBS"

BIN=$(find "$WORK/FEBio/build/bin" -maxdepth 1 -type f -name 'febio*' -perm -u+x | head -n 1)
[ -n "$BIN" ] || { echo "build finished but no febio executable found in FEBio/build/bin"; exit 1; }
cat > "$WORK/febio.sh" <<EOF
#!/usr/bin/env bash
export LD_LIBRARY_PATH="$MKL/lib:$WORK/FEBio/build/lib:\${LD_LIBRARY_PATH:-}"
exec "$BIN" "\$@"
EOF
chmod +x "$WORK/febio.sh"
echo "FEBio built: $BIN"
echo "run it with:  $WORK/febio.sh -i model.feb"
