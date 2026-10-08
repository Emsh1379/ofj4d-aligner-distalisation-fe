#!/usr/bin/env bash
# Build the FEBioAlignerPairs plug-in.  usage: FEBIO=$HOME/febio/FEBio bash pipeline/febio_pairs/build.sh
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
FEBIO=${FEBIO:-$HOME/febio/FEBio}
cmake -S "$HERE" -B "$HERE/build" -DCMAKE_BUILD_TYPE=Release -DFEBIO_ROOT="$FEBIO" -DFEBIO_LIB_DIR="$FEBIO/build/lib"
cmake --build "$HERE/build" -j "$(nproc)"
echo "plug-in: $HERE/build/libFEBioAlignerPairs.so"
