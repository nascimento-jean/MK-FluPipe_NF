#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SIF_DIR="$SCRIPT_DIR/sif"
TMP_DIR="$SCRIPT_DIR/.singularity-build"
mkdir -p "$SIF_DIR" "$TMP_DIR"

if command -v apptainer >/dev/null 2>&1; then
    SINGULARITY_BIN="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    SINGULARITY_BIN="singularity"
else
    echo "Neither apptainer nor singularity was found in PATH." >&2
    exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
    echo "docker was not found in PATH." >&2
    exit 1
fi

MK_FLU_TAR="$TMP_DIR/mk_flu_tools_local.tar"
MEDAKA_TAR="$TMP_DIR/medaka_tools_local.tar"

docker save -o "$MK_FLU_TAR" mk-flu-pipe/mk_flu_tools:local
docker save -o "$MEDAKA_TAR" mk-flu-pipe/medaka_tools:local

"$SINGULARITY_BIN" build "$SIF_DIR/mk_flu_tools_local.sif" "docker-archive://$MK_FLU_TAR"
"$SINGULARITY_BIN" build "$SIF_DIR/medaka_tools_local.sif" "docker-archive://$MEDAKA_TAR"

echo "Singularity/Apptainer images built successfully:"
echo "  $SIF_DIR/mk_flu_tools_local.sif"
echo "  $SIF_DIR/medaka_tools_local.sif"