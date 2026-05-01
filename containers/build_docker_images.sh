#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DOCKER_BUILDKIT=1 docker build \
  -t mk-flu-pipe/mk_flu_tools:local \
  "$SCRIPT_DIR/mk_flu_tools"

DOCKER_BUILDKIT=1 docker build \
  -t mk-flu-pipe/medaka_tools:local \
  "$SCRIPT_DIR/medaka_tools"

echo "Docker images built successfully:"
echo "  mk-flu-pipe/mk_flu_tools:local"
echo "  mk-flu-pipe/medaka_tools:local"