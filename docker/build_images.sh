#!/usr/bin/env bash
# Build the MutaXplorer images (one per conda environment) from envs/*.yml
set -euo pipefail
cd "$(dirname "$0")/.."
TAG=${1:-1.0.0}
for env in mutaxplorer cnvkit manta arriba; do
  name=$([[ $env == mutaxplorer ]] && echo mutaxplorer || echo mutaxplorer-$env)
  docker build -t "${name}:${TAG}" --build-arg ENV_FILE="envs/${env}.yml" -f Dockerfile .
done
