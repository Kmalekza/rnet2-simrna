#!/bin/bash
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Karim Malekzadeh
# Download the RibonanzaNet2 checkpoint (Kaggle model shujun717/ribonanzanet2, public, no login needed) into
# weights/ribonanzanet2/ (pairwise.yaml + pytorch_model_fsdp.bin, ~400 MB).
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
D=${RN2_DIR:-$REPO/weights/ribonanzanet2}
mkdir -p "$D"
if [ -f "$D/pytorch_model_fsdp.bin" ]; then echo "already present: $D"; exit 0; fi
curl -L --fail -o "$D/ribonanzanet2.tar.gz" https://www.kaggle.com/api/v1/models/shujun717/ribonanzanet2/pyTorch/alpha/1/download
tar -xzf "$D/ribonanzanet2.tar.gz" -C "$D" && rm "$D/ribonanzanet2.tar.gz"
ls -la "$D"
