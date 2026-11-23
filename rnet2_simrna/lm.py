# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Karim Malekzadeh
"""Base pairs from RibonanzaNet2: one attention head, symmetrised and APC-corrected, decoded greedily.

The head (layer 25, head 1; flat index 288 of 48 x 12) was chosen once, as the head whose attention best matched known base
pairs (average precision) on 612 training RNAs: SPOT-RNA PDB TR1 (112 after filtering) + 500 random bpRNA TR0. It was never
tuned on the test targets. Note: the training pool contains 4LVW, a 99% copy of 6Q57.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import torch

from .ribonanzanet import RibonanzaNet, load_config_from_yaml

TOKENS = {"A": 0, "C": 1, "G": 2, "U": 3}
CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
MIN_SEP = 4                  # |i - j| >= 4
PAIR_FRACTION = 0.3          # predict round(0.3 L) pairs
HEAD_LAYER, HEAD = 24, 0     # 0-based: layer 25, head 1  (flat index 24 * 12 + 0 = 288)


def load_model(ckpt_dir: str | Path, device: str) -> RibonanzaNet:
    """ckpt_dir holds pairwise.yaml and pytorch_model_fsdp.bin (Kaggle shujun717/ribonanzanet2, see scripts/get_weights.sh)."""
    ckpt_dir = Path(ckpt_dir)
    m = RibonanzaNet(load_config_from_yaml(str(ckpt_dir / "pairwise.yaml")))
    m.load_state_dict(torch.load(ckpt_dir / "pytorch_model_fsdp.bin", map_location="cpu"), strict=True)
    return m.to(device).eval()


@torch.no_grad()
def head_attention(model: RibonanzaNet, seq: str, device: str) -> torch.Tensor:
    """Attention map (L, L) of the chosen head. Layers after HEAD_LAYER do not affect it, so they are skipped."""
    src = torch.tensor([[TOKENS[c] for c in seq]], device=device); mask = torch.ones_like(src)
    x = model.encoder(src); z = model.outer_product_mean(x) + model.pos_encoder(x)
    for li, layer in enumerate(model.transformer_encoder):
        x, z, aw = layer([x, z, mask.clone(), True])
        if li == HEAD_LAYER:
            return aw[0][HEAD].float()
    raise RuntimeError("model has fewer layers than expected")


def sym_apc(a: torch.Tensor) -> torch.Tensor:
    a = 0.5 * (a + a.transpose(-1, -2))
    r = a.mean(-1, keepdim=True); c = a.mean(-2, keepdim=True); t = a.mean((-1, -2), keepdim=True)
    return a - r * c / (t + 1e-8)


def greedy(scores: np.ndarray, iu, seq: str) -> set[tuple[int, int]]:
    """Canonical pairs in score order, one partner per nucleotide, round(0.3 L) pairs. Crossing pairs are allowed."""
    L = len(seq); n = round(PAIR_FRACTION * L); used = np.zeros(L, bool); pred = set()
    for k in np.argsort(-scores, kind="stable"):
        if len(pred) >= n:
            break
        i, j = int(iu[0][k]), int(iu[1][k])
        if used[i] or used[j] or seq[i] + seq[j] not in CANON:
            continue
        pred.add((i, j)); used[i] = used[j] = True
    return pred


def predict_pairs(model: RibonanzaNet, seq: str, device: str) -> list[tuple[int, int]]:
    """0-based (i, j) pairs, i < j."""
    iu = np.triu_indices(len(seq), MIN_SEP)
    s = sym_apc(head_attention(model, seq, device))[torch.as_tensor(iu[0]), torch.as_tensor(iu[1])].cpu().numpy()
    return sorted(greedy(s, iu, seq))


def layers(pairs, L: int) -> list[list[tuple[int, int]]]:
    """Split pairs into non-crossing layers; the first layer is the main nested structure, the rest are pseudoknots."""
    rest = sorted(pairs); out = []
    while rest:
        lay, left = [], []
        for i, j in rest:
            (left if any(a < i < b < j or i < a < j < b for a, b in lay) else lay).append((i, j))
        out.append(lay); rest = left
    return out
