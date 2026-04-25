#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

IN_DIR="${PROJECT_DIR}/results/quantification"
OUT_DIR="${PROJECT_DIR}/results/visualization"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

TPM_FILE="${IN_DIR}/tpm_matrix.tsv"

[[ -f "$TPM_FILE" ]] || {
    echo "[ERROR] Missing TPM matrix: $TPM_FILE"
    echo "Run quantification.sh first."
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "[ERROR] python3 not found"
    exit 1
}

echo "Creating figures from: $TPM_FILE"

export TPM_FILE OUT_DIR

python3 <<'PYEOF'
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

tpm_file = os.environ["TPM_FILE"]
out_dir = os.environ["OUT_DIR"]

df = pd.read_csv(tpm_file, sep="\t")

if "Geneid" not in df.columns:
    print("[ERROR] Expected column Geneid not found", file=sys.stderr)
    sys.exit(1)

sample_cols = [c for c in df.columns if c not in ["Geneid", "Length"]]

for c in sample_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

df["mean_TPM"] = df[sample_cols].mean(axis=1)
top = df.sort_values("mean_TPM", ascending=False).head(20).copy()

top.to_csv(os.path.join(out_dir, "top20_tpm_used.tsv"), sep="\t", index=False)

plt.figure(figsize=(12, 7))
plt.bar(top["Geneid"], top["mean_TPM"])
plt.xticks(rotation=75, ha="right")
plt.ylabel("Mean TPM")
plt.xlabel("Gene")
plt.title("Top 20 expressed genes by mean TPM")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "barplot_top20_mean_tpm.png"), dpi=300)
plt.close()

heat = top.set_index("Geneid")[sample_cols].copy()
heat = np.log2(heat + 1)

fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(heat.values, aspect="auto")

ax.set_xticks(range(len(sample_cols)))
ax.set_xticklabels(sample_cols, rotation=45, ha="right")
ax.set_yticks(range(len(heat.index)))
ax.set_yticklabels(heat.index)

ax.set_xlabel("Sample")
ax.set_ylabel("Gene")
ax.set_title("Heatmap of top 20 expressed genes\nlog2(TPM + 1)")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("log2(TPM + 1)")

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "heatmap_top20_tpm.png"), dpi=300)
plt.close()

print("Saved:")
print("barplot_top20_mean_tpm.png")
print("heatmap_top20_tpm.png")
print("top20_tpm_used.tsv")
PYEOF

echo ""
echo "Visualization complete."
echo "Outputs saved in: $OUT_DIR"
