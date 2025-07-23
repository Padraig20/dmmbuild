import argparse
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt

pd.set_option("display.max_rows", 100)

def parse_profile_negln(path):
    with open(path) as f:
        lines = f.readlines()

    # gather header metadata until HMM column line
    meta = {}
    hmm_start = None
    for i, line in enumerate(lines):
        if line.startswith("HMM ") or line.startswith("DMM "):
            hmm_start = i
            break
        m = re.match(r"^(\S+)\s+(.*)$", line.rstrip())
        if m:
            key, val = m.group(1), m.group(2).strip()
            meta[key] = val

    if hmm_start is None:
        raise ValueError(f"'HMM ' block not found in {path}")

    trans_line = lines[hmm_start + 1]
    trans_labels = trans_line.strip().split()

    rows = []
    row_ids = []
    current_row_label = 1

    for line in lines[hmm_start+4::3]:

        s = line.strip()
        if not s:
            continue
        if line.startswith("//"):
            break

        vals = s.split()
        
        current_row_label = current_row_label + 1

        probs = []
        for v in vals:
            if v == "*":
                p = 0.0
            else:
                p = float(np.exp(-float(v)))
            probs.append(p)

        rows.append(probs)
        row_ids.append(current_row_label)

    df = pd.DataFrame(rows, columns=trans_labels, index=row_ids)
    if "COMPO" in df.index:
        df = df.drop(index="COMPO")
    try:
        df.index = df.index.astype(int)
    except Exception:
        pass
    return df, meta


def compare(dummer_file, hmmer_file, verbose=False):
    df_dummer, _ = parse_profile_negln(dummer_file)
    df_hmmer, _ = parse_profile_negln(hmmer_file)
    
    if verbose:
        print("DUMMER transition probabilities:")
        print(df_dummer.round(6))
        print("\nHMMER transition probabilities:")
        print(df_hmmer.round(6))
    
    cols = [c for c in df_hmmer.columns if c in df_dummer.columns]

    common_idx = df_dummer.index.intersection(df_hmmer.index)
    ddf = df_dummer.reindex(common_idx)[cols]
    hdf = df_hmmer.reindex(common_idx)[cols]

    diff = ddf - hdf
    mean_diff = diff.mean(axis=0, skipna=True)
    return mean_diff, diff, ddf, hdf

def plot_mean_diff(mean_diff):
    ax = mean_diff.plot(kind="bar")
    ax.axhline(0, linewidth=1)
    ax.set_ylabel("Mean Δp (DUMMER − HMMER)")
    ax.set_title("Mean Transition Probability Differences")
    plt.tight_layout()
    plt.show()

ap = argparse.ArgumentParser(
    description="Compare transition probabilities between DUMMER and HMMER .hmm profiles (-ln prob scale)."
)
ap.add_argument("dummer", help="DUMMER .hmm file")
ap.add_argument("hmmer", help="HMMER .hmm file")
ap.add_argument("--png", action="store_true", help="Plot bar plot PNG.")
ap.add_argument("--verbose", action="store_true", help="Print detailed output.")
args = ap.parse_args()

mean_diff, _, _, _ = compare(args.dummer, args.hmmer, verbose=args.verbose)

print("# Mean transition probability differences (DUMMER - HMMER)")
print(mean_diff.to_string(float_format=lambda x: f"{x:.6f}"))

if args.png:
    plot_mean_diff(mean_diff)
    
print("Done!")
