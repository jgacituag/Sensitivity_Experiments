"""
Generate a Saltelli (Sobol) design CSV ready for WRF runs and later sensitivity analysis.
"""
import os
import pandas as pd
from SALib.sample import sobol

# ==== CONFIGURE HERE ====
d = 5                 # Number of parameters (x1..xd)
max_total = 1000      # Total rows <= this (first-order); change to 2000, etc.
second_order = False  # True -> include second-order terms (bigger design)
seed = 42             # Fix for reproducibility; change/remove for a new draw
out_dir = "data"
out_path = os.path.join(out_dir, f"saltelli_d{d}_{'2nd' if second_order else '1st'}_le{max_total}.csv")
# ========================

def main():
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # SALib "problem" in [0,1] for each parameter
    problem = {
        "num_vars": d,
        "names": [f"x{i}" for i in range(1, d+1)],
        "bounds": [[0.0, 1.0]] * d
    }

    # Choose base N so total rows <= max_total (formula depends on first/second order)
    denom = (2 * d + 2) if second_order else (d + 2)
    N = max(1, max_total // denom)
    print(f"Chosen base N={N} to get total rows <= {max_total}")

    # Build the stacked Saltelli matrix
    X = sobol.sample(problem, N, calc_second_order=second_order, seed=seed)

    # Wrap in a DataFrame
    df = pd.DataFrame(X, columns=problem["names"])

    # Save to CSV
    df.to_csv(out_path, index=False)
    print(f"d={d}, second_order={second_order}, max_total={max_total}")
    print(f"Base N={N} -> total_rows={len(df)} (expected {N*denom})")
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    main()