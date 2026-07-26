# Plot analysis scaling results from native Rust Pegasus and OpenAeroStruct.
#
# Usage:
#   cargo run --release --example analysis_scaling \
#       > benchmark_scaling_out/results_rust.csv
#   python benchmark_analysis_oas.py
#   python plot_analysis_scaling.py

import numpy as np
import matplotlib.pyplot as plt

rust = np.genfromtxt("benchmark_scaling_out/results_rust.csv",
                      delimiter=",", skip_header=1)
oas = np.genfromtxt("benchmark_scaling_out/results_oas.csv",
                     delimiter=",", skip_header=1)

# Figure 1: Wall time per analysis vs mesh size
fig1, ax1 = plt.subplots(figsize=(7, 5))
ax1.plot(rust[:, 1], rust[:, 2], "o-", color="tab:blue", label="Pegasus (Rust)")
ax1.plot(oas[:, 1], oas[:, 2], "s-", color="tab:orange", label="OpenAeroStruct")
ax1.set_xlabel("Total panels (full-span)")
ax1.set_ylabel("Wall time per analysis [ms]")
ax1.set_title("Single Analysis Time vs Mesh Size")
ax1.legend()
ax1.grid(True, alpha=0.3)
fig1.tight_layout()
fig1.savefig("benchmark_scaling_out/analysis_time.png", dpi=150)

# Figure 2: Same data, log-log
fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.plot(rust[:, 1], rust[:, 2], "o-", color="tab:blue", label="Pegasus (Rust)")
ax2.plot(oas[:, 1], oas[:, 2], "s-", color="tab:orange", label="OpenAeroStruct")
ax2.set_xlabel("Total panels (full-span)")
ax2.set_ylabel("Wall time per analysis [ms]")
ax2.set_title("Single Analysis Time vs Mesh Size (log-log)")
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xscale("log")
ax2.set_yscale("log")
fig2.tight_layout()
fig2.savefig("benchmark_scaling_out/analysis_time_loglog.png", dpi=150)

plt.show()
print("Plots saved to benchmark_scaling_out/")
