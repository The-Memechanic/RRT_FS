import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import shutil

# Load data
df = pd.read_csv("rrt_test_results.csv")

# Convert runtime values to milliseconds
df["RuntimeMs"] = df["RuntimeUs"] / 1000.0

# Ensure proper data types
df["JumpSize"] = df["JumpSize"].astype(float)
df["DiskMultiplier"] = df["DiskMultiplier"].astype(float)
df["RecentWindow"] = df["RecentWindow"].astype(int)

# Set style
sns.set_theme(style="whitegrid")

# Constants for filtering
constant_diskmult = 2
constant_jump = 3
constant_recent = 5

output_folder = "plots"

# Delete the previous outputs first
if os.path.exists(output_folder):
    shutil.rmtree(output_folder)

# Then recreate the folder
os.makedirs(output_folder, exist_ok=True)


def filter_df(variable, track):
    if variable == "JumpSize":
        return df[(df["DiskMultiplier"] == constant_diskmult) &
                  (df["RecentWindow"] == constant_recent) &
                  (df["Track"] == track)]
    elif variable == "DiskMultiplier":
        return df[(df["JumpSize"] == constant_jump) &
                  (df["RecentWindow"] == constant_recent) &
                  (df["Track"] == track)]
    elif variable == "RecentWindow":
        return df[(df["JumpSize"] == constant_jump) &
                  (df["DiskMultiplier"] == constant_diskmult) &
                  (df["Track"] == track)]

def plot_with_mean(data, x, y, title, ylabel, filename, color=None):
    data = data.copy()
    data[x] = pd.to_numeric(data[x], errors='coerce')
    data = data.dropna(subset=[x, y])
    data = data.sort_values(by=x)

    plt.figure(figsize=(10, 6))

    # Scatterplot for the individual points (respects numeric x)
    sns.scatterplot(data=data, x=x, y=y, color="gray", alpha=0.4)

    # Lineplot for the mean line
    sns.lineplot(data=data, x=x, y=y, errorbar='ci', marker="o", color=color)

    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel(x)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, filename))
    plt.close()



for track in df["Track"].unique():
    # --- JumpSize plots ---
    df_jump = filter_df("JumpSize", track)
    plot_with_mean(df_jump, "JumpSize", "RuntimeMs",
                   f"Runtime vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})",
                   "Runtime (ms)", f"runtime_vs_jump_{track}.png")

    plot_with_mean(df_jump, "JumpSize", "AvgIterationTimeUs",
                   f"Avg Iteration Time vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})",
                   "Avg Iteration Time (us)", f"avg_iter_time_vs_jump_{track}.png", color="orange")

    plot_with_mean(df_jump, "JumpSize", "PathLength",
                   f"Path Length vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})",
                   "Path Length", f"pathlength_vs_jump_{track}.png", color="green")

    plot_with_mean(df_jump, "JumpSize", "MaxAngleDeg",
                   f"Max Angle vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})",
                   "Max Angle (degrees)", f"max_angle_vs_jump_{track}.png", color="purple")

    # --- DiskMultiplier plots ---
    df_disk = filter_df("DiskMultiplier", track)
    plot_with_mean(df_disk, "DiskMultiplier", "RuntimeMs",
                   f"Runtime vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})",
                   "Runtime (ms)", f"runtime_vs_diskmult_{track}.png")

    plot_with_mean(df_disk, "DiskMultiplier", "AvgIterationTimeUs",
                   f"Avg Iteration Time vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})",
                   "Avg Iteration Time (us)", f"avg_iter_time_vs_diskmult_{track}.png", color="orange")

    plot_with_mean(df_disk, "DiskMultiplier", "PathLength",
                   f"Path Length vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})",
                   "Path Length", f"pathlength_vs_diskmult_{track}.png", color="green")

    plot_with_mean(df_disk, "DiskMultiplier", "MaxAngleDeg",
                   f"Max Angle vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})",
                   "Max Angle (degrees)", f"max_angle_vs_diskmult_{track}.png", color="purple")

    # --- RecentWindow plots ---
    df_recent = filter_df("RecentWindow", track)
    plot_with_mean(df_recent, "RecentWindow", "RuntimeMs",
                   f"Runtime vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})",
                   "Runtime (ms)", f"runtime_vs_recentwindow_{track}.png")

    plot_with_mean(df_recent, "RecentWindow", "AvgIterationTimeUs",
                   f"Avg Iteration Time vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})",
                   "Avg Iteration Time (us)", f"avg_iter_time_vs_recentwindow_{track}.png", color="orange")

    plot_with_mean(df_recent, "RecentWindow", "PathLength",
                   f"Path Length vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})",
                   "Path Length", f"pathlength_vs_recentwindow_{track}.png", color="green")

    plot_with_mean(df_recent, "RecentWindow", "MaxAngleDeg",
                   f"Max Angle vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})",
                   "Max Angle (degrees)", f"max_angle_vs_recentwindow_{track}.png", color="purple")

# Best paths
best_paths = df.loc[df.groupby("Track")["PathLength"].idxmin(), ["Track", "PathLength", "PathFile"]].reset_index(drop=True)

print("Best path per track (shortest PathLength):")
print(best_paths)
print(f"Plots saved successfully in folder '{output_folder}'.")
