import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("rrt_test_results.csv")

# Set style
sns.set_theme(style="whitegrid")

# Ensure data types
df["JumpSize"] = df["JumpSize"].astype(float)
df["DiskMultiplier"] = df["DiskMultiplier"].astype(float)
df["RecentWindow"] = df["RecentWindow"].astype(int)

# Pick constant values (median) for other parameters
constant_diskmult = 2
constant_jump = 3
constant_recent = 5

# Create folder for plots if it doesn't exist
output_folder = "plots"
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

for track in df["Track"].unique():
    # --- JumpSize plots ---
    df_jump = filter_df("JumpSize", track)
    # Runtime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_jump, x="JumpSize", y="RuntimeMs", marker="o")
    plt.title(f"Runtime vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Runtime (ms)")
    plt.savefig(os.path.join(output_folder, f"runtime_vs_jump_{track}.png"))
    plt.close()

    # AvgIterationTime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_jump, x="JumpSize", y="AvgIterationTimeMs", marker="o", color="orange")
    plt.title(f"Avg Iteration Time vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Avg Iteration Time (ms)")
    plt.savefig(os.path.join(output_folder, f"avg_iter_time_vs_jump_{track}.png"))
    plt.close()

    # Path Length
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_jump, x="JumpSize", y="PathLength", marker="o", color="green")
    plt.title(f"Path Length vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Path Length")
    plt.savefig(os.path.join(output_folder, f"pathlength_vs_jump_{track}.png"))
    plt.close()

    # Max Angle
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_jump, x="JumpSize", y="MaxAngleDeg", marker="o", color="purple")
    plt.title(f"Max Angle vs JumpSize (DiskMultiplier={constant_diskmult}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Max Angle (degrees)")
    plt.savefig(os.path.join(output_folder, f"max_angle_vs_jump_{track}.png"))
    plt.close()

    # --- DiskMultiplier plots ---
    df_disk = filter_df("DiskMultiplier", track)
    # Runtime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_disk, x="DiskMultiplier", y="RuntimeMs", marker="o")
    plt.title(f"Runtime vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Runtime (ms)")
    plt.savefig(os.path.join(output_folder, f"runtime_vs_diskmult_{track}.png"))
    plt.close()

    # AvgIterationTime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_disk, x="DiskMultiplier", y="AvgIterationTimeMs", marker="o", color="orange")
    plt.title(f"Avg Iteration Time vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Avg Iteration Time (ms)")
    plt.savefig(os.path.join(output_folder, f"avg_iter_time_vs_diskmult_{track}.png"))
    plt.close()

    # Path Length
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_disk, x="DiskMultiplier", y="PathLength", marker="o", color="green")
    plt.title(f"Path Length vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Path Length")
    plt.savefig(os.path.join(output_folder, f"pathlength_vs_diskmult_{track}.png"))
    plt.close()

    # Max Angle
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_disk, x="DiskMultiplier", y="MaxAngleDeg", marker="o", color="purple")
    plt.title(f"Max Angle vs DiskMultiplier (JumpSize={constant_jump}, RecentWindow={constant_recent}, Track={track})")
    plt.ylabel("Max Angle (degrees)")
    plt.savefig(os.path.join(output_folder, f"max_angle_vs_diskmult_{track}.png"))
    plt.close()

    # --- RecentWindow plots ---
    df_recent = filter_df("RecentWindow", track)
    # Runtime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_recent, x="RecentWindow", y="RuntimeMs", marker="o")
    plt.title(f"Runtime vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})")
    plt.ylabel("Runtime (ms)")
    plt.savefig(os.path.join(output_folder, f"runtime_vs_recentwindow_{track}.png"))
    plt.close()

    # AvgIterationTime
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_recent, x="RecentWindow", y="AvgIterationTimeMs", marker="o", color="orange")
    plt.title(f"Avg Iteration Time vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})")
    plt.ylabel("Avg Iteration Time (ms)")
    plt.savefig(os.path.join(output_folder, f"avg_iter_time_vs_recentwindow_{track}.png"))
    plt.close()

    # Path Length
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_recent, x="RecentWindow", y="PathLength", marker="o", color="green")
    plt.title(f"Path Length vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})")
    plt.ylabel("Path Length")
    plt.savefig(os.path.join(output_folder, f"pathlength_vs_recentwindow_{track}.png"))
    plt.close()

    # Max Angle
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=df_recent, x="RecentWindow", y="MaxAngleDeg", marker="o", color="purple")
    plt.title(f"Max Angle vs RecentWindow (JumpSize={constant_jump}, DiskMultiplier={constant_diskmult}, Track={track})")
    plt.ylabel("Max Angle (degrees)")
    plt.savefig(os.path.join(output_folder, f"max_angle_vs_recentwindow_{track}.png"))
    plt.close()

# Find the best paths
best_paths = df.loc[df.groupby("Track")["PathLength"].idxmin(), ["Track", "PathLength", "PathFile"]]

best_paths = best_paths.reset_index(drop=True)


# Print to console
print("Best path per track (shortest PathLength):")
print(best_paths)


print(f"Plots saved successfully in folder '{output_folder}'.")
