import subprocess
import itertools
import csv
import os

# Parameter ranges
jump_sizes = [2.0, 3.0, 4.0, 5.0]
disk_multipliers = [1.0, 2.0, 3.0]
recent_windows = [3, 5, 7, 10]
track_names = ["track1", "track20"]

results = []
output_folder = "final_paths"
os.makedirs(output_folder, exist_ok=True)

path_id_counter = 1  # Start path IDs at 1

for track in track_names:
    for jump, disk_mul, recent in itertools.product(jump_sizes, disk_multipliers, recent_windows):
        cmd = ["./rrt", "--no-window", str(jump), str(disk_mul), str(recent), track]
        print(f"Running: Track={track}, Jump={jump}, DiskMul={disk_mul}, Recent={recent}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Parse output
        path_id = None
        for line in result.stdout.splitlines():
            if line.startswith("RESULT"):
                data = dict(item.split("=") for item in line.split(",")[1:])
                data = {k: float(v) if "." in v else int(v) for k, v in data.items()}
                data["Track"] = track  # Add track name
                results.append(data)
            elif line.startswith("PATH_VECTOR:"):
                path_id = f"path_{path_id_counter}.cpp"
                path_file = os.path.join(output_folder, path_id)
                with open(path_file, "w") as f:
                    f.write("#include \"../geometry.h\"\n\n")
                    f.write(line.split("PATH_VECTOR:")[1].strip())
                data["PathFile"] = path_id  # Add path ID to the result
                path_id_counter += 1  # Increment path ID

# Save to CSV
with open("rrt_test_results.csv", "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print("Testing complete. Results saved to rrt_test_results.csv and paths saved to final_paths folder.")