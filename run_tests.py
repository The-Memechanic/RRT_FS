import subprocess
import itertools
import csv

# Parameter ranges
jump_sizes = [2.0, 3.0, 4.0]
disk_multipliers = [1.0, 2.0, 3.0]
recent_windows = [3, 5, 10]
track_names = ["track20", "track1"]

results = []

for track in track_names:
    for jump, disk_mul, recent in itertools.product(jump_sizes, disk_multipliers, recent_windows):
        cmd = ["./rrt", "--no-window", str(jump), str(disk_mul), str(recent), track]
        print(f"Running: Track={track}, Jump={jump}, DiskMul={disk_mul}, Recent={recent}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Parse output
        for line in result.stdout.splitlines():
            if line.startswith("RESULT"):
                data = dict(item.split("=") for item in line.split(",")[1:])
                data = {k: float(v) if "." in v else int(v) for k, v in data.items()}
                data["Track"] = track  # Add track name
                results.append(data)

# Save to CSV
with open("rrt_test_results.csv", "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print("Testing complete. Results saved to rrt_test_results.csv")
