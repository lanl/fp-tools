#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def msd_plot(input_files, labels, output_dir, output_format="png", group_index=0):
    if len(input_files) != len(labels):
        raise ValueError("Number of input files and labels must match.")

    plt.figure()

    for file_path, label in zip(input_files, labels):
        times, values = [], []
        with open(file_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    times.append(float(parts[0]))
                    values.append(float(parts[1]))

        plt.loglog(times, values, label=label)

    plt.xlabel(r"$t$", fontsize=14)
    plt.ylabel("MSD", fontsize=14)
    plt.grid(True, which="both", ls="--", linewidth=0.3)
    plt.legend()
    #plt.tight_layout()

    # Use group index in filename to keep order and uniqueness
    plot_filename = f"figure_{group_index}.{output_format}"
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, format=output_format)
    plt.close()

def plot_mfpt(mfpt_file, output_dir, output_format="png"):
    """Plot Mean First Passage Time (MFPT) from data file."""
    r, avg_time = [], []
    with open(mfpt_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            r.append(float(parts[0]))
            avg_time.append(float(parts[1]))

    plt.figure(figsize=(6, 4))
    plt.loglog(r, avg_time, linestyle='-', color='tab:blue')
    plt.xlabel(r"$r$")
    plt.ylabel("Mean First Passage Time (MFPT)")
    plt.grid(True, linestyle='--', alpha=0.3)

    base_name = os.path.splitext(os.path.basename(mfpt_file))[0]
    plot_path = os.path.join(output_dir, f"{base_name}.{output_format}")
    #plt.tight_layout()
    plt.savefig(plot_path, format=output_format)

def plot_log_derivative(log_file, output_dir, output_format="png"):
    """Plot D(r) log-derivative data from file."""
    data = np.loadtxt(log_file, comments="#")
    r, dlog = data[:, 0], data[:, 1]

    plt.figure(figsize=(6, 4))
    #plt.xlim(0, 5)
    plt.plot(r, dlog, linestyle='-', color='tab:blue', label='D(r)')
    plt.xlabel(r"$r/a$")
    plt.ylabel(r"D(r)")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()

    base_name = os.path.splitext(os.path.basename(log_file))[0]
    plot_path = os.path.join(output_dir, f"{base_name}.{output_format}")
    plt.tight_layout()
    plt.savefig(plot_path, format=output_format)

def rdf_plot(input_files, labels, output_dir, output_format="png", group_index=0):
    import os
    import matplotlib.pyplot as plt

    plt.figure()

    for file_path, label in zip(input_files, labels):
        r, gr = [], []
        with open(file_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 2:
                    r.append(float(parts[0]))
                    gr.append(float(parts[1]))

        plt.plot(r, gr, label=label)

    plt.xlabel(r"$r$", fontsize=20)
    plt.ylabel(r"$g(r)$", fontsize=20)
    plt.grid(True, which="both", ls="--", alpha=0.3)
    plt.legend()
    #plt.tight_layout()

    plot_path = os.path.join(output_dir, f"figure_{group_index}.{output_format}")
    plt.savefig(plot_path, format=output_format)
    plt.close()

def plot_corr(filepath, output_dir, output_format="png"):
    """Plot a correlation function from a .dat file."""
    data = np.loadtxt(filepath, comments="#")
    lag, corr = data[:, 0], data[:, 1]

    plt.figure(figsize=(6, 4))
    #plt.xlim(0, 5)
    plt.plot(lag, corr, linestyle='-', color='tab:blue', label=os.path.basename(filepath).replace(".dat", ""))
    plt.xlabel("Lag")
    plt.ylabel("Correlation")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend()

    base_name = os.path.splitext(os.path.basename(filepath))[0]
    plot_path = os.path.join(output_dir, f"{base_name}.{output_format}")
    plt.tight_layout()
    plt.savefig(plot_path, format=output_format)
    plt.close()

def main():
    def usage():
        print("Usage:")
        print("  python3 fp_plot.py msd <group_index> <file1> <label1> [<file2> <label2> ...] <output_dir> <format>")
        print("  python3 fp_plot.py rdf <group_index> <file1> <label1> [<file2> <label2> ...] <output_dir> <format>")
        print("  python3 fp_plot.py mfpt <mfpt_file> <output_dir> [format]")
        print("  python3 fp_plot.py dlog <log_derivative_file> <output_dir> [format]")
        print("  python3 fp_plot.py corr <output_dir> <format>")


    if len(sys.argv) < 2:
        usage()
        return

    cmd = sys.argv[1]
    args = sys.argv[2:]

    if cmd == "msd":
        # Expect first argument after 'msd' to be group index
        if len(args) < 5 or (len(args) - 3) % 2 != 0:
            usage()
            return
        group_index = int(args[0])
        n = (len(args) - 3) // 2
        files = []
        labels = []
        for i in range(n):
            files.append(args[1 + 2*i])
            labels.append(args[2 + 2*i])
        output_dir = args[-2]
        output_format = args[-1]
        msd_plot(files, labels, output_dir, output_format, group_index)

    elif cmd == "mfpt" and len(args) >= 2:
        plot_mfpt(args[0], args[1], args[2] if len(args) > 2 else "png")

    elif cmd == "dlog" and len(args) >= 2:
        plot_log_derivative(args[0], args[1], args[2] if len(args) > 2 else "png")

    elif cmd == "rdf":
        # args: group_index, <file1> <label1> [file2 label2 ...], output_dir, format
        if len(args) < 5 or (len(args) - 3) % 2 != 0:
            usage()
            return
        group_index = int(args[0])
        n = (len(args) - 3) // 2
        files = []
        labels = []
        for i in range(n):
            files.append(args[1 + 2*i])
            labels.append(args[2 + 2*i])
        output_dir = args[-2]
        output_format = args[-1]
        rdf_plot(files, labels, output_dir, output_format, group_index)

    elif cmd == "corr" and len(args) >= 2:
        input_dir = args[0]
        output_format = args[1]
        corr_files = sorted(f for f in os.listdir(input_dir) if f.endswith("corr.dat"))
        for file in corr_files:
            plot_corr(os.path.join(input_dir, file), input_dir, output_format)


    else:
        print("Invalid command or insufficient arguments.\n")
        usage()

if __name__ == "__main__":
    main()
