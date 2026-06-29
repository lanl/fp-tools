<div align="center">
    <h1>Welcome to FP-Tools!</h1>
</div>

[Quick Start Guide](#quick-start) \|
[Input/Output Formats](#inputoutput-formats) \|
[License/Copyright Assertion](#license)

![Release](https://img.shields.io/badge/Release-0.1.0-purple)

![til](/docs/imgs/demo.svg)

`fp-tools` is an open-source toolkit for analyzing molecular dynamics (MD) simulations using first-passage time diagnostics. It is designed to provide a fast, flexible, and extensible analysis framework capable of working with large-scale trajectory data from formats such as LAMMPS, VASP, and XYZ.

## 📈 Features

`fp-tools` supports several trajectory-based statistical analyses. Enable any of the following by specifying them in your configuration file under `[analysis]`:

* **Mean Squared Displacement (MSD)**\
Computes the average squared displacement of particles as a function of lag time (i.e., the time interval between two positions along a trajectory, over which displacements are measured and averaged).

* **Radial Distribution Function (RDF)**\
Measures pairwise spatial correlations between atom types as a function of interparticle distance (highlights how particle density varies relative to a reference particle).

* **Mean First Passage Time (MFPT) and D(r)**\
Computes the average time it takes a particle to first reach a distance $r$ from its starting point. The logarithmic derivative $D(r)$ captures how this time scales with distance.

* **Correlation Functions**\
Includes built-in options for common time correlation functions as well as a black-box interface where users can provide two raw data files to compute cross-correlations or autocorrelations directly.

---

## 📦 Quick Start <a name="quick-start"></a>

### Requirements

- **CMake 3.24+**
- **Python 3.8+**
- **C++17 compatible compiler** (e.g., GCC 8.1+, Clang 5.0+)

#### Python Dependencies
These are required to run the plotting scripts:
- `numpy`
- `matplotlib`

  > **Recommended setup:** 
  > Create and activate a virtual environment with Conda or Mamba:
  >   ```bash
  >   conda create -n fptools
  >   conda activate fptools
  >   pip install -r requirements.txt
  >   ```  

  > **Or install directly with pip:** 
  >   ```bash
  >   pip install -r requirements.txt
  >   ```

### Build and Install

If you have not downloaded ```fp-tools``` yet, please clone it from GitHub via

```bash
git clone ssh://git@re-git.lanl.gov:10022/fp-tools/fp-tools.git
```

From the base of the ```fp-tools``` source directory, execute:

```bash
# Configure project and dependencies
cmake -S . -B build

# Build project
cmake --build build
```

### Add the Executable to Your PATH

To run ```fp-tools``` from anywhere in your terminal, add the directory containing the executable (e.g., `build/bin`) to your PATH environment variable.

For example, if the executable is in ```$HOME/src/fp-tools/build/bin```, add the following line to your `~/.bashrc`, `~/.zshrc`, or appropriate shell
configuration file and reload your terminal or source the file:

```bash
export PATH="$HOME/src/fp-tools/build/bin:$PATH"
```
### Basic Usage

Run the tool with a TOML configuration file:

```bash
fptools config.toml
```
If you have **not** added the executable to your PATH, run it by specifying the full or relative path to the executable. For example, if your executable is located in the ```build/bin``` directory inside the project root, run:

```bash
./build/bin/fptools config.toml
```

or with a full absolute path:

```bash
$HOME/src/fp-tools/build/bin/fptools config.toml
```

## 📝 Input/Output Formats <a name="inputoutput-formats"></a>

`fp-tools` uses a **TOML** configuration file to specify input files, system parameters, and analysis options. You can find a sample configuration in the `examples/` directory to help you get started.

### Input Trajectories

The tool currently supports the following trajectory formats:

- **LAMMPS** (`.dump` or `.lammpstrj`)
- **VASP** (`XDATCAR`)
- **XYZ** (`.xyz`)

Each frame/timestep **must** include atomic positions. Box dimensions can be provided separately if needed (as in the case of VASP's XDATCAR file format). Optionally, velocities can be included to enable calculations such as correlation functions.

### Output Files

By default, all outputs are saved in the current working directory. For each enabled analysis, `fp-tools` generates:

- A `.dat` file containing raw numerical data,
- An optional plot for quick visualization,
- A detailed log output documenting the analysis parameters and progress.

##  License

This project (O#: O5103) is licensed under the BSD-3 license. See the [LICENSE](LICENSE) file for details.
