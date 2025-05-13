# Autodock Deep Learning Toxicology GUI

## Installation and Setup

### Windows
- **Requirements**:
  - Autodock 4.2.6 must be installed.
  - Place the `autodock4.exe` file in the same directory as the program.
  - Include `python.exe` from `mgltools_win32_1.5.6` and the `mgltools` directory in the same directory as the program.
  - LaTeX must be installed and functional from the Anaconda prompt.
- **Notes**:
  - The starting package was built for Linux.
  - The Autodock and MGLTools setup differs for Windows.

### Linux
- **Requirements**:
  - Autodock4 must be installed and accessible via the command line (`autodock4` command).
  - Place the `mgltools_x86_64Linux2_1.5.6` directory in the same directory as the program.
  - Ensure `mpirun` is installed (it may come with `mpi4py`).
- **Installation Command**:
  ```bash
  sudo apt install autodock
  ```

### macOS
- **Requirements**:
  - Autodock4 must be installed and accessible via the command line (`autodock4` command).
  - LaTeX must be installed and functional.
  - `mpirun` must also be installed.
- **Notes**:
  - The Autodock4 binary for macOS must be downloaded from the Autodock4 website.
  - MGLTools is not required for macOS.

---

## Environment Setup

1. Use the `docking_gui.yml` file to set up Python and install required packages:
   ```bash
   conda env create --file docking_gui.yml
   ```
2. Activate the environment:
   ```bash
   conda activate deepdocktox
   ```
3. The models are availabe in the release. This git hub has one tag, click on that on the right to get the models
   https://github.com/mmccarthy1/Autodock_deeplearning_toxicology_GUI/releases/tag/v1.0
  
5. Run the program:
   ```bash
   python Autodock_deeplearning_toxicology_GUI.py
   ```
   - Ensure to use `--file` and not `-f` during installation to avoid issues with pip packages.

---

## Quick Start Instructions

1. Launch the program with:
   ```bash
   python Autodock_deeplearning_toxicology_GUI.py
   ```
2. Navigate to the **File** menu and select **1. Add Docking Project**. Enter a project name.
3. Populate the table:
   - In the **Compound** column, enter a molecule name, SMILES, or InChI code.
   - Ensure the **Name** column contains a corresponding unique name (no special characters or spaces, underscores are allowed).
   - Alternatively, paste a two-column list (comma-delimited) into the first two columns.
4. Select **Input molecules from table** from the **File** menu.
   - The program will search for compounds via PubChem and display "Job done" when complete.
5. Set up the docking project:
   - Choose **2. Setup Docking Project** from the **File** menu.
   - A `list-to-dock.txt` file will appear when completed.
6. Run the docking process:
   - Select **3. Run docking** from the **File** menu and specify the number of processors.
   - The process can take several hours and will display "Job Done" upon completion.
7. Process the docking results:
   - Choose **4. Process Docking**. This will run through 9 steps.
   - Note: Step 4 may seem complete but will continue in the background. Do not exit the program.
8. Once complete, navigate to the **Docking results and machine learning predictions** tab:
   - Click **Populate list** to display molecule names and receptor names.
   - Select molecule and receptor pairs to analyze.
   - Use "Remove from list" to refine the selections.
9. Generate a PDF report:
   - Requires LaTeX.
   - Click **Print report of displayed molecules** to save a PDF.

---

## Alternate Input Methods

1. **Generate 3D Structures**:
   - Use OpenBabel or RDKit to generate 3D structures from SMILES or InChI codes.
   - Replace Step 3 in the Quick Start and:
     - Place a text file in the `molecules` folder with SMILES or InChI codes.
     - Select **Generate 3D molecules** from the **File** menu and follow prompts.
   - The table will update with molecule information and 2D images.
2. **Use an SDF File**:
   - Skip Step 3 and select **Input molecules by 3D SDF file** from the **File** menu.
   - Continue with Step 5.

---

## Troubleshooting

- **OpenMPI Issues**:
  - If multiple docking runs are launched, this may indicate duplicate OpenMPI installations.
  - Uninstall `mpi4py` and reinstall it with:
    ```bash
    pip install mpi4py --no-deps
    ```

---

## Prerequisites

### Non-Python Installs
- **AutoDock**: 4.2.6 **NOTE**: this must be downloaded separately
- **MGLTools**: 1.5.6 (Linux and Windows only) **NOTE**: this must be downloaded separately
- **MPI**: v3 or v4 (OpenMPI, Intel OpenMPI, or Microsoft SDK version)
- **LaTeX**: Tested with LaTeX Live
- **CUDA**: v9.0
- **NVIDIA-Driver**: 570.124.04

### Python Packages
- **Python**: v3.8
- **Conda**: v4.9.2
- **binana.py** (from Durrant lab and University of Pittsburgh)
- **TensorFlow**: v2.5.0
- **Pandas**: v1.2.3
- **NumPy**: v1.19
- **Keras**: v2.4.3
- **mpi4py**: v3.0.3
- **scypy**: v1.10.1
- **OpenBabel**: v3.1.1 (verified version)
- **Plotly**: v4.9.3 (newer versions also tested)
- **Kaleido**
- **Matplotlib**
- **PyLaTeX**

---

