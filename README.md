# BORAT-EM

**BORAT-EM** (BlackOut RAy Tracer for Tracing Electromagnetic rays in plasmas) is a Python solver for electromagnetic (EM) wave propagation in ionized environments, based on geometrical optics ray tracing and high-frequency approximations.

---

## 📦 Installation

1. **Using Conda (recommended)**  
   Import the virtual environment provided in:
   ```
   env/BORAT-EM_Linux.yaml
   ```

2. **Or install dependencies manually**  
   Create your own virtual environment and install required packages:
   ```bash
   pip install pyvista scipy p-tqdm joblib
   ```

3. **TecIO Library (Optional for Tecplot output)**  
   Set the correct path to the TecIO shared library in `pytecio.py`.  
   If unavailable, deactivate TecIO output and rely on VTK output only.

---

## ▶️ Running BORAT

Main entry point:
```bash
python main.py
```

Edit the configuration file in `config/ExoMars.cfg` to change the test case or parameters.

---

## 🧭 Workflow Overview

1. **Read Configuration**
   - Loads all parameters from a `.cfg` file in the `config/` folder.

2. **Setup Domain**
   - Imports CFD mesh and associated fields.
   - Loads PEC and observation surfaces from `3DMeshes/`.

3. **Define EM Source**
   - Constructs antenna rays and ray tubes.

4. **Solve Eikonal Equation**
   - Propagates rays through the plasma medium.
   - Supports serial or parallel execution (joblib, multiprocessing, or p_tqdm).

5. **Save Ray Solutions**
   - Saves ray paths and tubes in `.vtk` and `.plk` formats.
   - Optionally exports `.szplt` via TecIO.

6. **Perform Aperture Integration**
   - Computes scattered field from rays using Physical Optics (PO) method.

7. **Post-process & Plot**
   - Generates polar plots of scattered field.
   - Outputs saved in `Output/<CaseName>/`.

---

## 📁 Folder Structure

```
BORAT-EM-open/
├── BORAT/              # Core ray-tracing modules
├── config/             # Configuration files (.cfg)
├── env/                # Conda environment YAML
├── 3DMeshes/           # VTK geometry files for boundaries (not provided)
├── CFD/                # CFD solutions (not provided)
├── Output/             # Simulation results
├── LICENSE
└── README.md           
```

---

## 📊 Visualization Options

### 🔹 ParaView
- Output: `.vtk` / `.vtm` in `Output/<CaseName>/`
- Visualize rays, tubes, and field data

### 🔹 Tecplot (if TecIO enabled)
- Output: `.szplt` in `Output/<CaseName>/`
- Requires TecIO library (`libtecio.so` or equivalent)

---

## 🔄 File Conversion Tips

1. **ParaView → Tecplot**
   - Export `.vtu` and convert using Tecplot import tools.

2. **Tecplot → ParaView**
   - Use `.plt` (2009 format) compatible with Paraview.

---

## 📝 License

This project is licensed under the [MIT License](LICENSE).

---

## 📣 Citation

This code supports an open-access scientific publication.  
If used in your work, please cite the corresponding article (DOI to be added).
