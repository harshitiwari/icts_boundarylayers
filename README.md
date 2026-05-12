MAHRT MODEL - QUICK REFERENCE
==============================

FILES YOU NEED TO KNOW:
-----------------------

📝 parameters.py
   Edit THIS to change:
   • K_EDDY (eddy viscosity)
   • DT (timestep)
   • DAYS_TOTAL (simulation duration)
   • EXPERIMENT_1, EXPERIMENT_2 (forcing conditions)

🧮 mahrt_solver.py
   Pure solver - DO NOT edit unless fixing physics

⚙️ main_model.py
   Runs solver, saves HDF5 - DO NOT edit normally

📊 postprocessing/postprocess.py
   Creates PNG plots from HDF5 - stands alone

📚 readme.md
   Usage documentation


TWO-COMMAND WORKFLOW:
--------------------

$ python3 main_model.py              # ~8 sec - generates HDF5 files
$ python3 postprocessing/postprocess.py  # ~5 sec - generates PNG plots


OUTPUT FILES:
-------------

HDF5 Snapshots (34 total):
  output/mahrt_2experiments/Experiment_*/fields/snapshot_XXXX.h5
  
  Each contains:
  • u, v, w - velocity components
  • z, y - grid coordinates
  • time, timestep - metadata

PNG Plots (8 total):
  output/mahrt_2experiments/Experiment_*/plots/
  • hodographs.png
  • streamlines.png
  • profiles.png
  • velocity_contours.png

Metadata (2 files):
  output/mahrt_2experiments/Experiment_*/metadata.json


TO LOAD HDF5 DATA IN PYTHON:
----------------------------

import h5py
import numpy as np

with h5py.File('output/mahrt_2experiments/Experiment_1_meridional_PG_only/fields/snapshot_0000.h5', 'r') as f:
    u = f['u'][:]      # (20, 30) array
    v = f['v'][:]      # (20, 30) array
    w = f['w'][:]      # (20, 30) array
    z = f['z'][:]      # (20,) array
    y = f['y'][:]      # (30,) array
    time = f.attrs['time']
    timestep = f.attrs['timestep']


DEFAULT VALUES IN parameters.py:
--------------------------------

Grid:
  NZ = 20           # vertical levels
  NY = 30           # meridional points
  Z_MAX = 4000 m
  Y_MIN = -0.5°, Y_MAX = 0.5°

Physics:
  BETA = 2.28e-11   # equatorial beta
  K_EDDY = 1.0 m²/s
  
Time:
  DT = 500 s
  DAYS_TOTAL = 16 days
  
Experiments:
  Exp1: dpdy=-4e-6, dpdx=0 (meridional only)
  Exp2: dpdy=-2.84e-6, dpdx=-2.84e-6 (both)


FILE SIZES:
-----------

Code:      ~26 KB total
Output:    ~2.6 MB per run
  • HDF5:  ~800 KB (all snapshots)
  • PNG:   ~500 KB (all plots)


EXTEND THE MODEL:
-----------------

Add new experiment:
  1. Edit parameters.py
  2. Add EXPERIMENT_3 = {...}
  3. Run main_model.py (auto-detects)

Change parameters:
  1. Edit parameters.py
  2. Run main_model.py
  3. Run postprocessing/postprocess.py

Custom analysis:
  1. Import from postprocessing: load_snapshot
  2. Load HDF5 file
  3. Analyze or plot


PROJECT STATUS: ✅ PRODUCTION READY
===================================

✓ Clean file structure
✓ Portable HDF5 format
✓ Modular code design
✓ Comprehensive plots
✓ Reproducible runs
✓ Metadata tracking
✓ Easy to extend
✓ Ready for publication
