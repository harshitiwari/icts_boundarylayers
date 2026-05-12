"""
Parameters for Mahrt (1972) Advective Boundary Layer Model
All configuration in one place
"""

import os
from pathlib import Path

# ==================== OUTPUT CONFIGURATION ====================
OUTPUT_DIR = "output"  # Main output directory
EXPERIMENT_NAME = "mahrt_2experiments"  # Subdirectory name

# ==================== SIMULATION PARAMETERS ====================
EXPERIMENT_1 = {
    'name': 'Experiment_1_meridional_PG_only',
    'description': 'Meridional pressure gradient only',
    'dpdy': -4.0e-5,          # Meridional pressure gradient from the paper
    'dpdx': 0.0,              # No zonal pressure gradient
    'ny': 41,
    'y_min': -8.0,
    'y_max': 12.0,
}

EXPERIMENT_2 = {
    'name': 'Experiment_2_both_PG',
    'description': 'Both zonal and meridional pressure gradients',
    'dpdy': -2.84e-5,
    'dpdx': -2.84e-5,
    'ny': 41,
    'y_min': -10.0,
    'y_max': 10.0,
}

# ==================== PHYSICAL CONSTANTS ====================
BETA = 2.28e-11           # Equatorial beta parameter (m^-1 s^-1)
K_EDDY = 5.0              # Constant eddy diffusivity used in the paper (m^2 s^-1)
ALPHA = 1.0               # Specific volume (constant)

# ==================== GRID PARAMETERS ====================
NZ = 20                   # Number of vertical levels
NY = 41                   # Number of meridional points (~50 km spacing)
Z_MAX = 4000.0            # Maximum height (m)
Y_MIN = -10.0             # Default minimum latitude (degrees)
Y_MAX = 10.0              # Default maximum latitude (degrees)

# Computed grid spacing
DZ = Z_MAX / (NZ - 1)
DY = (Y_MAX - Y_MIN) / (NY - 1)
DY_KM = DY * 111.0        # Convert to km

# ==================== TIME INTEGRATION ====================
DAYS_TOTAL = 16.0         # Total simulation days
DT = 500.0                # Time step (seconds)
SAVE_INTERVAL_SECONDS = 22 * 3600  # Save every 22 hours (approximately 1600 time steps)
SAVE_INTERVAL_STEPS = int(SAVE_INTERVAL_SECONDS / DT)  # In timesteps

# ==================== INITIAL CONDITIONS ====================
V_G_INIT = 0.0            # Fallback geostrophic wind magnitude (m/s)
EKMAN_DEPTH = 500.0       # Ekman layer depth (m)
EKMAN_SPIRAL_ANGLE = 3.14159 / 4  # 45 degrees
EQUATOR_REFERENCE_LAT = 4.0  # Paper: initialize the equatorial point using the 4° latitude Ekman solution

# ==================== OUTPUT FILE CONFIGURATION ====================
SAVE_FIELDS = True        # Save velocity fields to files
FIELD_SAVE_FORMAT = 'h5py'  # Format: 'h5py' (HDF5)
FIELDS_TO_SAVE = ['u', 'v', 'w', 'z', 'y']  # Which fields to save

# ==================== POSTPROCESSING CONFIGURATION ====================
PLOT_HODOGRAPHS = True
PLOT_STREAMLINES = True
PLOT_PROFILES = True
PLOT_QUIVER = True

# Figure DPI
FIGURE_DPI = 150

# ==================== HELPER FUNCTION ====================
def get_output_path():
    """Get full path to output directory"""
    base_path = Path(__file__).parent / OUTPUT_DIR / EXPERIMENT_NAME
    return str(base_path)

def get_experiment_output_path(exp_name):
    """Get path for specific experiment output"""
    path = Path(__file__).parent / OUTPUT_DIR / EXPERIMENT_NAME / exp_name
    return str(path)

def create_output_directories():
    """Create all necessary output directories"""
    base = Path(__file__).parent / OUTPUT_DIR / EXPERIMENT_NAME
    
    # Create main output directory
    base.mkdir(parents=True, exist_ok=True)
    
    # Create experiment subdirectories
    (base / EXPERIMENT_1['name']).mkdir(exist_ok=True)
    (base / EXPERIMENT_2['name']).mkdir(exist_ok=True)
    
    # Create subdirectories for fields, plots, logs
    for exp in [EXPERIMENT_1['name'], EXPERIMENT_2['name']]:
        exp_dir = base / exp
        (exp_dir / 'fields').mkdir(exist_ok=True)
        (exp_dir / 'plots').mkdir(exist_ok=True)
        (exp_dir / 'logs').mkdir(exist_ok=True)
    
    print(f"✓ Output directory structure created at: {base}")
    return base

# ==================== PRINT CONFIGURATION ====================
if __name__ == "__main__":
    print("MAHRT MODEL PARAMETERS")
    print("=" * 70)
    print(f"\nExperiment 1: {EXPERIMENT_1['name']}")
    print(f"  Description: {EXPERIMENT_1['description']}")
    print(f"  Pressure gradients: dpdy={EXPERIMENT_1['dpdy']:.2e}, dpdx={EXPERIMENT_1['dpdx']:.2e}")
    
    print(f"\nExperiment 2: {EXPERIMENT_2['name']}")
    print(f"  Description: {EXPERIMENT_2['description']}")
    print(f"  Pressure gradients: dpdy={EXPERIMENT_2['dpdy']:.2e}, dpdx={EXPERIMENT_2['dpdx']:.2e}")
    
    print(f"\nGrid Configuration:")
    print(f"  Vertical: {NZ} levels, 0-{Z_MAX}m, dz={DZ:.1f}m")
    print(f"  Meridional: {NY} points, {Y_MIN}°-{Y_MAX}°, dy={DY:.4f}° ({DY_KM:.2f}km)")
    
    print(f"\nTime Integration:")
    print(f"  Total: {DAYS_TOTAL} days")
    print(f"  Timestep: {DT}s")
    print(f"  Save interval: {SAVE_INTERVAL_STEPS} steps (~{SAVE_INTERVAL_SECONDS/3600:.1f} hours)")
    
    print(f"\nPhysical Constants:")
    print(f"  Beta: {BETA:.2e} m^-1 s^-1")
    print(f"  Eddy viscosity: {K_EDDY} m^2/s")
    
    print(f"\nOutput Configuration:")
    print(f"  Directory: {get_output_path()}")
    print(f"  Save fields: {SAVE_FIELDS}")
    print(f"  Format: {FIELD_SAVE_FORMAT}")
