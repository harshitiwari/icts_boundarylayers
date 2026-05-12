"""
Main execution script for Mahrt (1972) Model
Reads parameters, runs solver, saves outputs
"""

import numpy as np
from pathlib import Path
import json

from parameters import (
    EXPERIMENT_1, EXPERIMENT_2, DT, DAYS_TOTAL,
    SAVE_INTERVAL_STEPS, SAVE_FIELDS, FIELD_SAVE_FORMAT,
    create_output_directories, get_experiment_output_path
)
from mahrt_solver import MahrtSolver


def save_fields(output_dir, snapshot_num, u_data, v_data, w_data, z_data, y_data, time_val, timestep):
    """Save velocity fields to HDF5 file with time in filename"""
    import h5py
    
    output_dir = Path(output_dir) / 'fields'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Include time info in filename
    time_days = time_val / 86400
    filename = output_dir / f'snapshot_{snapshot_num:04d}_t_{time_days:5.1f}d.h5'
    
    with h5py.File(filename, 'w') as f:
        # Save velocity fields with compression
        f.create_dataset('u', data=u_data, compression='gzip')
        f.create_dataset('v', data=v_data, compression='gzip')
        f.create_dataset('w', data=w_data, compression='gzip')
        
        # Save grid and time info
        f.create_dataset('z', data=z_data)
        f.create_dataset('y', data=y_data)
        
        # Save metadata as attributes
        f.attrs['time'] = time_val
        f.attrs['timestep'] = timestep
    
    return str(filename)


def save_metadata(solver, output_dir, exp_config, history, elapsed_time):
    """Save metadata about the run"""
    metadata = {
        'experiment': exp_config['name'],
        'description': exp_config['description'],
        'parameters': {
            'dpdy': float(solver.dpdy),
            'dpdx': float(solver.dpdx),
            'beta': float(solver.beta),
            'K_eddy': float(solver.K),
            'dz': float(solver.dz),
            'dy': float(solver.dy),
        },
        'grid': {
            'nz': int(solver.nz),
            'ny': int(solver.ny),
            'z_max': float(solver.z[-1]),
            'y_min': float(solver.y[0]),
            'y_max': float(solver.y[-1]),
        },
        'integration': {
            'dt': float(DT),
            'total_steps': int(history['steps'][-1]),
            'total_days': float(DAYS_TOTAL),
            'snapshots_saved': len(history['steps']),
            'final_time_seconds': float(history['times'][-1]),
        },
        'timing': {
            'elapsed_seconds': float(elapsed_time),
        }
    }
    
    meta_file = Path(output_dir) / 'metadata.json'
    with open(meta_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    return str(meta_file)


def run_experiment(exp_config, exp_num):
    """Run a single experiment"""
    
    import time
    start_time = time.time()
    
    print("\n" + "="*70)
    print(f"RUNNING: {exp_config['name']}")
    print(f"Description: {exp_config['description']}")
    print("="*70)
    
    # Create solver
    solver = MahrtSolver(
        dpdy=exp_config['dpdy'],
        dpdx=exp_config['dpdx'],
        experiment_name=f"Exp{exp_num}",
        ny=exp_config.get('ny'),
        y_min=exp_config.get('y_min'),
        y_max=exp_config.get('y_max')
    )
    
    # Calculate total timesteps
    total_steps = int(DAYS_TOTAL * 24 * 3600 / DT)
    
    # Integrate
    history = solver.integrate(
        dt=DT,
        total_steps=total_steps,
        save_steps=SAVE_INTERVAL_STEPS,
        verbose=True
    )
    
    elapsed_time = time.time() - start_time
    
    # Get output directory
    output_dir = get_experiment_output_path(exp_config['name'])
    
    # Save fields
    if SAVE_FIELDS:
        print(f"\nSaving fields to {output_dir}/fields/...")
        for i, (u, v, w, t) in enumerate(zip(history['u_history'], history['v_history'], history['w_history'], history['times'])):
            filename = save_fields(output_dir, i, u, v, w, solver.z, solver.y, t, history['steps'][i])
            print(f"  Saved snapshot {i}: {Path(filename).name}")
    
    # Save metadata
    meta_file = save_metadata(solver, output_dir, exp_config, history, elapsed_time)
    print(f"Saved metadata: {Path(meta_file).name}")
    
    print(f"\nTiming: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"✓ {exp_config['name']} completed!")
    
    return solver, history


def main():
    """Main execution"""
    
    import shutil
    
    print("\n" + "="*70)
    print("MAHRT (1972) ADVECTIVE BOUNDARY LAYER MODEL")
    print("Numerical Solver - Main Execution")
    print("="*70)
    
    # Clean old output
    old_output = Path(__file__).parent / 'output'
    if old_output.exists():
        print(f"\n🗑️  Removing old output: {old_output}")
        shutil.rmtree(old_output)
    
    # Create output directories
    create_output_directories()
    
    # Run Experiment 1
    solver1, history1 = run_experiment(EXPERIMENT_1, 1)
    
    # Run Experiment 2
    solver2, history2 = run_experiment(EXPERIMENT_2, 2)
    
    print("\n" + "="*70)
    print("ALL EXPERIMENTS COMPLETED SUCCESSFULLY")
    print("="*70)
    print(f"\nOutput saved to: {get_experiment_output_path('.')}")
    print("\nNext step: Run postprocess.py to generate plots")
    print("="*70 + "\n")
    
    return solver1, history1, solver2, history2


if __name__ == "__main__":
    solver1, history1, solver2, history2 = main()
