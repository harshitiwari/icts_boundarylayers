"""
Postprocessing and Visualization for Mahrt (1972) Model
Generate plots from saved fields
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from parameters import (
    get_experiment_output_path, FIGURE_DPI, FIELD_SAVE_FORMAT,
    EXPERIMENT_1, EXPERIMENT_2
)


def load_snapshot(exp_dir, snapshot_num=0, format_type='h5py'):
    """Load a snapshot from disk (latest by default, or by number)"""
    import h5py
    import glob
    
    field_dir = Path(exp_dir) / 'fields'
    
    if format_type == 'h5py':
        # Find all h5 files and sort them
        h5_files = sorted(glob.glob(str(field_dir / '*.h5')))
        
        if not h5_files:
            raise FileNotFoundError(f"No HDF5 files found in {field_dir}")
        
        # Load the requested snapshot (default: last/final state)
        if snapshot_num < 0:
            snapshot_num = len(h5_files) + snapshot_num
        
        filename = h5_files[min(snapshot_num, len(h5_files)-1)]
        
        with h5py.File(filename, 'r') as f:
            return {
                'u': f['u'][:],
                'v': f['v'][:],
                'w': f['w'][:],
                'z': f['z'][:],
                'y': f['y'][:],
                'time': float(f.attrs['time']),
                'timestep': int(f.attrs['timestep']),
                'filename': Path(filename).name,
            }
    
    # Legacy NPZ support
    elif format_type == 'npz':
        filename = field_dir / f'snapshot_{snapshot_num:04d}.npz'
        data = np.load(filename)
        return {
            'u': data['u'],
            'v': data['v'],
            'w': data['w'],
            'z': data['z'],
            'y': data['y'],
            'time': float(data['time']),
            'timestep': int(data['timestep']),
        }


def _select_lat_indices(y_values, targets=(-6.0, 0.0, 6.0)):
    """Pick representative latitude indices closest to the requested targets."""
    indices = []
    for target in targets:
        idx = int(np.argmin(np.abs(y_values - target)))
        if idx not in indices:
            indices.append(idx)
    return indices


def _compute_terms_from_snapshot(u, v, z, y, beta, k_eddy, dpdx=0.0, dpdy=0.0):
    """Recompute the momentum terms used by the paper's depth criterion."""
    z = np.asarray(z)
    y = np.asarray(y)
    dz = float(z[1] - z[0]) if len(z) > 1 else 1.0
    dy = float((y[1] - y[0]) * 111000.0) if len(y) > 1 else 1.0
    y_meters = y * 111000.0

    w = np.zeros_like(v)
    for iz in range(len(z) - 1):
        dvdy = np.zeros(len(y))
        dvdy[1:-1] = (v[iz, 2:] - v[iz, :-2]) / (2 * dy)
        dvdy[0] = (v[iz, 1] - v[iz, 0]) / dy
        dvdy[-1] = (v[iz, -1] - v[iz, -2]) / dy
        w[iz + 1, :] = w[iz, :] - dvdy * dz

    adv_u = np.zeros_like(u)
    adv_v = np.zeros_like(v)
    for iz in range(len(z)):
        for iy in range(len(y)):
            if v[iz, iy] > 0 and iy > 0:
                dudy = (u[iz, iy] - u[iz, iy - 1]) / dy
                dvdy = (v[iz, iy] - v[iz, iy - 1]) / dy
            elif v[iz, iy] < 0 and iy < len(y) - 1:
                dudy = (u[iz, iy + 1] - u[iz, iy]) / dy
                dvdy = (v[iz, iy + 1] - v[iz, iy]) / dy
            else:
                dudy = 0.0
                dvdy = 0.0

            if w[iz, iy] > 0 and iz > 0:
                dudz = (u[iz, iy] - u[iz - 1, iy]) / dz
                dvdz = (v[iz, iy] - v[iz - 1, iy]) / dz
            elif w[iz, iy] < 0 and iz < len(z) - 1:
                dudz = (u[iz + 1, iy] - u[iz, iy]) / dz
                dvdz = (v[iz + 1, iy] - v[iz, iy]) / dz
            else:
                dudz = 0.0
                dvdz = 0.0

            adv_u[iz, iy] = -v[iz, iy] * dudy - w[iz, iy] * dudz
            adv_v[iz, iy] = -v[iz, iy] * dvdy - w[iz, iy] * dvdz

    diff_u = np.zeros_like(u)
    diff_v = np.zeros_like(v)
    for iz in range(1, len(z) - 1):
        diff_u[iz, :] = k_eddy * (u[iz + 1, :] - 2 * u[iz, :] + u[iz - 1, :]) / (dz ** 2)
        diff_v[iz, :] = k_eddy * (v[iz + 1, :] - 2 * v[iz, :] + v[iz - 1, :]) / (dz ** 2)

    coriolis_u = np.zeros_like(u)
    coriolis_v = np.zeros_like(v)
    for iz in range(len(z)):
        coriolis_u[iz, :] = beta * y_meters * v[iz, :]
        coriolis_v[iz, :] = -beta * y_meters * u[iz, :]

    pressure_u = np.full_like(u, -dpdx)
    pressure_v = np.full_like(v, -dpdy)

    return {
        'adv_u': adv_u,
        'adv_v': adv_v,
        'diff_u': diff_u,
        'diff_v': diff_v,
        'cor_u': coriolis_u,
        'cor_v': coriolis_v,
        'press_u': pressure_u,
        'press_v': pressure_v,
    }


def _paper_depth_profile(u, v, z, y, beta=2.28e-11, k_eddy=5.0, threshold=0.20, dpdx=0.0, dpdy=0.0):
    """Paper-style boundary-layer depth: first level where diffusion exceeds a fraction of max local term."""
    terms = _compute_terms_from_snapshot(u, v, z, y, beta, k_eddy, dpdx=dpdx, dpdy=dpdy)
    diff_mag = np.sqrt(terms['diff_u'] ** 2 + terms['diff_v'] ** 2)
    adv_mag = np.sqrt(terms['adv_u'] ** 2 + terms['adv_v'] ** 2)
    cor_mag = np.sqrt(terms['cor_u'] ** 2 + terms['cor_v'] ** 2)
    press_mag = np.sqrt(terms['press_u'] ** 2 + terms['press_v'] ** 2)

    depth = np.zeros(len(y))
    ratio = np.zeros((len(z), len(y)))
    for iy in range(len(y)):
        local_max = np.maximum.reduce([
            diff_mag[:, iy],
            adv_mag[:, iy],
            cor_mag[:, iy],
            press_mag[:, iy],
        ])
        local_max = np.where(local_max <= 0.0, np.nan, local_max)
        ratio[:, iy] = diff_mag[:, iy] / local_max
        above = np.where(ratio[:, iy] >= threshold)[0]
        depth[iy] = z[above[0]] if len(above) else z[-1]

    return depth, ratio, terms


def plot_hodographs(exp_dir, exp_name):
    """Plot wind hodographs at 3 latitudes in a single plot"""
    
    # Load final snapshot
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    z = data['z']
    y = data['y']
    
    # Select 3 latitudes: south of equator, near equator, north of equator
    lat_indices = _select_lat_indices(y)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # blue, orange, green
    lat_labels = [f'{y[i]:.2f}°' for i in lat_indices]
    
    fig, ax = plt.subplots(figsize=(11, 10))
    
    for idx, label, color in zip(lat_indices, lat_labels, colors):
        u_col = u[:, idx]
        v_col = v[:, idx]
        
        # Plot the hodograph spiral
        ax.plot(u_col, v_col, color=color, linewidth=3.0, label=label, alpha=0.85, zorder=5)
        
        # Mark surface (circle) and top (square)
        ax.plot(u_col[0], v_col[0], 'o', color=color, markersize=12, 
               markeredgecolor='black', markeredgewidth=1.5, zorder=10, alpha=0.9)
        ax.plot(u_col[-1], v_col[-1], 's', color=color, markersize=12, 
               markeredgecolor='black', markeredgewidth=1.5, zorder=10, alpha=0.9)
        
        # Add height annotations at regular intervals
        for i in range(0, len(z), max(1, len(z)//6)):
            ax.plot(u_col[i], v_col[i], 'k.', markersize=6, zorder=8, alpha=0.6)
            ax.annotate(f'{z[i]:.0f}m', (u_col[i], v_col[i]),
                       fontsize=8, xytext=(8, 8), textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor=color, alpha=0.2),
                       color=color, fontweight='bold')
    
    ax.set_xlabel('Zonal Wind u (m/s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Meridional Wind v (m/s)', fontsize=12, fontweight='bold')
    ax.set_title(f'{exp_name}: Wind Hodographs at Three Latitudes\n(Circles=Surface, Squares=Top)', 
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)
    ax.set_aspect('equal')
    ax.legend(fontsize=11, loc='best', title='Latitude', title_fontsize=12, framealpha=0.95)
    
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'hodographs.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_term_profiles(exp_dir, exp_name):
    """Plot signed components of momentum terms vs latitude at three heights.
    
    Shows u and v components separately for advection, diffusion, Coriolis, and pressure terms.
    Signed values reveal direction and magnitude of each force contribution.
    """
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    z = data['z']
    y = data['y']
    
    if exp_name == EXPERIMENT_1['name']:
        dpdx = EXPERIMENT_1['dpdx']
        dpdy = EXPERIMENT_1['dpdy']
    else:
        dpdx = EXPERIMENT_2['dpdx']
        dpdy = EXPERIMENT_2['dpdy']
    
    terms = _compute_terms_from_snapshot(u, v, z, y, beta=2.28e-11, k_eddy=5.0, 
                                         dpdx=dpdx, dpdy=dpdy)
    
    # Select 3 representative heights
    z_indices = [0, len(z)//2, -1]  # surface, mid, top
    z_labels = [f'{z[i]:.0f}m' for i in z_indices]
    
    # Create figure with 4 rows (one per term: adv, diff, cor, press) and 2 cols (u and v)
    fig, axes = plt.subplots(4, 2, figsize=(14, 12))
    fig.suptitle(f'{exp_name}: Signed Force Components vs Latitude\n(Positive = acceleration in that direction)', 
                 fontsize=14, fontweight='bold', y=0.995)
    
    term_names = ['Advection', 'Diffusion', 'Coriolis', 'Pressure']
    term_keys = [('adv_u', 'adv_v'), ('diff_u', 'diff_v'), ('cor_u', 'cor_v'), ('press_u', 'press_v')]
    
    colors = ['blue', 'green', 'red']  # for z_indices
    
    for row, (term_name, (key_u, key_v)) in enumerate(zip(term_names, term_keys)):
        # Column 0: u component
        ax = axes[row, 0]
        for z_idx, color, z_label in zip(z_indices, colors, z_labels):
            u_term = terms[key_u][z_idx, :]
            ax.plot(y, u_term, color=color, linewidth=2.0, marker='o', 
                   markersize=4, label=z_label, alpha=0.8)
        ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_ylabel(f'{term_name}\nu-component (m/s²)', fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
        if row == 0:
            ax.set_title('Zonal (u) Component', fontsize=11, fontweight='bold')
        
        # Column 1: v component
        ax = axes[row, 1]
        for z_idx, color, z_label in zip(z_indices, colors, z_labels):
            v_term = terms[key_v][z_idx, :]
            ax.plot(y, v_term, color=color, linewidth=2.0, marker='s', 
                   markersize=4, label=z_label, alpha=0.8)
        ax.axhline(0, color='k', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_ylabel(f'{term_name}\nv-component (m/s²)', fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
        if row == 0:
            ax.set_title('Meridional (v) Component', fontsize=11, fontweight='bold')
    
    # x-label on bottom row
    for ax in axes[-1, :]:
        ax.set_xlabel('Latitude (degrees)', fontsize=11, fontweight='bold')
    
    # Mark equator
    for ax in axes.flat:
        if np.min(y) < 0.0 < np.max(y):
            ax.axvline(0.0, color='gray', linestyle=':', linewidth=1.0, alpha=0.6)
    
    plt.tight_layout()
    out = Path(exp_dir) / 'plots' / 'terms_vs_latitude.png'
    plt.savefig(out, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {out.name}")
    return str(out)


def plot_streamlines(exp_dir, exp_name):
    """Plot streamlines in y-z plane"""
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    w = data.get('w', np.zeros_like(v))
    z = data['z']
    y = data['y']

    fig, ax = plt.subplots(figsize=(12, 8))

    # For streamplot: x axis = latitude (y), y axis = height (z)
    # streamplot expects arrays shaped (ny, nx) where nx=len(x), ny=len(y)
    # Our arrays are (nz, ny) so pass x,y 1D and u=meridional(v), v_comp=vertical(w)
    speed = np.sqrt(v**2 + w**2)

    strm = ax.streamplot(y, z, v, w,
                         color=speed, cmap='viridis', density=1.5)

    cbar = plt.colorbar(strm.lines, ax=ax, label='Wind Speed (m/s)')
    
    ax.set_xlabel('Latitude (degrees)', fontsize=12)
    ax.set_ylabel('Height (m)', fontsize=12)
    ax.set_title(f'{exp_name}: Streamlines in y-z Plane', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'streamlines.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_vertical_profiles(exp_dir, exp_name):
    """Plot vertical profiles at different latitudes"""
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    w = data['w']
    z = data['z']
    y = data['y']
    
    # Select 3 latitudes: south of equator, near equator, north of equator
    lat_indices = _select_lat_indices(y)
    colors = ['blue', 'red', 'green']
    labels = [f'{y[i]:.2f}°' for i in lat_indices]
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: u profiles
    ax = axes[0, 0]
    for idx, color, label in zip(lat_indices, colors, labels):
        ax.plot(u[:, idx], z, marker='o', linewidth=2.5, color=color, label=label, markersize=4)
    ax.set_xlabel('Zonal Wind u (m/s)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Height (m)', fontsize=11, fontweight='bold')
    ax.set_title('Vertical Profile: Zonal Wind', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: v profiles
    ax = axes[0, 1]
    for idx, color, label in zip(lat_indices, colors, labels):
        ax.plot(v[:, idx], z, marker='s', linewidth=2.5, color=color, label=label, markersize=4)
    ax.set_xlabel('Meridional Wind v (m/s)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Height (m)', fontsize=11, fontweight='bold')
    ax.set_title('Vertical Profile: Meridional Wind', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Speed profile
    ax = axes[1, 0]
    speed = np.sqrt(u**2 + v**2)
    for idx, color, label in zip(lat_indices, colors, labels):
        ax.plot(w[:, idx], z, marker='^', linewidth=2.5, color=color, label=label, markersize=4)
    ax.set_xlabel('Wind Speed (m/s)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Height (m)', fontsize=11, fontweight='bold')
    ax.set_title('Vertical Profile: Total Wind Speed', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Direction profile
    ax = axes[1, 1]
    angle = np.degrees(np.arctan2(v, u))
    for idx, color, label in zip(lat_indices, colors, labels):
        ax.plot(angle[:, idx], z, marker='d', linewidth=2.5, color=color, label=label, markersize=4)
    ax.set_xlabel('Wind Direction (degrees from East)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Height (m)', fontsize=11, fontweight='bold')
    ax.set_title('Vertical Profile: Wind Direction\n(0°=East, 90°=North)', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'profiles.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_quiver(exp_dir, exp_name):
    """Plot wind vectors as a single quiver diagram for the full domain"""
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    w = data.get('w', np.zeros_like(v))
    z = data['z']
    y = data['y']
    
    fig, ax = plt.subplots(figsize=(13, 8))
    
    # Build mesh for latitude (y) and height (z)
    Y, Z = np.meshgrid(y, z)
    # Use meridional (v) and vertical (w) components for y-z plane
    speed = np.sqrt(v**2 + w**2)
    
    cf = ax.contourf(Y, Z, speed, levels=16, cmap='YlOrRd', alpha=0.75)
    plt.colorbar(cf, ax=ax, label='Wind Speed (m/s)')
    
    # Subsample for readability
    skip_z = max(1, len(z) // 10)
    skip_y = max(1, len(y) // 18)
    Q = ax.quiver(
        Y[::skip_z, ::skip_y], Z[::skip_z, ::skip_y],
        v[::skip_z, ::skip_y], w[::skip_z, ::skip_y],
        scale=14, scale_units='inches', width=0.005,
        headwidth=4.5, headlength=6.0, headaxislength=5.5,
        pivot='tail', alpha=0.9, color='navy'
    )
    ax.quiverkey(Q, 0.92, 0.95, 1, '1 m/s', labelpos='E', coordinates='axes')
    
    # Mark the equator and domain center
    if np.min(y) < 0.0 < np.max(y):
        ax.axvline(0.0, color='white', linestyle='--', linewidth=1.5, alpha=0.9, label='Equator')
    ax.set_xlabel('Latitude (degrees)', fontsize=11)
    ax.set_ylabel('Height (m)', fontsize=11)
    ax.set_title(f'{exp_name}: Single Quiver Plot', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.25, linewidth=0.5)
    ax.legend(fontsize=9, loc='upper right')
    
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'quiver.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_velocity_contours(exp_dir, exp_name):
    """Plot velocity field contours"""
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    z = data['z']
    y = data['y']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    Y, Z = np.meshgrid(y, z)
    
    # Zonal wind
    cf1 = ax1.contourf(Y, Z, u, levels=15, cmap='RdBu_r')
    ax1.contour(Y, Z, u, levels=8, colors='k', alpha=0.3, linewidths=0.5)
    plt.colorbar(cf1, ax=ax1, label='Zonal Wind u (m/s)')
    ax1.set_title('Zonal Wind Component', fontsize=12)
    ax1.set_xlabel('Latitude (degrees)')
    ax1.set_ylabel('Height (m)')
    
    # Meridional wind
    cf2 = ax2.contourf(Y, Z, v, levels=15, cmap='RdBu_r')
    ax2.contour(Y, Z, v, levels=8, colors='k', alpha=0.3, linewidths=0.5)
    plt.colorbar(cf2, ax=ax2, label='Meridional Wind v (m/s)')
    ax2.set_title('Meridional Wind Component', fontsize=12)
    ax2.set_xlabel('Latitude (degrees)')
    ax2.set_ylabel('Height (m)')
    
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'velocity_contours.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_paper_depth(exp_dir, exp_name):
    """Plot the paper-style boundary-layer depth criterion versus latitude"""

    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    z = data['z']
    y = data['y']

    if exp_name == EXPERIMENT_1['name']:
        dpdx = EXPERIMENT_1['dpdx']
        dpdy = EXPERIMENT_1['dpdy']
    else:
        dpdx = EXPERIMENT_2['dpdx']
        dpdy = EXPERIMENT_2['dpdy']

    depth, ratio, _ = _paper_depth_profile(u, v, z, y, dpdx=dpdx, dpdy=dpdy)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(y, depth, color='black', linewidth=2.5)
    ax.axvline(0.0, color='gray', linestyle='--', linewidth=1.2, alpha=0.8)
    ax.set_xlabel('Latitude (degrees)', fontsize=11)
    ax.set_ylabel('Boundary-layer depth (m)', fontsize=11)
    ax.set_title(f'{exp_name}: Paper-style depth criterion', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plot_file = Path(exp_dir) / 'plots' / 'paper_depth.png'
    plt.tight_layout()
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def plot_ekman_spiral_3d(exp_dir, exp_name):
    """Create a 3D plot showing the Ekman spiral at selected latitudes"""
    
    from mpl_toolkits.mplot3d import Axes3D
    
    data = load_snapshot(exp_dir, snapshot_num=-1, format_type=FIELD_SAVE_FORMAT)
    u = data['u']
    v = data['v']
    z = data['z']
    y = data['y']
    
    # Select 3 latitudes to show spiral pattern
    lat_indices = _select_lat_indices(y)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # blue, orange, green
    labels = [f'{y[i]:.1f}°' for i in lat_indices]
    
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot spirals for each latitude
    for lat_idx, color, label in zip(lat_indices, colors, labels):
        u_profile = u[:, lat_idx]
        v_profile = v[:, lat_idx]
        z_profile = z
        
        # Plot the main spiral curve
        ax.plot(u_profile, v_profile, z_profile, 
                color=color, linewidth=3.5, label=label, marker='o', 
                markersize=6, markevery=max(1, len(z_profile)//8), 
                alpha=0.9, zorder=10)
        
        # Add arrow at top showing upwind direction
        if len(z_profile) >= 2:
            z_top_idx = -1
            z_prev_idx = -2
            ax.quiver(u_profile[z_prev_idx], v_profile[z_prev_idx], z_profile[z_prev_idx],
                     u_profile[z_top_idx] - u_profile[z_prev_idx],
                     v_profile[z_top_idx] - v_profile[z_prev_idx],
                     z_profile[z_top_idx] - z_profile[z_prev_idx],
                     color=color, arrow_length_ratio=0.25, linewidth=3.0, 
                     alpha=1.0, zorder=20)
        
        # Mark surface wind point
        ax.scatter([u_profile[0]], [v_profile[0]], [z_profile[0]], 
                  color=color, s=150, marker='s', alpha=0.8, 
                  edgecolors='black', linewidths=1.5, zorder=15)
    
    # Labels and formatting
    ax.set_xlabel('Zonal Wind u (m/s)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Meridional Wind v (m/s)', fontsize=12, fontweight='bold')
    ax.set_zlabel('Height (m)', fontsize=12, fontweight='bold')
    ax.set_title(f'{exp_name}: 3D Ekman Spiral (Wind Hodograph)\nVelocity vectors spiral with height due to Coriolis', 
                 fontsize=14, fontweight='bold', pad=25)
    
    # Adjust viewing angle for best 3D visualization
    ax.view_init(elev=25, azim=60)
    
    ax.legend(fontsize=10, loc='upper left', title='Latitude', title_fontsize=11, framealpha=0.95)
    ax.grid(True, alpha=0.4, linewidth=0.8)
    
    # Adjust margins for better view
    plt.tight_layout()
    plot_file = Path(exp_dir) / 'plots' / 'ekman_spiral_3d.png'
    plt.savefig(plot_file, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Saved: {plot_file.name}")
    return str(plot_file)


def postprocess_experiment(exp_dir, exp_name):
    """Generate all plots for an experiment"""
    
    print(f"\nPostprocessing: {exp_name}")
    print(f"  Output directory: {exp_dir}")
    
    # Create plots directory
    Path(exp_dir, 'plots').mkdir(exist_ok=True)
    
    # Generate all plots
    plot_hodographs(exp_dir, exp_name)
    plot_streamlines(exp_dir, exp_name)
    plot_vertical_profiles(exp_dir, exp_name)
    plot_quiver(exp_dir, exp_name)
    plot_velocity_contours(exp_dir, exp_name)
    plot_paper_depth(exp_dir, exp_name)
    plot_term_profiles(exp_dir, exp_name)
    plot_ekman_spiral_3d(exp_dir, exp_name)
    
    print(f"✓ {exp_name} postprocessing complete!")


def main():
    """Main postprocessing execution"""
    
    print("\n" + "="*70)
    print("MAHRT (1972) MODEL - POSTPROCESSING")
    print("="*70)
    
    base_dir = Path(__file__).parent.parent / 'output' / 'mahrt_2experiments'
    
    if not base_dir.exists():
        print(f"ERROR: Output directory not found: {base_dir}")
        print("Run main_model.py first to generate output files!")
        return
    
    # Postprocess both experiments
    exp1_dir = base_dir / EXPERIMENT_1['name']
    exp2_dir = base_dir / EXPERIMENT_2['name']
    
    postprocess_experiment(str(exp1_dir), EXPERIMENT_1['name'])
    postprocess_experiment(str(exp2_dir), EXPERIMENT_2['name'])
    
    print("\n" + "="*70)
    print("POSTPROCESSING COMPLETE")
    print(f"Plots saved in: {base_dir}")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
