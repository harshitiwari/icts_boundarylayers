"""
Mahrt (1972) Advective Boundary Layer Model - CORE SOLVER
"""

import numpy as np
from parameters import (
    BETA, K_EDDY, ALPHA, NZ, NY, Z_MAX, Y_MIN, Y_MAX,
    DZ, DY, V_G_INIT, EKMAN_DEPTH, EKMAN_SPIRAL_ANGLE, EQUATOR_REFERENCE_LAT
)


class MahrtSolver:
    """2D numerical solver for Mahrt (1972) advective boundary layer model"""
    
    def __init__(self, dpdy, dpdx, experiment_name="Exp", ny=None, y_min=None, y_max=None, nz=None, z_max=None):
        """
        Initialize solver
        
        Parameters:
        -----------
        dpdy : float
            Meridional pressure gradient (m/s^2)
        dpdx : float
            Zonal pressure gradient (m/s^2)
        experiment_name : str
            Name of experiment for logging
        """
        
        self.experiment_name = experiment_name
        
        # Physical parameters
        self.beta = BETA
        self.K = K_EDDY
        self.alpha = ALPHA
        self.dpdy = dpdy
        self.dpdx = dpdx
        
        # Grid setup
        self.nz = nz if nz is not None else NZ
        self.ny = ny if ny is not None else NY
        self.z_max = z_max if z_max is not None else Z_MAX
        self.y_min = y_min if y_min is not None else Y_MIN
        self.y_max = y_max if y_max is not None else Y_MAX
        self.z = np.linspace(0, self.z_max, self.nz)
        self.y = np.linspace(self.y_min, self.y_max, self.ny)
        
        # Grid spacing (recalculated to match linspace exactly)
        self.dz = self.z[1] - self.z[0] if len(self.z) > 1 else DZ
        self.dy_deg = self.y[1] - self.y[0] if len(self.y) > 1 else DY
        self.dy = self.dy_deg * 111000.0  # Convert to meters
        
        # Initialize velocity fields
        self.u = np.zeros((self.nz, self.ny))
        self.v = np.zeros((self.nz, self.ny))
        self.w = np.zeros((self.nz, self.ny))
        
        # Time tracking
        self.time = 0.0
        self.timestep = 0
        
        # Initialize with Ekman spiral
        self._initialize_ekman()
        
        print(f"[{self.experiment_name}] Solver initialized")
        print(f"  Grid: nz={self.nz}, ny={self.ny}")
        print(f"  Domain: y=[{self.y[0]:.3f}°, {self.y[-1]:.3f}°], z=[0, {self.z[-1]:.0f}]m")
        print(f"  Pressure gradients: dpdy={self.dpdy:.2e}, dpdx={self.dpdx:.2e} m/s²")
    
    def _initialize_ekman(self):
        """Initialize velocity field with analytical Ekman solution"""
        D_e = EKMAN_DEPTH
        theta_e = EKMAN_SPIRAL_ANGLE

        # Precompute y-dependent Coriolis f (with paper equator rule)
        y_meters = self.y * 111000.0
        f = self.beta * y_meters
        # apply paper rule: if exact equator exists, set its f from +4 deg reference
        eq_indices = np.where(np.isclose(self.y, 0.0))[0]
        if eq_indices.size > 0:
            eq_idx = eq_indices[0]
            f[eq_idx] = self.beta * (EQUATOR_REFERENCE_LAT * 111000.0)

        # Loop over vertical levels only (nz is small); operations over ny are vectorized
        for iz, z_val in enumerate(self.z):
            if z_val <= D_e:
                zeta = z_val / D_e
                decay = np.exp(-np.pi * zeta)
                angle = theta_e * (1 - decay)
            else:
                decay = 0.0
                angle = theta_e

            # geostrophic wind at each latitude (vectorized)
            # avoid divide-by-zero by leaving f as set above
            u_geo = -self.dpdy / f
            v_geo = self.dpdx / f
            geostrophic_speed = np.hypot(u_geo, v_geo)

            # where geostrophic_speed == 0, result stays 0
            mask = geostrophic_speed > 0.0

            e_g_u = np.zeros_like(u_geo)
            e_g_v = np.zeros_like(v_geo)
            e_g_u[mask] = u_geo[mask] / geostrophic_speed[mask]
            e_g_v[mask] = v_geo[mask] / geostrophic_speed[mask]

            e_p_u = e_g_v
            e_p_v = -e_g_u

            along = geostrophic_speed * (1.0 - decay * np.cos(angle))
            cross = geostrophic_speed * decay * np.sin(angle)

            self.u[iz, :] = along * e_g_u + cross * e_p_u
            self.v[iz, :] = along * e_g_v + cross * e_p_v
    
    def _calculate_vertical_velocity(self):
        """Calculate w from continuity equation: ∂v/∂y + ∂w/∂z = 0"""
        # Vectorized computation of dv/dy for all levels
        dvdy = np.zeros_like(self.v)
        # interior centered differences
        dvdy[:, 1:-1] = (self.v[:, 2:] - self.v[:, :-2]) / (2 * self.dy)
        # boundaries one-sided
        dvdy[:, 0] = (self.v[:, 1] - self.v[:, 0]) / self.dy
        dvdy[:, -1] = (self.v[:, -1] - self.v[:, -2]) / self.dy

        # integrate vertically: w[0,:]=0 and w[iz+1] = -sum(dvdy[0:iz+1])*dz
        self.w = np.zeros_like(self.v)
        if self.nz > 1:
            self.w[1:, :] = -np.cumsum(dvdy[:-1, :] * self.dz, axis=0)
    
    def _advection_u(self):
        """Calculate advection terms for u: -v*∂u/∂y - w*∂u/∂z"""
        adv = np.zeros_like(self.u)

        # Y-direction upwind differences (vectorized)
        # left_diff = u - u shifted right (u[:,i]-u[:,i-1]) with zero at i=0
        left = self.u - np.concatenate([self.u[:, :1], self.u[:, :-1]], axis=1)
        right = np.concatenate([self.u[:, 1:], self.u[:, -1:]], axis=1) - self.u
        dudy = np.where(self.v > 0, left / self.dy, np.where(self.v < 0, right / self.dy, 0.0))

        # Z-direction upwind differences (vectorized)
        up = self.u - np.concatenate([self.u[:1, :], self.u[:-1, :]], axis=0)
        down = np.concatenate([self.u[1:, :], self.u[-1:, :]], axis=0) - self.u
        dudz = np.where(self.w > 0, up / self.dz, np.where(self.w < 0, down / self.dz, 0.0))

        adv = -self.v * dudy - self.w * dudz
        return adv
    
    def _advection_v(self):
        """Calculate advection terms for v: -v*∂v/∂y - w*∂v/∂z"""
        # vectorized analogous to _advection_u
        left = self.v - np.concatenate([self.v[:, :1], self.v[:, :-1]], axis=1)
        right = np.concatenate([self.v[:, 1:], self.v[:, -1:]], axis=1) - self.v
        dvdy = np.where(self.v > 0, left / self.dy, np.where(self.v < 0, right / self.dy, 0.0))

        up = self.v - np.concatenate([self.v[:1, :], self.v[:-1, :]], axis=0)
        down = np.concatenate([self.v[1:, :], self.v[-1:, :]], axis=0) - self.v
        dvdz = np.where(self.w > 0, up / self.dz, np.where(self.w < 0, down / self.dz, 0.0))

        adv = -self.v * dvdy - self.w * dvdz
        return adv
    
    def _diffusion_u(self):
        """Calculate vertical diffusion: K*∂²u/∂z²"""
        diff = np.zeros_like(self.u)
        if self.nz > 2:
            d2 = (self.u[2:, :] - 2 * self.u[1:-1, :] + self.u[:-2, :]) / (self.dz ** 2)
            diff[1:-1, :] = self.K * d2
        return diff
    
    def _diffusion_v(self):
        """Calculate vertical diffusion: K*∂²v/∂z²"""
        diff = np.zeros_like(self.v)
        if self.nz > 2:
            d2 = (self.v[2:, :] - 2 * self.v[1:-1, :] + self.v[:-2, :]) / (self.dz ** 2)
            diff[1:-1, :] = self.K * d2
        return diff
    
    def _coriolis_beta_term(self):
        """Calculate Coriolis terms: β*y*v for u, -β*y*u for v"""
        y_meters = self.y * 111000.0
        # broadcasting over vertical levels
        coriolis_u = self.beta * self.v * y_meters[np.newaxis, :]
        coriolis_v = -self.beta * self.u * y_meters[np.newaxis, :]
        return coriolis_u, coriolis_v
    
    def step(self, dt):
        """Perform one time step"""
        # Calculate vertical velocity from continuity
        self._calculate_vertical_velocity()
        
        # Calculate all terms
        adv_u = self._advection_u()
        adv_v = self._advection_v()
        diff_u = self._diffusion_u()
        diff_v = self._diffusion_v()
        cor_u, cor_v = self._coriolis_beta_term()
        
        # Update velocities using forward differencing
        u_new = self.u + dt * (adv_u + cor_u + diff_u - self.dpdx)
        v_new = self.v + dt * (adv_v + cor_v + diff_v - self.dpdy)
        
        # Apply boundary conditions
        # No-slip at surface
        u_new[0, :] = 0.0
        v_new[0, :] = 0.0
        
        # Zero vertical shear at top
        u_new[-1, :] = u_new[-2, :]
        v_new[-1, :] = v_new[-2, :]
        
        # Zero horizontal shear at lateral boundaries
        u_new[:, 0] = u_new[:, 1]
        u_new[:, -1] = u_new[:, -2]
        v_new[:, 0] = v_new[:, 1]
        v_new[:, -1] = v_new[:, -2]
        
        self.u = u_new
        self.v = v_new
        
        self.timestep += 1
        self.time += dt
    
    def integrate(self, dt, total_steps, save_steps, verbose=True):
        """
        Integrate forward in time
        
        Parameters:
        -----------
        dt : float
            Timestep (seconds)
        total_steps : int
            Total number of steps
        save_steps : int
            Save every N steps
        verbose : bool
            Print progress
        """
        
        # Storage for saved snapshots
        saved_steps = []
        u_history = []
        v_history = []
        w_history = []
        times = []
        
        if verbose:
            print(f"\n[{self.experiment_name}] Starting integration...")
            print(f"  Duration: {total_steps * dt / 86400:.1f} days ({total_steps} steps)")
            print(f"  Save interval: every {save_steps} steps")
        
        for step in range(total_steps):
            self.step(dt)
            
            # Save if needed
            if (step + 1) % save_steps == 0:
                saved_steps.append(step + 1)
                u_history.append(self.u.copy())
                v_history.append(self.v.copy())
                w_history.append(self.w.copy())
                times.append(self.time)
                
                if verbose:
                    elapsed_days = self.time / 86400
                    print(f"  Step {step + 1:5d}/{total_steps} ({elapsed_days:6.2f} days) - SAVED")
        
        if verbose:
            print(f"[{self.experiment_name}] Integration complete!")
            print(f"  Total saved snapshots: {len(saved_steps)}")
        
        return {
            'u_history': u_history,
            'v_history': v_history,
            'w_history': w_history,
            'times': times,
            'steps': saved_steps,
        }
    
    def get_current_state(self):
        """Get current velocity fields"""
        return {
            'u': self.u.copy(),
            'v': self.v.copy(),
            'w': self.w.copy(),
            'z': self.z.copy(),
            'y': self.y.copy(),
            'time': self.time,
            'timestep': self.timestep,
        }
