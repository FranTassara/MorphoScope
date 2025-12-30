#!/usr/bin/env python3
"""
Synthetic Volume Generator for MorphoScope
==========================================

Generates synthetic 3D volumes with known geometries to:
1. Validate the spread quantification algorithm
2. Help biologists understand what spread measures
3. Demonstrate how different variables affect spread values

Volume specifications:
- Z: 20 slices, 1 µm per slice (20 µm total depth)
- XY: 800x800 pixels, 0.1 µm per pixel (80x80 µm field of view)

Author: Francisco Tassara (MorphoScope project - Ceriani Lab)
"""

import numpy as np
from pathlib import Path
import tifffile
from typing import Tuple, Optional
from dataclasses import dataclass


@dataclass
class VolumeSpecs:
    """Volume specifications matching typical confocal imaging parameters."""
    # Dimensions in pixels/slices
    nx: int = 800
    ny: int = 800
    nz: int = 20
    
    # Physical sizes in micrometers
    pixel_size_um: float = 0.1  # µm per pixel in XY
    z_step_um: float = 1.0      # µm per slice in Z
    
    @property
    def size_x_um(self) -> float:
        return self.nx * self.pixel_size_um
    
    @property
    def size_y_um(self) -> float:
        return self.ny * self.pixel_size_um
    
    @property
    def size_z_um(self) -> float:
        return self.nz * self.z_step_um
    
    def um_to_pixels_xy(self, um: float) -> float:
        """Convert micrometers to pixels in XY."""
        return um / self.pixel_size_um
    
    def um_to_slices_z(self, um: float) -> float:
        """Convert micrometers to slices in Z."""
        return um / self.z_step_um


class SyntheticVolumeGenerator:
    """
    Generate synthetic 3D volumes with known geometries.
    
    All positions and sizes can be specified in micrometers.
    The generator handles conversion to pixel coordinates internally.
    """
    
    def __init__(self, specs: Optional[VolumeSpecs] = None):
        self.specs = specs or VolumeSpecs()
        
        # Create coordinate grids in micrometers
        x_um = np.arange(self.specs.nx) * self.specs.pixel_size_um
        y_um = np.arange(self.specs.ny) * self.specs.pixel_size_um
        z_um = np.arange(self.specs.nz) * self.specs.z_step_um
        
        # Meshgrid: Z, Y, X order for proper array indexing
        self.Z, self.Y, self.X = np.meshgrid(z_um, y_um, x_um, indexing='ij')
    
    def _create_empty_volume(self, dtype=np.float32) -> np.ndarray:
        """Create an empty volume array."""
        return np.zeros((self.specs.nz, self.specs.ny, self.specs.nx), dtype=dtype)
    
    def _normalize_and_convert(self, volume: np.ndarray, max_intensity: int = 255) -> np.ndarray:
        """Normalize volume and convert to uint8."""
        if volume.max() > 0:
            volume = volume / volume.max() * max_intensity
        return volume.astype(np.uint8)
    
    # =========================================================================
    # Basic geometric primitives
    # =========================================================================
    
    def create_sphere(
        self,
        center_um: Tuple[float, float, float],
        radius_um: float,
        intensity: float = 1.0,
        smooth_edge: bool = True,
        edge_width_um: float = 1.0
    ) -> np.ndarray:
        """
        Create a sphere.
        
        Parameters
        ----------
        center_um : tuple (x, y, z) in micrometers
        radius_um : radius in micrometers
        intensity : peak intensity (0-1)
        smooth_edge : if True, apply Gaussian falloff at edge
        edge_width_um : width of smooth edge transition
        
        Returns
        -------
        volume : 3D numpy array
        """
        cx, cy, cz = center_um
        
        # Distance from center
        dist = np.sqrt((self.X - cx)**2 + (self.Y - cy)**2 + (self.Z - cz)**2)
        
        if smooth_edge:
            # Smooth edge using sigmoid-like function
            volume = intensity * np.exp(-((dist - radius_um).clip(min=0) / edge_width_um)**2)
        else:
            volume = np.where(dist <= radius_um, intensity, 0.0)
        
        return volume.astype(np.float32)
    
    def create_ellipsoid(
        self,
        center_um: Tuple[float, float, float],
        radii_um: Tuple[float, float, float],
        intensity: float = 1.0,
        smooth_edge: bool = True,
        edge_width_um: float = 1.0,
        rotation_deg: float = 0.0
    ) -> np.ndarray:
        """
        Create an ellipsoid, optionally rotated in the XY plane.
        
        Parameters
        ----------
        center_um : tuple (x, y, z) in micrometers
        radii_um : tuple (rx, ry, rz) semi-axes in micrometers
        rotation_deg : rotation angle in XY plane (degrees)
        
        Returns
        -------
        volume : 3D numpy array
        """
        cx, cy, cz = center_um
        rx, ry, rz = radii_um
        
        # Translate to center
        dx = self.X - cx
        dy = self.Y - cy
        dz = self.Z - cz
        
        # Apply rotation in XY plane
        if rotation_deg != 0:
            theta = np.radians(rotation_deg)
            dx_rot = dx * np.cos(theta) + dy * np.sin(theta)
            dy_rot = -dx * np.sin(theta) + dy * np.cos(theta)
            dx, dy = dx_rot, dy_rot
        
        # Normalized distance (1.0 at surface)
        norm_dist = np.sqrt((dx/rx)**2 + (dy/ry)**2 + (dz/rz)**2)
        
        if smooth_edge:
            volume = intensity * np.exp(-((norm_dist - 1.0).clip(min=0) * min(rx, ry, rz) / edge_width_um)**2)
        else:
            volume = np.where(norm_dist <= 1.0, intensity, 0.0)
        
        return volume.astype(np.float32)
    
    def create_cylinder(
        self,
        start_um: Tuple[float, float, float],
        end_um: Tuple[float, float, float],
        radius_um: float,
        intensity: float = 1.0,
        smooth_edge: bool = True,
        edge_width_um: float = 0.5
    ) -> np.ndarray:
        """
        Create a cylinder between two points.
        
        Parameters
        ----------
        start_um, end_um : endpoints (x, y, z) in micrometers
        radius_um : cylinder radius in micrometers
        
        Returns
        -------
        volume : 3D numpy array
        """
        p1 = np.array(start_um)
        p2 = np.array(end_um)
        
        # Cylinder axis
        axis = p2 - p1
        length = np.linalg.norm(axis)
        axis_norm = axis / length
        
        # Vector from p1 to each point
        points = np.stack([self.X, self.Y, self.Z], axis=-1)
        v = points - p1
        
        # Project onto axis
        t = np.dot(v, axis_norm)
        
        # Perpendicular distance from axis
        proj = np.outer(t.ravel(), axis_norm).reshape(points.shape)
        perp = v - proj
        dist_from_axis = np.linalg.norm(perp, axis=-1)
        
        # Inside cylinder: 0 <= t <= length and dist <= radius
        inside_length = (t >= 0) & (t <= length)
        
        if smooth_edge:
            radial_falloff = np.exp(-((dist_from_axis - radius_um).clip(min=0) / edge_width_um)**2)
            volume = intensity * radial_falloff * inside_length
        else:
            inside_radius = dist_from_axis <= radius_um
            volume = np.where(inside_length & inside_radius, intensity, 0.0)
        
        return volume.astype(np.float32)
    
    def create_gaussian_blob(
        self,
        center_um: Tuple[float, float, float],
        sigma_um: Tuple[float, float, float],
        intensity: float = 1.0,
        rotation_deg: float = 0.0
    ) -> np.ndarray:
        """
        Create a 3D Gaussian blob (more realistic for fluorescence).
        
        Parameters
        ----------
        center_um : tuple (x, y, z) in micrometers
        sigma_um : tuple (σx, σy, σz) standard deviations in micrometers
        rotation_deg : rotation angle in XY plane (degrees)
        
        Returns
        -------
        volume : 3D numpy array
        """
        cx, cy, cz = center_um
        sx, sy, sz = sigma_um
        
        dx = self.X - cx
        dy = self.Y - cy
        dz = self.Z - cz
        
        if rotation_deg != 0:
            theta = np.radians(rotation_deg)
            dx_rot = dx * np.cos(theta) + dy * np.sin(theta)
            dy_rot = -dx * np.sin(theta) + dy * np.cos(theta)
            dx, dy = dx_rot, dy_rot
        
        volume = intensity * np.exp(-0.5 * ((dx/sx)**2 + (dy/sy)**2 + (dz/sz)**2))
        
        return volume.astype(np.float32)
    
    # =========================================================================
    # Neuronal-like structures
    # =========================================================================
    
    def create_branching_structure(
        self,
        root_um: Tuple[float, float, float],
        main_length_um: float,
        main_direction: Tuple[float, float, float] = (1, 0, 0),
        branch_length_um: float = 10.0,
        n_branches: int = 5,
        branch_spread_um: float = 5.0,
        fiber_radius_um: float = 1.0,
        intensity: float = 1.0
    ) -> np.ndarray:
        """
        Create a branching axon-like structure.
        
        Parameters
        ----------
        root_um : starting point (x, y, z)
        main_length_um : length of main axon
        main_direction : direction vector of main axon
        branch_length_um : typical length of branches
        n_branches : number of branches
        branch_spread_um : spread of branch endpoints from main axis
        fiber_radius_um : radius of all fibers
        
        Returns
        -------
        volume : 3D numpy array
        """
        volume = self._create_empty_volume()
        
        # Normalize direction
        d = np.array(main_direction, dtype=float)
        d = d / np.linalg.norm(d)
        
        root = np.array(root_um)
        tip = root + d * main_length_um
        
        # Main axon
        volume += self.create_cylinder(root, tip, fiber_radius_um, intensity)
        
        # Create branches at regular intervals along main axon
        for i in range(n_branches):
            # Position along main axon
            t = (i + 1) / (n_branches + 1)
            branch_start = root + d * main_length_um * t
            
            # Random perpendicular direction for branch
            np.random.seed(i * 42)  # Reproducible
            perp = np.array([d[1], -d[0], 0])  # Perpendicular in XY
            if np.linalg.norm(perp) < 0.1:
                perp = np.array([0, d[2], -d[1]])
            perp = perp / np.linalg.norm(perp)
            
            # Add some Z variation
            z_offset = (np.random.rand() - 0.5) * 2
            branch_dir = perp + np.array([0, 0, z_offset * 0.3])
            branch_dir = branch_dir / np.linalg.norm(branch_dir)
            
            # Alternate sides
            if i % 2 == 0:
                branch_dir = -branch_dir
            
            branch_end = branch_start + branch_dir * branch_length_um
            
            # Clip to volume bounds
            branch_end = np.clip(branch_end, [1, 1, 1], 
                                 [self.specs.size_x_um-1, self.specs.size_y_um-1, self.specs.size_z_um-1])
            
            volume += self.create_cylinder(branch_start, branch_end, fiber_radius_um * 0.7, intensity * 0.9)
        
        return np.clip(volume, 0, intensity).astype(np.float32)
    
    def create_fasciculated_bundle(
        self,
        center_um: Tuple[float, float, float],
        length_um: float,
        direction: Tuple[float, float, float] = (1, 0, 0),
        n_fibers: int = 4,
        bundle_radius_um: float = 2.0,
        fiber_radius_um: float = 0.8,
        intensity: float = 1.0
    ) -> np.ndarray:
        """
        Create a fasciculated (bundled/compact) structure.
        
        This simulates the "dusk" state of s-LNv projections where
        axons are tightly bundled together.
        
        Parameters
        ----------
        center_um : center of the bundle
        length_um : length of the bundle
        n_fibers : number of fibers in bundle
        bundle_radius_um : how tightly packed the fibers are
        fiber_radius_um : radius of individual fibers
        
        Returns
        -------
        volume : 3D numpy array
        """
        volume = self._create_empty_volume()
        
        d = np.array(direction, dtype=float)
        d = d / np.linalg.norm(d)
        
        center = np.array(center_um)
        
        # Create fibers arranged in a tight bundle
        for i in range(n_fibers):
            angle = 2 * np.pi * i / n_fibers
            
            # Perpendicular offset
            if abs(d[2]) < 0.9:
                perp1 = np.cross(d, [0, 0, 1])
            else:
                perp1 = np.cross(d, [1, 0, 0])
            perp1 = perp1 / np.linalg.norm(perp1)
            perp2 = np.cross(d, perp1)
            
            offset = (perp1 * np.cos(angle) + perp2 * np.sin(angle)) * bundle_radius_um
            
            start = center - d * length_um / 2 + offset
            end = center + d * length_um / 2 + offset
            
            volume += self.create_cylinder(start, end, fiber_radius_um, intensity)
        
        return np.clip(volume, 0, intensity).astype(np.float32)
    
    def create_defasciculated_structure(
        self,
        center_um: Tuple[float, float, float],
        base_length_um: float,
        direction: Tuple[float, float, float] = (1, 0, 0),
        n_fibers: int = 4,
        spread_angle_deg: float = 30.0,
        fiber_radius_um: float = 0.8,
        intensity: float = 1.0
    ) -> np.ndarray:
        """
        Create a defasciculated (spread out) structure.
        
        This simulates the "dawn" state of s-LNv projections where
        axons spread out from each other.
        
        Parameters
        ----------
        center_um : starting point where fibers diverge
        base_length_um : length of each fiber
        n_fibers : number of fibers
        spread_angle_deg : angle of spread from main axis
        fiber_radius_um : radius of fibers
        
        Returns
        -------
        volume : 3D numpy array
        """
        volume = self._create_empty_volume()
        
        d = np.array(direction, dtype=float)
        d = d / np.linalg.norm(d)
        
        center = np.array(center_um)
        spread_rad = np.radians(spread_angle_deg)
        
        # Create perpendicular basis
        if abs(d[2]) < 0.9:
            perp1 = np.cross(d, [0, 0, 1])
        else:
            perp1 = np.cross(d, [1, 0, 0])
        perp1 = perp1 / np.linalg.norm(perp1)
        perp2 = np.cross(d, perp1)
        
        for i in range(n_fibers):
            angle = 2 * np.pi * i / n_fibers
            
            # Direction with spread
            spread_dir = (d * np.cos(spread_rad) + 
                         (perp1 * np.cos(angle) + perp2 * np.sin(angle)) * np.sin(spread_rad))
            spread_dir = spread_dir / np.linalg.norm(spread_dir)
            
            end = center + spread_dir * base_length_um
            
            volume += self.create_cylinder(center, end, fiber_radius_um, intensity)
        
        return np.clip(volume, 0, intensity).astype(np.float32)
    
    # =========================================================================
    # Pedagogical demonstration volumes
    # =========================================================================
    
    def demo_spread_concept(self) -> dict:
        """
        Generate volumes that demonstrate what spread measures.
        
        Returns a dictionary of named volumes.
        """
        volumes = {}
        center = (40, 40, 10)  # Center of volume in µm
        
        # 1. Sphere (equal spread in all axes)
        volumes['01_sphere_r10um'] = self.create_sphere(center, radius_um=10.0)
        
        # 2. Sphere with different radius (spread scales with size)
        volumes['02_sphere_r5um'] = self.create_sphere(center, radius_um=5.0)
        volumes['03_sphere_r15um'] = self.create_sphere(center, radius_um=15.0)
        
        # 3. Ellipsoid elongated in X (high X spread, low Y/Z spread)
        volumes['04_ellipsoid_x_elongated'] = self.create_ellipsoid(
            center, radii_um=(20, 5, 5)
        )
        
        # 4. Ellipsoid elongated in Y
        volumes['05_ellipsoid_y_elongated'] = self.create_ellipsoid(
            center, radii_um=(5, 20, 5)
        )
        
        # 5. Ellipsoid elongated in Z
        volumes['06_ellipsoid_z_elongated'] = self.create_ellipsoid(
            center, radii_um=(5, 5, 8)  # Limited by 20 slices
        )
        
        # 6. Same geometry, different intensity (spread should be same)
        volumes['07_sphere_low_intensity'] = self.create_sphere(
            center, radius_um=10.0, intensity=0.3
        )
        volumes['08_sphere_high_intensity'] = self.create_sphere(
            center, radius_um=10.0, intensity=1.0
        )
        
        return volumes
    
    def demo_pca_rotation(self) -> dict:
        """
        Generate volumes that demonstrate PCA axis alignment.
        
        The same ellipsoid at different rotations should give
        the same spread values after PCA alignment.
        """
        volumes = {}
        center = (40, 40, 10)
        radii = (20, 8, 5)
        
        for angle in [0, 30, 45, 60, 90]:
            volumes[f'ellipsoid_rot{angle:03d}deg'] = self.create_ellipsoid(
                center, radii_um=radii, rotation_deg=angle
            )
        
        return volumes
    
    def demo_fasciculation(self) -> dict:
        """
        Generate volumes demonstrating fasciculation vs defasciculation.
        
        This is the key biological phenomenon: s-LNv axons are
        - Fasciculated (bundled) at dusk → low Y/Z spread
        - Defasciculated (spread) at dawn → high Y/Z spread
        """
        volumes = {}
        center = (40, 40, 10)
        
        # Fasciculated (dusk-like, compact bundle)
        volumes['01_fasciculated_tight'] = self.create_fasciculated_bundle(
            center, length_um=50, n_fibers=5, bundle_radius_um=2.0
        )
        
        # Partially defasciculated
        volumes['02_partially_spread'] = self.create_defasciculated_structure(
            center, base_length_um=30, n_fibers=5, spread_angle_deg=15
        )
        
        # Fully defasciculated (dawn-like, spread out)
        volumes['03_defasciculated_spread'] = self.create_defasciculated_structure(
            center, base_length_um=30, n_fibers=5, spread_angle_deg=35
        )
        
        # Branching structure (more realistic)
        volumes['04_branching_compact'] = self.create_branching_structure(
            root_um=(10, 40, 10), main_length_um=50, 
            n_branches=4, branch_spread_um=3, branch_length_um=8
        )
        
        volumes['05_branching_spread'] = self.create_branching_structure(
            root_um=(10, 40, 10), main_length_um=50,
            n_branches=6, branch_spread_um=8, branch_length_um=15
        )
        
        return volumes
    
    def demo_axonal_volume(self) -> dict:
        """
        Generate volumes demonstrating axonal volume concept.
        
        Axonal volume = sum of all intensity values × voxel size
        This represents total fluorophore content.
        """
        volumes = {}
        center = (40, 40, 10)
        
        # Same size, different "density" (intensity)
        volumes['01_dense_structure'] = self.create_sphere(
            center, radius_um=10, intensity=1.0
        )
        volumes['02_sparse_structure'] = self.create_sphere(
            center, radius_um=10, intensity=0.5
        )
        
        # Different sizes (different total material)
        volumes['03_small_full'] = self.create_sphere(
            center, radius_um=7, intensity=1.0
        )
        volumes['04_large_full'] = self.create_sphere(
            center, radius_um=12, intensity=1.0
        )
        
        # Same spread, different volume (hollow vs solid concept)
        # Thick shell
        shell = self.create_sphere(center, radius_um=12, intensity=1.0)
        inner = self.create_sphere(center, radius_um=8, intensity=1.0)
        volumes['05_thick_shell'] = np.clip(shell - inner * 0.8, 0, 1).astype(np.float32)
        
        # Solid
        volumes['06_solid_sphere'] = self.create_sphere(
            center, radius_um=12, intensity=1.0
        )
        
        return volumes
    
    def generate_all_demos(self, output_dir: str = "synthetic_volumes"):
        """
        Generate all demonstration volumes and save as TIFF stacks.
        
        Parameters
        ----------
        output_dir : output directory path
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Generate all demo sets
        demo_sets = {
            'spread_concept': self.demo_spread_concept(),
            'pca_rotation': self.demo_pca_rotation(),
            'fasciculation': self.demo_fasciculation(),
            'axonal_volume': self.demo_axonal_volume(),
        }
        
        # Save each volume
        for category, volumes in demo_sets.items():
            category_dir = output_path / category
            category_dir.mkdir(exist_ok=True)
            
            for name, volume in volumes.items():
                # Convert to uint8 for saving
                vol_uint8 = self._normalize_and_convert(volume)
                
                # Save as TIFF stack
                filepath = category_dir / f"{name}.tif"
                tifffile.imwrite(filepath, vol_uint8, imagej=True,
                                metadata={'spacing': self.specs.z_step_um,
                                         'unit': 'um'})
                
                print(f"Saved: {filepath}")
        
        # Save metadata file
        self._save_metadata(output_path)
        
        print(f"\n✓ All synthetic volumes saved to: {output_path.absolute()}")
        print(f"\nVolume specifications:")
        print(f"  - XY: {self.specs.nx}×{self.specs.ny} pixels ({self.specs.size_x_um}×{self.specs.size_y_um} µm)")
        print(f"  - Z: {self.specs.nz} slices ({self.specs.size_z_um} µm)")
        print(f"  - Pixel size: {self.specs.pixel_size_um} µm")
        print(f"  - Z step: {self.specs.z_step_um} µm")
    
    def _save_metadata(self, output_path: Path):
        """Save metadata file describing all generated volumes."""
        metadata = f"""
Synthetic Volumes for MorphoScope Validation
=============================================

Generated volumes to demonstrate and validate the spread quantification algorithm.

Volume Specifications
---------------------
- XY dimensions: {self.specs.nx} x {self.specs.ny} pixels
- Z dimensions: {self.specs.nz} slices
- Pixel size (XY): {self.specs.pixel_size_um} um
- Z step: {self.specs.z_step_um} um
- Physical size: {self.specs.size_x_um} x {self.specs.size_y_um} x {self.specs.size_z_um} um

Demo Categories
---------------

1. spread_concept/
   Demonstrates what spread measures:
   - Spheres of different sizes: spread scales with size
   - Ellipsoids elongated in different axes: anisotropic spread
   - Same geometry at different intensities: spread unchanged
   
2. pca_rotation/
   Demonstrates PCA axis alignment:
   - Same ellipsoid at different rotations
   - After PCA alignment, all should give similar spread values
   
3. fasciculation/
   Demonstrates the biological phenomenon:
   - fasciculated (dusk state): axons bundled tightly, low Y/Z spread
   - defasciculated (dawn state): axons spread out, high Y/Z spread
   
4. axonal_volume/
   Demonstrates axonal volume measurement:
   - Same size, different intensity: different volume
   - Different sizes: different volume
   - Shell vs solid: spread similar, volume different

Usage with MorphoScope
----------------------
1. Load any .tif file from these directories
2. The entire image can be used as ROI (no selection needed)
3. Compare spread values between related volumes
4. Use to verify algorithm correctness against known geometries

Expected Results
----------------
spread_concept/:
  - Larger spheres: larger spread values
  - X-elongated ellipsoid: largest X spread
  - Different intensity: same spread (intensity-weighted mean)

pca_rotation/:
  - All rotations should give similar spread after PCA alignment
  - X spread should always be largest (PCA aligns to maximum variance)

fasciculation/:
  - Fasciculated: low Y/Z spread, similar X spread
  - Defasciculated: high Y/Z spread
  - Branching structures: intermediate values

axonal_volume/:
  - Higher intensity: higher axonal volume
  - Larger structure: higher axonal volume
  - Shell vs solid: similar spread, different volume
"""
        
        with open(output_path / "README.txt", 'w', encoding='utf-8') as f:
            f.write(metadata)


def main():
    """Generate all synthetic demonstration volumes."""
    print("=" * 60)
    print("Synthetic Volume Generator for MorphoScope")
    print("=" * 60)
    print()
    
    # Create generator with default specs
    generator = SyntheticVolumeGenerator()
    
    # Generate all demos
    generator.generate_all_demos("W:/Fran/Script Plasticidad/GUI/Test/synthetic_volumes")
    
    print("\n" + "=" * 60)
    print("Generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
