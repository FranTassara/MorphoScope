"""
Core Image Processor for Structural Plasticity Analysis

Implements the quantification algorithm described in Petsakou et al., 2015:
"Circadian rhythms in Rho1 activity regulate neuronal plasticity and network hierarchy"

This module provides tools to analyze 3D microscopy images of neuronal projections,
quantifying their structural plasticity by measuring axonal spread in three dimensions.

Key Features:
- PCA-based image alignment for consistent orientation
- Calculation of local and global spread metrics
- Comprehensive fluorescence quantification

Author: Francisco Tassara
Date: 2025-11-12
Based on: Petsakou, Sapsis & Blau, Cell 2015
"""

import numpy as np
from scipy import ndimage as ndi
from typing import Tuple, Dict, Optional
import logging
import matplotlib.pyplot as plt
from pathlib import Path

logger = logging.getLogger(__name__)


class ImageProcessor:
    """
    Core processor for structural plasticity quantification following Petsakou et al., 2015.
    
    Coordinate System Convention
    ----------------------------
    After PCA rotation and standardization (90° rotation):
    - X axis: HORIZONTAL direction (maximum projection length)
    - Y axis: VERTICAL direction (perpendicular to X, typically narrower)
    - Z axis: Depth direction (confocal slices, dorsal to ventral)
    
    Image Format
    ------------
    All methods expect images in (Y, X, Z) format = (rows, cols, slices)
    This matches numpy's natural indexing: image[y, x, z]
    
    Notes
    -----
    - The algorithm applies PCA to find the principal axis of the projection
    - An additional 90° rotation standardizes the coordinate system
    - Spread metrics are calculated along the mean curve of the projection
    """
    
    def __init__(self, voxel_size_x: float, voxel_size_y: float, voxel_size_z: float):
        """
        Initialize processor with voxel dimensions.
        
        Parameters
        ----------
        voxel_size_x : float
            Voxel size in X dimension (µm) - horizontal after standardization
        voxel_size_y : float
            Voxel size in Y dimension (µm) - vertical after standardization
        voxel_size_z : float
            Voxel size in Z dimension (µm/slice) - depth direction
        
        Examples
        --------
        >>> processor = ImageProcessor(voxel_size_x=0.124, voxel_size_y=0.124, voxel_size_z=1.0)
        """
        self.voxel_size_x = voxel_size_x
        self.voxel_size_y = voxel_size_y
        self.voxel_size_z = voxel_size_z
        
        logger.info(f"ImageProcessor initialized - Voxel sizes: "
                   f"X={voxel_size_x:.4f}µm, Y={voxel_size_y:.4f}µm, Z={voxel_size_z:.4f}µm")
    
    
    def calculate_fluorescence(
        self, 
        image_3d: np.ndarray, 
        mask_area_pixels: float
    ) -> Tuple[float, float]:
        """
        Calculate normalized fluorescence metrics.
        
        Fluorescence is computed as total intensity normalized by the ROI area,
        providing both pixel-based and physical (µm²) measurements.
        
        Parameters
        ----------
        image_3d : np.ndarray
            3D image stack (any shape compatible with processing)
        mask_area_pixels : float
            Region of interest (ROI) area in pixels (2D projection area)
        
        Returns
        -------
        Tuple[float, float]
            - fluorescence_per_pixel: Total intensity / ROI area (pixels)
            - fluorescence_per_um2: Total intensity / ROI area (µm²)
        
        Notes
        -----
        The fluorescence metric accounts for differences in imaging parameters
        and allows comparison across different experimental conditions.
        """
        total_intensity = np.sum(image_3d)
        
        # Fluorescence normalized by pixel count
        fluor_px = total_intensity / mask_area_pixels
        
        # Fluorescence normalized by physical area (µm²)
        area_um2 = mask_area_pixels * (self.voxel_size_x * self.voxel_size_y)
        fluor_um = total_intensity / area_um2
        
        logger.debug(f"Fluorescence - Per pixel: {fluor_px:.2f}, Per µm²: {fluor_um:.2f}")
        
        return round(fluor_px, 2), round(fluor_um, 2)
    
    
    def calculate_pca_rotation_angle(self, image_yxz: np.ndarray) -> float:
        """
        Calculate optimal rotation angle using Principal Component Analysis.
        
        This method finds the principal axis of the Z-projected image to align
        the projection along its maximum length direction. The first principal
        component (eigenvector with largest eigenvalue) defines the rotation angle.
        
        Parameters
        ----------
        image_yxz : np.ndarray
            3D image with shape (Y, X, Z) = (rows, cols, slices)
        
        Returns
        -------
        float
            Rotation angle in degrees. Positive values indicate counter-clockwise rotation.
        
        Notes
        -----
        - The Z-projection is intensity-weighted for better alignment
        - Covariance matrix is computed from the normalized projection
        - Returns 0.0 if projection is empty or invalid
        
        References
        ----------
        Lee, J.M. (2007). Introduction to Smooth Manifolds. Springer.
        """
        # Create Z-projection (sum along depth axis)
        projection_z = np.sum(image_yxz, axis=2)  # Shape: (Y, X) = (rows, cols)
        
        # Normalize projection for PCA
        sum_projection = np.sum(projection_z)
        if sum_projection == 0:
            logger.warning("Empty Z-projection - returning 0° rotation angle")
            return 0.0
        
        projection_z_norm = projection_z / sum_projection
        
        # Create coordinate grids
        nrows, ncols = projection_z.shape  # (Y, X)
        xx = np.arange(ncols)  # X coordinates (columns)
        yy = np.arange(nrows)  # Y coordinates (rows)
        
        # Create meshgrid with 'ij' indexing for proper row/column correspondence
        YY, XX = np.meshgrid(yy, xx, indexing='ij')  # Both shape (Y, X)
        
        # Calculate weighted centroid
        Mx = np.sum(XX * projection_z_norm)
        My = np.sum(YY * projection_z_norm)
        
        # Calculate covariance matrix components
        Cxx = np.sum(projection_z_norm * (XX - Mx)**2)
        Cyy = np.sum(projection_z_norm * (YY - My)**2)
        Cxy = np.sum(projection_z_norm * (XX - Mx) * (YY - My))
        
        # Construct covariance matrix
        cov_matrix = np.array([[Cxx, Cxy], 
                               [Cxy, Cyy]])
        
        # Eigenvalue decomposition
        eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
        
        # First eigenvector (largest eigenvalue) defines principal axis
        # Calculate rotation angle from eigenvector components
        phi = np.arctan(eigenvectors[1, 0] / eigenvectors[0, 0]) * 180 / np.pi
        
        logger.info(f"PCA rotation angle: {phi:.2f}°")
        logger.debug(f"Eigenvalues: {eigenvalues}")
        
        return phi
    
    
    def rotate_image(self, image_yxz: np.ndarray, angle: float) -> np.ndarray:
        """
        Rotate 3D image by specified angle in the XY plane.
        
        This method applies the rotation slice-by-slice to maintain the Z-axis alignment.
        Bilinear interpolation (order=1) is used for smooth rotation while preserving
        image properties.
        
        Parameters
        ----------
        image_yxz : np.ndarray
            3D image with shape (Y, X, Z) = (rows, cols, slices)
        angle : float
            Rotation angle in degrees (positive = counter-clockwise)
        
        Returns
        -------
        np.ndarray
            Rotated 3D image. Shape may change if reshape=True is used.
            Output is in float64 format to preserve intensity values.
        
        Notes
        -----
        - Each Z-slice is rotated independently
        - The image is automatically resized to fit the rotated content
        - Intensity is approximately conserved (small losses due to interpolation)
        
        See Also
        --------
        scipy.ndimage.rotate : Underlying rotation function
        """
        n_slices = image_yxz.shape[2]
        
        # Determine output shape from first slice rotation
        rotated_first = ndi.rotate(
            image_yxz[:, :, 0], 
            angle=angle, 
            reshape=True,  # Adjust canvas size to fit rotated content
            order=1        # Bilinear interpolation
        )
        
        # Preallocate output array
        rotated_image = np.zeros(
            (rotated_first.shape[0], rotated_first.shape[1], n_slices), 
            dtype=image_yxz.dtype
        )
        
        # Rotate each Z-slice
        logger.debug(f"Rotating {n_slices} slices by {angle:.2f}°")
        for zi in range(n_slices):
            rotated_image[:, :, zi] = ndi.rotate(
                image_yxz[:, :, zi], 
                angle=angle, 
                reshape=True,
                order=1
            )
        
        # Log rotation statistics
        intensity_loss = (1 - np.sum(rotated_image) / np.sum(image_yxz)) * 100
        logger.debug(f"Rotation complete - New shape: {rotated_image.shape}")
        logger.debug(f"Intensity loss: {intensity_loss:.2f}%")
        logger.debug(f"Non-zero voxels: {np.count_nonzero(rotated_image)}")
        
        return rotated_image
    
    
    def plot_z_projection(
        self,
        image_yxz: np.ndarray,
        title: str = "Z-Projection",
        save_path: Optional[str] = None
    ) -> None:
        """
        Create and display/save a Z-projection visualization.
        
        This method generates a maximum intensity projection along the Z-axis
        for visual inspection of the image orientation and quality.
        
        Parameters
        ----------
        image_yxz : np.ndarray
            3D image with shape (Y, X, Z) = (rows, cols, slices)
        title : str, optional
            Plot title (default: "Z-Projection")
        save_path : Optional[str], optional
            If provided, save figure to this path instead of displaying it
        
        Notes
        -----
        - Uses 'hot' colormap for better visualization of intensity variations
        - Includes intensity statistics in the plot
        - Grid overlay helps with spatial reference
        """
        # Create Z-projection (sum along depth axis)
        projection = np.sum(image_yxz, axis=2)
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Display projection with hot colormap
        im = ax.imshow(
            projection, 
            cmap='hot', 
            origin='lower',  # Origin at bottom-left
            interpolation='nearest'
        )
        
        # Set labels and title
        ax.set_title(
            f'{title}\nShape: {image_yxz.shape}', 
            fontsize=14, 
            fontweight='bold'
        )
        ax.set_xlabel('X axis (pixels) - HORIZONTAL', fontsize=12)
        ax.set_ylabel('Y axis (pixels) - VERTICAL', fontsize=12)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Integrated Intensity', fontsize=12)
        
        # Add grid for spatial reference
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        
        # Add statistics text box
        stats_text = (
            'Statistics:\n'
            f'Min: {projection.min():.1f}\n'
            f'Max: {projection.max():.1f}\n'
            f'Mean: {projection.mean():.1f}\n'
            f'Total: {np.sum(projection):.0f}'
        )
        
        ax.text(
            0.02, 0.98, stats_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        )
        
        plt.tight_layout()
        
        # Save or display
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            logger.info(f"Z-projection saved to: {save_path}")
            plt.close(fig)
        else:
            plt.show()
    
    
    def calculate_local_spreads(
        self, 
        image_yxz_rotated: np.ndarray,
        after_90deg_rotation: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calculate local spread metrics along the principal axis.
        
        This method computes intensity-weighted variance (spread) in the Y and Z
        directions for each position along the X axis. The spreads quantify how
        the projection is distributed around its mean curve.
        
        Parameters
        ----------
        image_yxz_rotated : np.ndarray
            Rotated 3D image with shape (Y, X, Z) = (rows, cols, slices)
        after_90deg_rotation : bool, optional
            If True, assumes the standard coordinate convention where X is horizontal (columns).
            If False, uses MATLAB convention where X corresponds to rows.
            Default: False
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
            - MMsum: Total intensity at each X position
            - MMyy: Local variance in Y direction at each X position
            - MMzz: Local variance in Z direction at each X position
            - xx00: X coordinate array (1-indexed)
        
        Notes
        -----
        - Variances are intensity-weighted to account for non-uniform distributions
        - NaN values indicate positions with no signal
        - The method adapts to different coordinate conventions via `after_90deg_rotation`
        
        Mathematical Details
        --------------------
        For each position x_i:
        - MMsum[i] = Σ intensity at x_i
        - MMyy[i] = Σ (y - mean_y)² × intensity / MMsum[i]  (weighted variance)
        - MMzz[i] = Σ (z - mean_z)² × intensity / MMsum[i]  (weighted variance)
        """
        nrows, ncols, depth = image_yxz_rotated.shape
        
        if after_90deg_rotation:
            # Standard convention: X=horizontal (columns), Y=vertical (rows)
            logger.debug("Processing with standard convention: X=horizontal(cols), Y=vertical(rows)")
            
            xx00 = np.arange(1, ncols + 1)  # X coordinates (columns), 1-indexed
            yy00 = np.arange(1, nrows + 1)  # Y coordinates (rows), 1-indexed
            zz00 = np.arange(1, depth + 1)  # Z coordinates (slices), 1-indexed
            
            # Preallocate arrays
            MMsum = np.zeros(ncols)
            MMyy = np.zeros(ncols)
            MMzz = np.zeros(ncols)
            
            # Iterate over columns (X axis)
            for i in range(ncols):
                # Extract all rows for column i: shape (nrows, depth) = (Y, Z)
                MM = image_yxz_rotated[:, i, :].T  # Transpose to (Z, Y)
                
                MMsum[i] = np.sum(MM)
                
                if MMsum[i] > 0:
                    # Calculate weighted means
                    MMz = np.sum(zz00 * np.sum(MM, axis=1)) / MMsum[i]
                    MMy = np.sum(yy00 * np.sum(MM, axis=0)) / MMsum[i]
                    
                    # Calculate weighted variances
                    MMyy[i] = np.sum((yy00 - MMy)**2 * np.sum(MM, axis=0)) / MMsum[i]
                    MMzz[i] = np.sum((zz00 - MMz)**2 * np.sum(MM, axis=1)) / MMsum[i]
                else:
                    MMyy[i] = np.nan
                    MMzz[i] = np.nan
                    
        else:
            # MATLAB convention: iterate over first dimension (rows)
            logger.debug("Processing with MATLAB convention: iterating over rows")
            
            xx00 = np.arange(1, nrows + 1)  # MATLAB "X" = rows, 1-indexed
            yy00 = np.arange(1, ncols + 1)  # MATLAB "Y" = columns, 1-indexed
            zz00 = np.arange(1, depth + 1)  # Z = slices, 1-indexed
            
            # Preallocate arrays
            MMsum = np.zeros(nrows)
            MMyy = np.zeros(nrows)
            MMzz = np.zeros(nrows)
            
            # Iterate over rows (MATLAB "X" axis)
            for i in range(nrows):
                MM = image_yxz_rotated[i, :, :].T  # Shape: (Z, ncols)
                
                MMsum[i] = np.sum(MM)
                
                if MMsum[i] > 0:
                    # Calculate weighted means
                    MMz = np.sum(zz00 * np.sum(MM, axis=1)) / MMsum[i]
                    MMy = np.sum(yy00 * np.sum(MM, axis=0)) / MMsum[i]
                    
                    # Calculate weighted variances
                    MMyy[i] = np.sum((yy00 - MMy)**2 * np.sum(MM, axis=0)) / MMsum[i]
                    MMzz[i] = np.sum((zz00 - MMz)**2 * np.sum(MM, axis=1)) / MMsum[i]
                else:
                    MMyy[i] = np.nan
                    MMzz[i] = np.nan
        
        logger.debug(f"Local spreads calculated - {len(xx00)} positions")
        logger.debug(f"Total intensity: {np.sum(MMsum):.2f}")
        
        return MMsum, MMyy, MMzz, xx00
    
    
    def calculate_global_spreads(
        self, 
        MMsum: np.ndarray, 
        MMyy: np.ndarray, 
        MMzz: np.ndarray, 
        xx00: np.ndarray
    ) -> Dict[str, float]:
        """
        Calculate global spread metrics from local spreads.
        
        This method aggregates the local spread measurements into global metrics
        that quantify the overall spatial extent of the neuronal projection in
        all three dimensions.
        
        Parameters
        ----------
        MMsum : np.ndarray
            Total intensity at each X position (from calculate_local_spreads)
        MMyy : np.ndarray
            Local Y-variances at each X position
        MMzz : np.ndarray
            Local Z-variances at each X position
        xx00 : np.ndarray
            X coordinate array
        
        Returns
        -------
        Dict[str, float]
            Dictionary containing:
            - spread_x_pixel, spread_y_pixel, spread_z_pixel: Individual axis spreads (pixels)
            - spread_xy_pixel, spread_xyz_pixel: Combined spreads (pixels)
            - spread_x_um, spread_y_um, spread_z_um: Individual axis spreads (µm)
            - spread_xy_um, spread_xyz_um: Combined spreads (µm)
            - axonal_volume: Total integrated intensity
            - MMsum, MMyy, MMzz: Arrays for further analysis/plotting
        
        Notes
        -----
        - All spreads are intensity-weighted averages of local spreads
        - The 3D spread (spread_xyz) is the product of spreads in all three axes
        - Axonal volume represents the total fluorescence signal
        
        Mathematical Details
        --------------------
        - spread_x = sqrt(Σ (x - mean_x)² × intensity / total_intensity)
        - spread_y = sqrt(Σ local_variance_y × intensity / total_intensity)
        - spread_z = sqrt(Σ local_variance_z × intensity / total_intensity)
        - spread_xyz = spread_x × spread_y × spread_z
        
        References
        ----------
        Petsakou et al., Cell 2015 - Supplementary Methods
        """
        total_sum = np.sum(MMsum)
        
        if total_sum == 0:
            logger.error("Total sum is zero - cannot calculate spreads")
            return self._get_empty_results()
        
        # Calculate X-spread (standard deviation of intensity distribution along X)
        MMx = np.sum(xx00 * MMsum) / total_sum
        MMxx = np.sum(MMsum * (xx00 - MMx)**2) / total_sum
        spread_x_px = np.sqrt(MMxx)
        
        # Calculate Y-spread (intensity-weighted average of local Y-variances)
        Total_var_y = np.nansum(MMyy * MMsum) / total_sum
        spread_y_px = np.sqrt(Total_var_y)
        
        # Calculate Z-spread (intensity-weighted average of local Z-variances)
        Total_var_z = np.nansum(MMzz * MMsum) / total_sum
        spread_z_px = np.sqrt(Total_var_z)
        
        logger.debug(f"Global spreads (pixels) - X: {spread_x_px:.2f}, "
                    f"Y: {spread_y_px:.2f}, Z: {spread_z_px:.2f}")
        
        # Convert spreads to physical units (µm)
        spread_x_um = spread_x_px * self.voxel_size_x
        spread_y_um = spread_y_px * self.voxel_size_y
        spread_z_um = spread_z_px * self.voxel_size_z
        
        # Calculate combined spreads
        spread_xy_px = spread_x_px * spread_y_px
        spread_xyz_px = spread_x_px * spread_y_px * spread_z_px
        
        spread_xy_um = spread_x_um * spread_y_um
        spread_xyz_um = spread_x_um * spread_y_um * spread_z_um
        
        # Axonal volume = total integrated intensity
        axonal_volume = total_sum
        
        logger.info(f"Global spreads (µm) - X: {spread_x_um:.2f}, "
                   f"Y: {spread_y_um:.2f}, Z: {spread_z_um:.2f}")
        logger.info(f"3D Spread: {spread_xyz_um:.2f} µm³")
        logger.info(f"Axonal volume: {axonal_volume:.2f}")
        
        # Compile results
        results = {
            'spread_x_pixel': round(spread_x_px, 2),
            'spread_y_pixel': round(spread_y_px, 2),
            'spread_z_pixel': round(spread_z_px, 2),
            'spread_xy_pixel': round(spread_xy_px, 2),
            'spread_xyz_pixel': round(spread_xyz_px, 2),
            'spread_x_um': round(spread_x_um, 2),
            'spread_y_um': round(spread_y_um, 2),
            'spread_z_um': round(spread_z_um, 2),
            'spread_xy_um': round(spread_xy_um, 2),
            'spread_xyz_um': round(spread_xyz_um, 2),
            'axonal_volume': round(axonal_volume, 2),
            # Keep distributions for downstream analysis
            'MMsum': MMsum,
            'MMyy': MMyy,
            'MMzz': MMzz
        }
        
        return results
    
    
    def process_image(
        self, 
        image_3d: np.ndarray,
        mask_area_pixels: float
    ) -> Dict[str, float]:
        """
        Execute the complete image processing pipeline.
        
        This is the main entry point for processing 3D microscopy images.
        It orchestrates all steps: fluorescence calculation, PCA alignment,
        rotation, and spread quantification.
        
        Pipeline Steps:
        1. Calculate fluorescence metrics
        2. Determine optimal rotation angle via PCA
        3. Rotate image to align with principal axis
        4. Apply 90° standardization rotation
        5. Calculate local and global spreads
        
        Parameters
        ----------
        image_3d : np.ndarray
            3D image stack. Should be in (Y, X, Z) format = (rows, cols, slices).
            If image is in (Z, Y, X) format, it will be automatically transposed.
        mask_area_pixels : float
            Region of interest area in pixels (2D projection area)
        
        Returns
        -------
        Dict[str, float]
            Complete set of metrics including:
            - All spread measurements (pixels and µm)
            - Fluorescence metrics
            - Rotation angles applied
            - Axonal volume
            - Distribution arrays (MMsum, MMyy, MMzz)
        
        Raises
        ------
        ValueError
            If image contains negative values (indicates preprocessing error)
        
        Notes
        -----
        - Image format is automatically detected and corrected if needed
        - The pipeline preserves total intensity (small losses due to interpolation)
        - All intermediate steps are logged for debugging
        
        Examples
        --------
        >>> processor = ImageProcessor(0.124, 0.124, 1.0)
        >>> results = processor.process_image(image_3d, mask_area=5000)
        >>> print(f"3D spread: {results['spread_xyz_um']:.2f} µm³")
        """
        logger.info("="*60)
        logger.info("Starting structural plasticity analysis pipeline")
        logger.info(f"Input shape: {image_3d.shape}")
        
        # Ensure correct image format: (Y, X, Z) = (rows, cols, slices)
        if image_3d.shape[0] < image_3d.shape[2]:
            # Likely (Z, Y, X) format - transpose to (Y, X, Z)
            image_yxz = image_3d.transpose(1, 2, 0)
            logger.info(f"Transposed from (Z,Y,X) to (Y,X,Z): {image_yxz.shape}")
        else:
            image_yxz = image_3d
            logger.info(f"Using image as-is (assumed Y,X,Z format): {image_yxz.shape}")
        
        # Validate image (no negative values should exist)
        if np.any(image_yxz < 0):
            raise ValueError(f"Image contains {np.sum(image_yxz < 0)} negative values!")
        
        # Step 1: Calculate fluorescence metrics
        logger.info("Step 1/4: Calculating fluorescence...")
        fluor_px, fluor_um = self.calculate_fluorescence(image_yxz, mask_area_pixels)
        
        # Step 2: Calculate optimal rotation angle using PCA
        logger.info("Step 2/4: Calculating PCA rotation angle...")
        angle = self.calculate_pca_rotation_angle(image_yxz)
        
        # Step 3: Rotate image to align with principal axis
        logger.info("Step 3/4: Rotating image...")
        image_rotated = self.rotate_image(image_yxz, angle)
        
        # Validate rotated image
        if np.any(image_rotated < 0):
            raise ValueError(f"Rotated image contains {np.sum(image_rotated < 0)} negative values!")
        
        # Log intensity conservation
        intensity_loss = (1 - np.sum(image_rotated) / np.sum(image_yxz)) * 100
        logger.debug(f"Intensity loss after rotation: {intensity_loss:.2f}%")
        
        # Step 3.5: Apply 90° standardization rotation
        # This ensures X is horizontal (elongated) and Y is vertical (narrow)
        logger.info("Applying 90° standardization rotation (X=horizontal, Y=vertical)...")
        image_rotated_90 = np.rot90(image_rotated, k=1, axes=(0, 1))
        
        logger.debug(f"Shape after 90° rotation: {image_rotated_90.shape}")
        logger.debug(f"Final coordinate convention: X=horizontal (elongated), Y=vertical (narrow)")
        
        # Step 4: Calculate spread metrics
        logger.info("Step 4/4: Calculating spread metrics...")
        MMsum, MMyy, MMzz, xx00 = self.calculate_local_spreads(
            image_rotated_90, 
            after_90deg_rotation=True
        )
        spread_results = self.calculate_global_spreads(MMsum, MMyy, MMzz, xx00)
        
        logger.info(f"Final spreads - X: {spread_results['spread_x_um']:.2f}µm (horizontal), "
                   f"Y: {spread_results['spread_y_um']:.2f}µm (vertical), "
                   f"Z: {spread_results['spread_z_um']:.2f}µm")
        
        # Add fluorescence and metadata to results
        spread_results['fluorescence_px'] = fluor_px
        spread_results['fluorescence_um'] = fluor_um
        spread_results['rotation_angle'] = round(angle, 2)
        spread_results['additional_rotation'] = 90  # Standardization rotation
        
        logger.info("Processing completed successfully")
        logger.info("="*60)
        
        return spread_results
    
    
    def _get_empty_results(self) -> Dict[str, float]:
        """
        Return empty results dictionary for error cases.
        
        Returns
        -------
        Dict[str, float]
            Dictionary with all metrics set to 0.0 or empty arrays
        """
        return {
            'spread_x_pixel': 0.0,
            'spread_y_pixel': 0.0,
            'spread_z_pixel': 0.0,
            'spread_xy_pixel': 0.0,
            'spread_xyz_pixel': 0.0,
            'spread_x_um': 0.0,
            'spread_y_um': 0.0,
            'spread_z_um': 0.0,
            'spread_xy_um': 0.0,
            'spread_xyz_um': 0.0,
            'axonal_volume': 0.0,
            'fluorescence_px': 0.0,
            'fluorescence_um': 0.0,
            'rotation_angle': 0.0,
            'MMsum': np.array([]),
            'MMyy': np.array([]),
            'MMzz': np.array([])
        }


def validate_parameters(
    voxel_size_x: float, 
    voxel_size_y: float, 
    voxel_size_z: float
) -> Tuple[bool, str]:
    """
    Validate voxel size parameters for physical plausibility.
    
    Parameters
    ----------
    voxel_size_x, voxel_size_y, voxel_size_z : float
        Voxel dimensions in micrometers (µm)
    
    Returns
    -------
    Tuple[bool, str]
        - is_valid: True if parameters pass all checks
        - error_message: Description of validation failure (empty if valid)
    
    Notes
    -----
    Validation checks:
    - All values must be positive
    - XY dimensions should be < 10µm (typical for high-resolution confocal)
    - Z dimension should be < 50µm (typical z-step size)
    
    Examples
    --------
    >>> is_valid, msg = validate_parameters(0.124, 0.124, 1.0)
    >>> print(is_valid)  # True
    >>> is_valid, msg = validate_parameters(-1, 0.124, 1.0)
    >>> print(msg)  # "All voxel sizes must be positive"
    """
    if voxel_size_x <= 0 or voxel_size_y <= 0 or voxel_size_z <= 0:
        return False, "All voxel sizes must be positive"
    
    if voxel_size_x > 10 or voxel_size_y > 10:
        return False, "XY voxel sizes > 10µm seem unusual. Please verify your calibration."
    
    if voxel_size_z > 50:
        return False, "Z voxel size > 50µm seems unusual. Please verify your z-step."
    
    return True, ""