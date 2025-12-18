# MorphoScope

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17901990.svg)](https://doi.org/10.5281/zenodo.17901990)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python-based GUI application for quantifying structural plasticity in neuronal projections, implementing the methodology from Petsakou et al., 2015 (*Cell*).

## Overview

This tool provides automated 3D quantification of neuronal structural plasticity from confocal microscopy images. It calculates:
- **3D spread** measurements (X, Y, Z axes)
- **Axonal volume** (integrated intensity)
- **Geometric volume** (physical volume of occupied voxels)
- **Fluorescence density** (normalized intensity)
- **PCA-based rotation** for objective axis alignment

The algorithm is based on the research published in:
> **Petsakou A, Sapsis T & Blau J (2015).** *Circadian rhythms in Rho1 activity regulate neuronal plasticity and network hierarchy.* Cell, 162(4):823-835. [DOI: 10.1016/j.cell.2015.07.010](https://doi.org/10.1016/j.cell.2015.07.010)

## Features

- **Multi-format support**: TIFF, CZI, and JPEG image stacks
- **Interactive ROI selection**: Polygon-based region of interest drawing
- **Real-time preview**: Visualize raw and filtered images
- **Advanced filtering**: Gaussian blur, median filter, and thresholding
- **Dual-channel analysis**: Separate channels for morphology and protein quantification
- **Batch processing**: Process multiple images sequentially
- **CSV export**: Export all measurements to spreadsheet-compatible format

## Installation

### Option 1: Standalone Executable (Windows)

Download the latest release from the [Releases](https://github.com/FranTassara/MorphoScope/releases/tag/v1.0.0) page:
```bash
# Download Morphoscope.exe
# Double-click to run - no installation needed!
```

### Option 2: Python Installation

#### Prerequisites

- Python 3.8 or higher
- PySide6 (Qt for Python)
- Scientific computing libraries

#### Setup

1. **Clone the repository**
```bash
git clone https://github.com/FranTassara/MorphoScope.git
cd MorphoScope
```

2. **Create virtual environment** (recommended)
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

4. **Run the application**
```bash
python MorphoScope.py
```

## Quick Start

### Basic Workflow

1. **Load Images**: Click "Load Image" and select your confocal stack (TIFF/CZI)
2. **Set Voxel Dimensions**: Enter X, Y, Z voxel sizes in micrometers
3. **Apply Filters** (optional): 
   - Gaussian blur for noise reduction
   - Median filter for salt-and-pepper noise
   - Threshold to remove background
4. **Draw ROI**: Click "Select ROI" and draw polygon around region of interest
5. **Select Z-Range**: Choose which slices to analyze (default: all)
6. **Process**: Click "Process" to run the analysis
7. **Export Results**: Save measurements to CSV file

### Image Requirements

- **Format**: Confocal Z-stacks (TIFF, CZI, or LSM)
- **Recommended resolution**: 1024×1024 pixels
- **Bit depth**: 8-bit or 16-bit
- **Z-step size**: Consistent spacing (typically 1 μm)
- **Imaging**: Single channel or multi-channel (specify channel for analysis)

## Methodology

### Algorithm Overview

The analysis pipeline follows these steps:

1. **Image Loading & Preprocessing**
   - Load 3D image stack
   - Apply optional filters (Gaussian, median, threshold)
   - Extract user-defined ROI

2. **PCA-Based Rotation**
   - Calculate centroid of intensity-weighted Z-projection
   - Compute covariance matrix
   - Rotate image to align maximum spread with X-axis (objective coordinate system)
   - Apply additional 90° standardization rotation

3. **Mean 3D Curve Calculation**
   - For each X position, calculate intensity-weighted mean Y and Z positions
   - This defines the "spine" of the projection

4. **Local Spread Calculation**
   - For each X position, compute typical deviation (weighted standard deviation) in Y and Z
   - Weighting by fluorescence intensity (C_f) ensures bright regions contribute more

5. **Global Spread Metrics**
   - X-spread: Intensity-weighted standard deviation along X-axis
   - Y-spread: Weighted average of local Y typical deviations
   - Z-spread: Weighted average of local Z typical deviations
   - **3D spread**: Product of X, Y, and Z spreads (σ̄ₓ × σ̄ᵧ × σ̄_z)

6. **Volume and Fluorescence Metrics**
   - **Axonal volume (M)**: Sum of all intensity values × voxel volume (total fluorescent material)
   - **Geometric volume**: Physical volume of occupied voxels in µm³
   - **Fluorescence density**: Intensity normalized by occupied area

### Mathematical Formulation

Where C_f(x, y, z) represents fluorescence intensity at each voxel.

#### 3D Mean Curve

For each x position x_j, the intensity-weighted mean Y and Z coordinates define the central path:

$$P_y(x_j) = \frac{\sum_i y_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

$$P_z(x_j) = \frac{\sum_i z_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

This 3D curve (x, P_y(x), P_z(x)) represents the "backbone" of the axonal projection.

#### Local Spreads (Typical Deviations)

The local typical deviations measure spread around the mean curve at each x position:

$$\sigma_y(x_j) = \sqrt{\frac{\sum_i (y_i - P_y(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

$$\sigma_z(x_j) = \sqrt{\frac{\sum_i (z_i - P_z(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

These are intensity-weighted standard deviations at each X slice.

#### Global Spreads

**Y-axis and Z-axis global spreads** — weighted averages of local spreads:

$$\bar{\sigma}_y = \frac{\sum_i \sigma_y(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

$$\bar{\sigma}_z = \frac{\sum_i \sigma_z(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

**X-axis global spread** — typical deviation of the entire distribution along X:

$$\bar{\sigma}_x = \sqrt{\frac{\sum_i (x_i - m_x)^2 \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}}$$

where the intensity-weighted center of mass along X is:

$$m_x = \frac{\sum_i x_i \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

**3D Spread** — final metric combining all three axes:

$$\text{3D Spread} = \bar{\sigma}_x \times \bar{\sigma}_y \times \bar{\sigma}_z$$

#### Axonal Volume (Integrated Intensity)

Total fluorescent material, equivalent to "M" in Petsakou et al.:

$$M = \sum_i C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z$$

This represents the total amount of fluorescent signal, weighted by voxel volume.

#### Geometric Volume

Physical volume occupied by the neuronal projection:

$$V_{geometric} = N_{occupied} \cdot \delta x \cdot \delta y \cdot \delta z$$

Where N_occupied is the number of non-zero voxels. This metric quantifies the actual 3D space occupied by the projection, independent of fluorescence intensity.

#### Fluorescence Density

The fluorescence metric quantifies the concentration of fluorescent material within the projection. This implementation uses a maximum intensity Z-projection approach combined with automatic background thresholding:

$$\text{Fluorescence}_{px} = \frac{\sum_{(x,y) \in \Omega} \max_z(C_f(x,y,z))}{N_{occupied}}$$

$$\text{Fluorescence}_{\mu m^2} = \frac{\sum_{(x,y) \in \Omega} \max_z(C_f(x,y,z))}{N_{occupied} \cdot \delta x \cdot \delta y}$$

Where:
- max_z(C_f(x,y,z)) is the maximum intensity Z-projection
- $\Omega$ is the set of pixels where intensity > $T_{triangle}$ (Triangle threshold).
- N_occupied is the number of pixels in $\Omega$ (pixels above threshold).
- δx, δy are voxel dimensions in µm

**Note**: This differs from the original MATLAB implementation which normalized by the full ROI cuboid volume. The current approach provides a more robust measure that:
1. Uses maximum projection to reduce background contributions
2. Counts only occupied pixels, making it independent of ROI shape
3. Reports both per-pixel (AU/pixel) and per-area (AU/µm²) values

#### Symbol Reference

| Symbol | Description |
|--------|-------------|
| C_f(x,y,z) | Fluorescence intensity at voxel (x,y,z) |
| δx, δy, δz | Voxel dimensions (spacing between points) in µm |
| M | Axonal volume (total integrated intensity × voxel volume) |
| P_y, P_z | Intensity-weighted mean positions (3D curve) |
| σ_y, σ_z | Local typical deviations at each x position |
| σ̄_x, σ̄_y, σ̄_z | Global spreads along each axis |
| N_occupied | Number of non-zero voxels |

### Key Differences from Standard Metrics

⚠️ **Important**: These are not simple standard deviations:

1. **Intensity-weighted**: Bright regions contribute more to the measurements
2. **Asymmetric treatment**: X-axis uses global variance, while Y and Z use averaged local variances
3. **Curved structures**: The 3D mean curve accounts for non-linear axonal paths
4. **Normalized by mass**: Global spreads are normalized by total fluorescence (M)

This approach is specifically designed for **curved, non-uniform neuronal projections** where simple bounding box or centroid-based measurements would be inaccurate.

## Output Format

Results are exported to CSV with the following columns:

| Column | Description | Units |
|--------|-------------|-------|
| Image filename | Source filename | — |
| Spread x [pixel] | X-axis spread (horizontal, maximum length direction) | pixels |
| Spread y [pixel] | Y-axis spread (vertical, perpendicular to X) | pixels |
| Spread z [pixel] | Z-axis spread (depth, dorsal-ventral) | pixels |
| Spread x*y [pixel²] | Product of X and Y spreads | pixels² |
| Spread x*y*z [pixel³] | 3D spread (product of all three) | pixels³ |
| Spread x [µm] | X-axis spread in physical units | µm |
| Spread y [µm] | Y-axis spread in physical units | µm |
| Spread z [µm] | Z-axis spread in physical units | µm |
| Spread x*y [µm²] | Product of X and Y spreads | µm² |
| Spread x*y*z [µm³] | 3D spread in physical units | µm³ |
| Axonal Volume (integrated intensity) | Sum of all intensity values × voxel volume | AU × µm³ |
| Geometric volume [µm³] | Physical volume of occupied voxels | µm³ |
| Fluorescence_px [AU/pixel] | Mean intensity per occupied pixel | AU/pixel |
| Fluorescence_um [AU/µm²] | Mean intensity per physical area | AU/µm² |
| Observation | User notes | — |

### Metric Interpretation

- **3D Spread (σ̄_x × σ̄_y × σ̄_z)**: Measures spatial dispersion of the projection. Higher values indicate more spread/defasciculated projections. *Note: This is NOT a physical volume measurement.*

- **Axonal Volume (M)**: Total fluorescent content. Useful for comparing total protein/marker levels between conditions.

- **Geometric Volume**: Actual physical space occupied. Independent of fluorescence intensity.

- **Fluorescence Density**: Concentration of fluorescent signal. Useful for ensuring comparable staining (PDF) between experimental and control samples.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

- **Original Algorithm**: Afroditi Petsakou, Themistoklis P. Sapsis, Justin Blau
- **Python Implementation**: Francisco Joaquín Tassara

## Acknowledgments

- Original MATLAB implementation by Petsakou et al. (2015)
- Open-source Python scientific computing community (NumPy, SciPy, PySide6)

---

*For questions or issues, please open an issue on GitHub.*
