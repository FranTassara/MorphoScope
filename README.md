# MorphoScope

A Python-based GUI application for quantifying structural plasticity in neuronal projections, implementing the methodology from Petsakou et al., 2015 (*Cell*).

## Overview

This tool provides automated 3D quantification of neuronal structural plasticity from confocal microscopy images. It calculates:
- **3D spread** measurements (X, Y, Z axes)
- **Axonal volume**
- **Fluorescence intensity**
- **PCA-based rotation** for objective axis alignment

The algorithm is based on the research published in:
> **Petsakou A, Sapsis T & Blau J (2015).** *Circadian rhythms in Rho1 activity regulate neuronal plasticity and network hierarchy.* Cell, 162(4):823-835.

## Features

- **Multi-format support**: TIFF, CZI, and JPEG image stacks
- **Interactive ROI selection**: Polygon-based region of interest drawing
- **Real-time preview**: Visualize raw and filtered images
- **Advanced filtering**: Gaussian blur, median filter, and thresholding
- **Batch processing**: Process multiple images sequentially
- **CSV export**: Export all measurements to spreadsheet-compatible format

## Screenshots

### Main Interface
*[Add screenshot of your main interface]*

### ROI Selection
*[Add screenshot showing ROI polygon drawing]*

### Results Export
*[Add screenshot of results table]*

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
   - Calculate centroid of intensity-weighted projection
   - Compute covariance matrix
   - Rotate image to align maximum spread with X-axis (objective coordinate system)

3. **Mean 3D Curve Calculation**
   - For each X position, calculate intensity-weighted mean Y and Z positions
   - This defines the "spine" of the projection

4. **Local Spread Calculation**
   - For each X position, compute typical deviation (weighted variance) in Y and Z
   - Weighting by fluorescence intensity (Cf) ensures bright regions contribute more

5. **Global Spread Metrics**
   - X-spread: Intensity-weighted standard deviation along X-axis
   - Y-spread: Weighted average of local Y typical deviations
   - Z-spread: Weighted average of local Z typical deviations
   - **3D spread**: Product of X, Y, and Z spreads (σ̄ₓ × σ̄ᵧ × σ̄ᵧ)

6. **Additional Metrics**
   - **Axonal volume (M)**: Sum of all intensity values (total fluorescent material)
   - **Fluorescence**: [COMPLETAR]

### Mathematical Formulation

Where Cf(x, y, z) represents fluorescence intensity at each voxel.

#### 3D Mean Curve

For each x position xⱼ, the intensity-weighted mean Y and Z coordinates define the central path:

$$P_y(x_j) = \frac{\sum_i y_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

$$P_z(x_j) = \frac{\sum_i z_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

This 3D curve (x, Py(x), Pz(x)) represents the "backbone" of the axonal projection.

#### Local Spreads (Typical Deviations)

The local typical deviations measure spread around the mean curve at each x position:

$$\sigma_y(x_j) = \sqrt{\frac{\sum_i (y_i - P_y(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

$$\sigma_z(x_j) = \sqrt{\frac{\sum_i (z_i - P_z(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

These are intensity-weighted standard deviations at each X slice.

#### Global Spreads

**Y-axis and Z-axis global spreads** - weighted averages of local spreads:

$$\bar{\sigma}_y = \frac{\sum_i \sigma_y(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

$$\bar{\sigma}_z = \frac{\sum_i \sigma_z(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

**X-axis global spread** - typical deviation of the entire distribution along X:

$$\bar{\sigma}_x = \sqrt{\frac{\sum_i (x_i - m_x)^2 \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}}$$

where the intensity-weighted center of mass along X is:

$$m_x = \frac{\sum_i x_i \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

**Axonal Volume (M)** - total fluorescent material:

$$M = \sum_i C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z$$

**3D Spread** - final metric combining all three axes:

$$\text{3D Spread} = \bar{\sigma}_x \times \bar{\sigma}_y \times \bar{\sigma}_z$$

Where:
- δx, δy, δz are voxel dimensions (spacing between points)
- Cf(x, y, z) is fluorescence intensity at each voxel
- M is the total axonal volume (sum of all intensities × voxel volume)

#### Fluorescence Density [COMPLETAR]

$$\text{Fluorescence} = \frac{M}{V_{cuboid}}$$

where Vcuboid is the volume of the bounding box containing the ROI.

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
| Image Name | Source filename | - |
| Spread X | X-axis spread | pixels / μm |
| Spread Y | Y-axis spread | pixels / μm |
| Spread Z | Z-axis spread | pixels / μm |
| Spread XY | Product of X and Y spreads | pixels² / μm² |
| Spread XYZ | 3D spread (product of all three) | pixels³ / μm³ |
| Axonal Volume | Sum of all intensity values | a.u. |
| Fluorescence | Normalized fluorescence | pixels / μm² |
| Observation | User notes | - |


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
- **Python Implementation**: [Francisco Joaquín Tassara]

## Acknowledgments

- Original MATLAB implementation by Petsakou et al. (2015)
- Open-source Python scientific computing community

---
