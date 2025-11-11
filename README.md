# MorphoScope

A Python-based GUI application for quantifying structural plasticity in neuronal projections, implementing the methodology from Petsakou et al., 2015 (*Cell*).

## 📋 Overview

This tool provides automated 3D quantification of neuronal structural plasticity from confocal microscopy images. It calculates:
- **3D spread** measurements (X, Y, Z axes)
- **Axonal volume**
- **Fluorescence intensity**
- **PCA-based rotation** for objective axis alignment

The algorithm is based on the research published in:
> **Petsakou A, Sapsis T & Blau J (2015).** *Circadian rhythms in Rho1 activity regulate neuronal plasticity and network hierarchy.* Cell, 162(4):823-835.

## ✨ Features

- **Multi-format support**: TIFF, CZI, and JPEG image stacks
- **Interactive ROI selection**: Polygon-based region of interest drawing
- **Real-time preview**: Visualize raw and filtered images
- **Advanced filtering**: Gaussian blur, median filter, and thresholding
- **Batch processing**: Process multiple images sequentially
- **CSV export**: Export all measurements to spreadsheet-compatible format

## 🖼️ Screenshots

### Main Interface
*[Add screenshot of your main interface]*

### ROI Selection
*[Add screenshot showing ROI polygon drawing]*

### Results Export
*[Add screenshot of results table]*

## 📦 Installation

### Prerequisites

- Python 3.8 or higher
- PySide6 (Qt for Python)
- Scientific computing libraries

### Setup

1. **Clone the repository**
```bash
git clone https://github.com/YOUR_USERNAME/structural-plasticity-analyzer.git
cd structural-plasticity-analyzer
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
python main.py
```

## 📚 Dependencies

```
PySide6>=6.5.0
numpy>=1.24.0
scipy>=1.10.0
matplotlib>=3.7.0
tifffile>=2023.0.0
pylibCZIrw>=3.0.0
scikit-image>=0.20.0
pyqtgraph>=0.13.0
shapely>=2.0.0
```

## 🚀 Quick Start

### Basic Workflow

1. **Load Images**: Click "Load Image" and select your confocal stack (TIFF/CZI)
2. **Set Voxel Dimensions**: Enter X, Y, Z voxel sizes in micrometers
3. **Draw ROI**: Click "Select ROI" and draw polygon around region of interest
4. **Apply Filters** (optional): 
   - Gaussian blur for noise reduction
   - Median filter for salt-and-pepper noise
   - Threshold to remove background
5. **Select Z-Range**: Choose which slices to analyze (default: all)
6. **Process**: Click "Process" to run the analysis
7. **Export Results**: Save measurements to CSV file

### Image Requirements

- **Format**: Confocal Z-stacks (TIFF, CZI, or series of JPEGs)
- **Recommended resolution**: 1024×1024 pixels
- **Bit depth**: 8-bit or 16-bit
- **Z-step size**: Consistent spacing (typically 1 μm)
- **Imaging**: Single channel or multi-channel (specify channel for analysis)

## 📖 Methodology

### Algorithm Overview

The analysis pipeline follows these steps:

1. **Image Loading & Preprocessing**
   - Load 3D image stack
   - Apply optional filters (Gaussian, median, threshold)
   - Extract user-defined ROI

2. **PCA-Based Rotation**
   - Calculate centroid of Z-projection
   - Compute covariance matrix
   - Rotate image to align maximum spread with X-axis

3. **Local Spread Calculation**
   - For each X position, calculate intensity-weighted mean Y and Z positions
   - Compute variance around these means

4. **Global Spread Metrics**
   - X-spread: Standard deviation of X distribution
   - Y-spread: Weighted average of local Y variances
   - Z-spread: Weighted average of local Z variances
   - 3D spread: Product of X, Y, and Z spreads

5. **Additional Metrics**
   - Axonal volume: Sum of all intensity values
   - Fluorescence: Normalized intensity per unit area

### Mathematical Formulation

#### 3D Mean Curve
For each x position, the mean Y and Z coordinates are:

$$P_y(x_j) = \frac{\sum_i y_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

$$P_z(x_j) = \frac{\sum_i z_i \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}$$

#### Local Spreads
The local standard deviations at each x position:

$$\sigma_y(x_j) = \sqrt{\frac{\sum_i (y_i - P_y(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

$$\sigma_z(x_j) = \sqrt{\frac{\sum_i (z_i - P_z(x_j))^2 \cdot C_f(x_j, y_i, z_i)}{\sum_i C_f(x_j, y_i, z_i)}}$$

#### Global Spreads
Weighted averages across all x positions:

$$\bar{\sigma}_y = \frac{\sum_i \sigma_y(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

$$\bar{\sigma}_z = \frac{\sum_i \sigma_z(x_i) \cdot C_f(x_i, y_i, z_i) \cdot \delta x \cdot \delta y \cdot \delta z}{M}$$

Where M is the total axonal volume (sum of all intensity values).

## 📊 Output Format

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

## 🔬 Original MATLAB Implementation

This Python implementation modernizes the original MATLAB script from Petsakou et al. (2015) with:
- ✅ User-friendly GUI
- ✅ No hardcoded filters (user-controlled)
- ✅ Native support for modern file formats (CZI, multi-page TIFF)
- ✅ Batch processing capabilities
- ✅ Improved error handling and logging

The original `spread.m` MATLAB script is included in the `/matlab` directory for reference.

## 📝 Citation

If you use this tool in your research, please cite the original methodology:

```bibtex
@article{petsakou2015circadian,
  title={Circadian rhythms in Rho1 activity regulate neuronal plasticity and network hierarchy},
  author={Petsakou, Afroditi and Sapsis, Themistoklis P and Blau, Justin},
  journal={Cell},
  volume={162},
  number={4},
  pages={823--835},
  year={2015},
  publisher={Elsevier}
}
```

And optionally cite this software:

```bibtex
@software{structural_plasticity_analyzer,
  title={Structural Plasticity Analyzer},
  author={YOUR_NAME},
  year={2025},
  url={https://github.com/YOUR_USERNAME/structural-plasticity-analyzer}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Original Algorithm**: Afroditi Petsakou, Themistoklis P. Sapsis, Justin Blau
- **Python Implementation**: [Your Name]

## 🙏 Acknowledgments

- Original MATLAB implementation by Petsakou et al. (2015)
- NYU Biology Department for the foundational research
- Open-source Python scientific computing community

## 📧 Contact

For questions or support, please:
- Open an issue on GitHub
- Email: [your.email@example.com]

## 🔗 Related Resources

- [Original Cell paper (2015)](https://doi.org/10.1016/j.cell.2015.07.010)
- [Circadian clock research at NYU](https://as.nyu.edu/biology.html)
- [Python image processing with scikit-image](https://scikit-image.org/)

---

**Note**: This is a research tool. Results should be validated against the original methodology and appropriate controls.
