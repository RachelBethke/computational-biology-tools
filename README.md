# Computational Biology Tools for Analysis

A Python toolkit for computational biology analysis, providing interactive interfaces for various biological data analysis tools.

## Key Features
- Interactive GUIs for complex biological analyses
- Real-time visualization of results
- Modular design for multiple analysis types
- Cross-platform compatibility

## Current Implementation: Conservation Genetics Module
The first module focuses on conservation genetics analysis, featuring:
- Loading and analysis of haplotype data
- Calculation of key genetic diversity metrics:
  - Nucleotide diversity (π)
  - Watterson's theta
  - Tajima's D
- Interactive visualization of allele frequencies
- GUI interface for easy data exploration

## Getting Started

### Prerequisites
Ensure you have the following installed:
- Python 3.7 or higher
- `tkinter` (may need to be installed separately for Linux: `sudo apt-get install python3-tk`)

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/RachelBethke/computational-biology-tools.git
   cd computational-biology-tools
   ```
2. Install dependencies:
   ```bash
   # Install the package
   pip install -e .

   # Run the GUI
   python -m genetics
   ```

### Setting Up a Virtual Environment (Recommended)
I used a virtual environment with this program so I can keep the dependencies isolated, however this isn't strictly necessary to use it.
1. Create a virtual environment:
   ```bash
   python -m venv venv
   ```
2. Activate the virtual environment:
   - **Windows**: `venv\Scripts\activate`
   - **macOS/Linux**: `source venv/bin/activate`

Following these steps, you can launch the GUI from your terminal:
```bash
python -m genetics
```

## Current Development Plans

### Phase 1: Conservation Genetics (Current)
- [x] Core diversity metrics
- [x] Core calculation testing
- [x] Basic GUI implementation
- [x] Single plot visualization
- [ ] Block selection interface
- [ ] Additional conservation metrics (?)
- [ ] Visualization options and selection
- [ ] GUI testing
- [ ] Data export capabilities

### Interlude: Toolkit Structure
- [ ] Restructure project for multiple tools
- [ ] Tool selection GUI
- [ ] Universal data export capabilities
- [ ] Shared visualization tools

### Phase 2: Population Dynamics
- [ ] Model options:
  - [ ] Exponential Growth Model
  - [ ] Logistic Growth Model
  - [ ] Predator-Prey Dynamics (Lotka-Volterra Model)
- [ ] Interactive parameter adjustment through GUI
- [ ] Real-time population visualization
- [ ] Support for multiple population models
- [ ] Customizable simulation length and time steps

### Phase 3: Genetic Drift Simulator
- [ ] Core simulation engine:
  - [ ] Wright-Fisher implementation
  - [ ] Parameter controls
  - [ ] Basic testing framework
- [ ] Basic GUI:
  - [ ] Population size input
  - [ ] Initial frequency selection
  - [ ] Run controls
- [ ] Results visualization:
  - [ ] Trajectory plotting
  - [ ] Multi-run comparison
  - [ ] Bottleneck analysis

### Beyond: Possible Future Modules
- Molecular Evolution
- Reference Genome Mapping and Sequence Analysis
- Ancestry and Admixture
- Comparative Genomics
- Brain Network Analysis
- Gene Expression Analysis

## Requirements
- Python 3.7+
- `numpy`
- `pandas`
- `matplotlib`
- `tkinter`

## Acknowledgements
Many of the calculations in the `core.py` files are based on programs I developed for Cornell's BIOCB 2010: Intro to Computational Biology. The course provided the foundations and practical exercises that inspired this project. While the code for these calculations is original, I was given starter code to help me on my way. These set ups that were given to me shaped the design of each practical and pushed me towards the correct implementation. I have made every effort to ensure that all equations are either my own or properly cited when necessary.

## Project Status
This project is currently in active development, with the conservation genetics module as the primary focus. Features and modules will be added incrementally, with regular updates and improvements to existing functionality.