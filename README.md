# Standard Curve Calculator

A Python tool that calculates concentrations from your measurement data using standard curves. Perfect for analytical chemistry, biochemistry, and molecular biology applications where you need to determine unknown concentrations using standard curves.

## Core Features

- **Standard Curve Analysis**:
  - Fits your standard data points to create a calibration curve
  - Multiple curve fitting options:
    - Linear regression (y = mx + b)
    - 4PL (Four Parameter Logistic: y = d + (a-d)/(1 + (x/c)^b))
    - 5PL (Five Parameter Logistic: adds asymmetry parameter)
    - Auto-selection of best fit
  - Statistical validation of curve fit (R², RMSE)

- **Concentration Calculation**:
  - Interpolates unknown concentrations from your measurements
  - Handles multiple replicates
  - Provides statistical confidence intervals
  - Flags out-of-range values
  - Includes measurement uncertainty

- **Quality Control**:
  - Validates standard curve linearity/fit
  - Checks replicate consistency
  - Monitors curve regression quality
  - Alerts for extrapolation

## Requirements

Core dependencies:
```
pandas
numpy
scipy
matplotlib
openpyxl
```

## Installation

1. Clone repository:
```bash
git clone https://github.com/yourusername/standard-curve-calculator.git
cd standard-curve-calculator
```

2. Install dependencies:
```bash
pip install -r requirements.txt

For conda users:
conda env create -f standard-curve-calculator.yml
conda activate standard-curve-calculator
```

## Basic Workflow

1. **Input**: Provide your data with:
   - Standard concentrations and their measurements
   - Unknown sample measurements

2. **Processing**:
   - Fits standard curve
   - Calculates concentrations
   - Performs statistical analysis

3. **Output**:
   - Calculated concentrations
   - Standard curve plot
   - Statistical report
   - Quality checks

## Usage

### Quick Start

```bash
python std_curve_calc.py -i your_data.xlsx -m auto
```

### Detailed Options

```bash
python std_curve_calc.py --input data.xlsx \
                        --sheet "Data" \
                        --method 4 \
                        --output results.xlsx \
                        --verbose
```

### Parameters

- `--input`, `-i`: Your data file (Excel/CSV)
- `--method`, `-m`: Curve type:
  - `linear`: Linear regression
  - `4`: 4PL curve
  - `5`: 5PL curve
  - `auto`: Automatic selection
- `--output`, `-o`: Results file path
- `--verbose`, `-v`: Show calculation details

## Input Format

The tool accepts two formats for your data:

### Basic Format (Excel/CSV)

Your input file should have two sections separated by an empty row:
1. Standards section (top)
2. Samples section (bottom)

```
Standards    0       0.5     1       2       4         # Concentrations
            0.001   0.125   0.250   0.500   1.000     # Replicate 1
            0.002   0.128   0.255   0.505   0.995     # Replicate 2
            0.001   0.122   0.248   0.495   0.998     # Replicate 3

Samples        Measurement
Sample_A1      0.350           # Will be grouped as replicates
Sample_A2      0.355           # of "Sample"
Sample-Rep1    0.348           # Same here
Unknown1       0.750           # Single measurement
TestB_1        0.480           # Will be grouped with TestB_2
TestB_2        0.485           # as replicates of "TestB"
```

Key points:
- First row: Standard concentrations
- Following rows before empty line: Standard measurements (replicates)
- After empty line: Sample measurements
- Headers can be customized but structure must be maintained

### Technical Replicate Detection

#### For Standards
- All rows before the first empty row are considered technical replicates
- Number of replicate rows determines `n` for statistical calculations
- Each standard concentration must have same number of replicates
- Tool automatically calculates mean, SD, and CV% for each concentration

#### For Samples
Tool automatically detects replicates in two ways:

1. **Identical Names**: Samples with exactly the same name are treated as replicates
2. **Pattern Matching**: Samples following specific naming patterns are grouped as replicates

Example of identical names:
```
Sample Name    Measurement    → Result
Sample         0.350         → Grouped together
Sample         0.355         → as replicates
Sample         0.348         → of "Sample"
```

```
Original Names    →   Grouped As
Sample_1              Sample
Sample_2              Sample
Sample-A              Sample
Sample-B              Sample
Sample(1)             Sample
Sample(2)             Sample
Sample_Rep1           Sample
Sample_Rep2           Sample
Test 1                Test
Test 2                Test
```

Supported patterns:
- `name_number`: Sample_1, Sample_2
- `name-number`: Sample-1, Sample-2
- `name(number)`: Sample(1), Sample(2)
- `name number`: Sample 1, Sample 2
- `name_Repnumber`: Sample_Rep1, Sample_Rep2
- `name-Repnumber`: Sample-Rep1, Sample-Rep2

Important Notes:
- **Replicate Detection**: 
  - Identical names: "Sample" and "Sample" are always treated as replicates
  - Pattern matching: "Sample_1" and "Sample_2" are grouped as "Sample"
- **Case Sensitivity**: "Sample_1" and "sample_1" are treated as different groups
- **Special Characters**: "IL-6_1" and "IL-6_2" group as "IL-6"
- **Numbers Only**: "1", "2", "3" do not group together
- **Unintended Grouping Warning**: Be careful with sample naming as identical names will be treated as replicates even if not intended

Example replicate statistics output:
```
Original Names        Grouped As    N    Mean    SD     CV%
Sample_1,Sample_2     Sample       2    0.352   0.004   1.1
Test1,Test2          Test         2    0.482   0.003   0.6
Unknown1             Unknown1     1    0.750   N/A     N/A
```

### Quality Control
The tool performs automated QC checks on replicates:
- Warns if standards have missing replicates
- Flags samples with high CV% (default threshold: 15%)
- Notes samples with single measurements
- Alerts if replicate groups are uneven

## Output Examples

1. **Concentration Results**:
```
Sample     Measurement    Calc.Conc    95%CI     Status
Unknown1   1.250         25.3 ng/mL   ±2.1      PASS
Unknown2   0.750         12.1 ng/mL   ±1.5      PASS
Unknown3   2.100         250.5 ng/mL  ±15.2     Warning: High
```

2. **Standard Curve Information**:
```
Curve Type: 4PL
Equation: y = 3.25 + (0.001-3.25)/(1 + (x/15.2)^1.12)
R² = 0.9998
Working Range: 1-1000 ng/mL
```

3. **Statistical Analysis**:
- R² and adjusted R²
- Residual analysis
- Confidence intervals
- Prediction intervals

## Applications

Perfect for concentration calculations in:
- Protein assays (BCA, Bradford)
- Immunoassays (ELISA)
- Chemical analysis
- Biomarker quantification
- Quality control testing

## Mathematical Methods

### Linear Regression
- y = mx + b
- Used for linear standard curves
- Direct concentration calculation

### 4PL (Four Parameter Logistic)
- y = d + (a-d)/(1 + (x/c)^b)
- a: Maximum asymptote
- b: Slope factor
- c: EC50 (mid-point)
- d: Minimum asymptote

### 5PL (Five Parameter Logistic)
- Adds asymmetry parameter
- Better for asymmetric curves
- More complex but sometimes necessary
