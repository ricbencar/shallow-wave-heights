# Shallow Wave Heights calculator

## Overview
This program computes local shallow-foreshore wave-height distribution parameters using a **Composed Weibull distribution model** based on the work of **H. Groenendijk** in:

> *"Shallow foreshore wave height statistics"* - Master's Thesis, Delft University of Technology, 1998.

The program is designed to calculate key wave parameters, perform cubic spline interpolation on empirical data, and generate a comprehensive report.

## Features
- Reads **local significant spectral wave height (Hm0)** and **local depth (d)** from the command line or user input.
- Computes key wave parameters including:
  - Free-surface variance `m0`
  - Mean square wave height `Hrms`
  - Dimensional and dimensionless transitional wave heights (`Htr`, `H̃_tr`)
- Uses **natural cubic spline interpolation** on empirical data from *Groenendijk's Table 7.1*.
- Computes final **dimensional wave heights** based on interpolated values.
- Generates a **detailed report** with all calculations and results, saved to `report.txt`.

## Installation
To compile and run the program, follow these steps:

### Prerequisites
- A C++ compiler supporting **C++17** (e.g., GCC, Clang, MSVC)
- OpenMP for parallel execution (optional but recommended)

### Compilation
Run the following command to compile the program:
```sh
 g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
     -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
     -o shallow-wave-heights shallow-wave-heights.cpp
```

### Running the Program
The program accepts input via **command-line arguments** or **interactive input**:

#### Command-line usage:
```sh
./shallow-wave-heights <Hm0> <d>
```
Example:
```sh
./shallow-wave-heights 1.5 5.0
```

#### Interactive mode:
If no command-line arguments are provided, the program prompts the user for input:
```sh
Enter local significant spectral wave height Hm0 (m): 1.5
Enter local depth d (m): 5.0
```

## Computation Details
### 1. Input Parameters
The program requires two inputs:
- **Hm0**: Local significant spectral wave height (m)
- **d**: Local depth (m)

### 2. Calculated Parameters
The program calculates the following values:

#### **Free-surface variance**:
```math
m0 = \left(\frac{Hm0}{4}\right)^2
```
#### **Mean square wave height**:
```math
Hrms = 3 \times \sqrt{m0}
```
#### **Dimensional transitional wave height**:
```math
Htr = \frac{0.12 \times d}{\sqrt{m0}}
```
#### **Dimensionless transitional parameter**:
```math
H̃_tr = \frac{Htr}{Hrms}
```

### 3. Wave Height Interpolation
The program uses **natural cubic spline interpolation** to determine **nondimensional wave heights** from empirical data:

**Interpolated wave height ratios (H/Hrms):**
- **H̃₁**: Single highest wave
- **H̃₂**: Second highest wave
- **H̃(1/3)**: Significant wave height (1/3 highest waves)
- **H̃(1/10)**: Wave height from top 10%
- **H̃(1/50), H̃(1/100), H̃(1/1000)**: Extreme wave heights

### 4. Final Dimensional Wave Heights
The final wave heights are computed using:
```math
H_{final} = (H/Hrms) \times Hrms
```
for each of the interpolated values.

## Output
The program generates a **detailed report** that includes:
1. **Input Parameters**
2. **Computed Wave Parameters**
3. **Dimensionless Wave Heights (ratios H/Hrms)**
4. **Final Dimensional Wave Heights (in meters)**
5. **Diagnostic Ratios Among Wave Heights**

## References
- **H. Groenendijk**, *"Shallow foreshore wave height statistics"*, Master's Thesis, Delft University of Technology, 1998. Available at: https://repository.tudelft.nl/record/uuid:fe03dda9-40d9-4046-87fb-459f01fcd3d3
  
