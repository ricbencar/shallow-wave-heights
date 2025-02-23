# Shallow Wave Heights Program

## Program Overview

This program computes local shallow-foreshore wave-height distribution parameters using a Composed Weibull distribution model, as described in the thesis by H. Groenendijk, Master's Thesis, Delft University of Technology, 1998. 

### Key Steps:

1. **Input Parameters**: Reads the local significant spectral wave height (\(H_{m0}\), in meters) and the local depth (\(d\), in meters) from command line or user prompt.
2. **Free-Surface Variance Calculation**: Computes the free-surface variance \(m_0\) using the deep-water assumption:
   \[
   m_0 = \left(\frac{H_{m0}}{4}\right)^2
   \]
3. **Mean Square Wave Height Calculation**: Computes the mean square wave height as:
   \[
   H_{rms} = 3 \sqrt{m_0}
   \]
4. **Dimensional and Dimensionless Transitional Wave Height**: Computes the dimensional transitional wave height:
   \[
   H_{tr} = 0.12 \frac{d}{\sqrt{m_0}}
   \]
   and then calculates the dimensionless transitional parameter:
   \[
   \tilde{H}_{tr} = \frac{H_{tr}}{H_{rms}}
   \]
5. **Natural Cubic Spline Interpolation**: Interpolates values using Groenendijk’s Table 7.1, which comprises various dimensionless wave heights as ratios relative to \(H_{rms}\).
6. **Final Wave Heights Calculation**: Computes final (dimensional) wave heights by multiplying dimensionless values by \(H_{rms}\).
7. **Report Generation**: Produces a detailed report that includes:
   - Input parameters (\(H_{m0}\) and \(d\))
   - Calculated parameters (\(m_0\), \(H_{rms}\), and \(\tilde{H}_{tr}\))
   - Dimensionless wave heights
   - Final wave heights (in meters)
   - Ratios among wave heights

The report is printed to the console and also written to `report.txt`.

## Compilation Instructions

To compile the program, use the following command:
```bash
g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ -o shallow-wave-heights shallow-wave-heights.cpp
