/***********************************************************************
 * Program: shallow-water-waves.cpp
 * 
 * Description:
 *   This program computes local shallow-foreshore wave-height distribution 
 *   parameters using a Composed Weibull distribution model as described in:
 * 
 *      "Shallow foreshore wave height statistics" by H. Groenendijk,
 *      Master's Thesis, Delft University of Technology, 1998.
 * 
 *   The program performs the following steps:
 *     1. Reads the local significant spectral wave height (Hm0, in m) and the
 *        local depth (d, in m) from the command line or by prompting the user.
 *     2. Computes the free-surface variance m0 using the deep-water assumption:
 *            m0 = (Hm0 / 4)^2
 *        (Note: For shallow water, a more appropriate local model for m0 may be needed.)
 *     3. Computes the mean square wave height as:
 *            Hrms = 3 * sqrt(m0)
 *     4. Computes the dimensional transitional wave height:
 *            Htr = 0.12 * d / sqrt(m0)
 *        and then calculates the dimensionless transitional parameter as:
 *            H̃_tr = Htr / Hrms
 *     5. Uses natural cubic spline interpolation on Groenendijk’s Table 7.1,
 *        where each column (H̃₁, H̃₂, H̃₁/₃, H̃₁/₁₀, H̃₁/₅₀, H̃₁/₁₀₀, H̃₁/₁₀₀₀)
 *        is expressed as a ratio relative to Hrms.
 *     6. Computes the final (dimensional) wave heights by multiplying the 
 *        nondimensional values by Hrms.
 *     7. Builds a detailed report that includes:
 *            - The input parameters (Hm0 and d)
 *            - Calculated parameters (m0, Hrms, and H̃_tr)
 *            - The dimensionless wave heights (ratios H/Hrms)
 *            - The final wave heights in meters
 *            - Ratios among wave heights
 *        The report is printed to the command line and written to "report.txt".
 *
 * Compile:
 *   g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic \
 *       -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++ \
 *       -o shallow-wave-heights shallow-wave-heights.cpp
 *
 * References:
 *   H. Groenendijk, "Shallow foreshore wave height statistics", Master's Thesis,
 *   Delft University of Technology, 1998.
 *   Available at: http://example.com/groenendijk1998_shallow_foreshore_wave_height_statistics.pdf
 *
 ***********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

using namespace std;

// ---------------------------------------------------------------------
// Function: computeSplineSecondDerivatives
//
// Purpose:
//   Given arrays x and y (data points), this routine computes the second 
//   derivatives m[] for natural cubic spline interpolation. Natural 
//   boundary conditions (m[0]=m[n-1]=0) are applied.
//
// Parameters:
//   x - vector of independent variable values (must be sorted in ascending order)
//   y - vector of dependent variable values corresponding to x
//
// Returns:
//   A vector<double> containing the second derivatives at each data point.
//
// References:
//   Standard texts on numerical analysis for cubic spline interpolation.
// ---------------------------------------------------------------------
vector<double> computeSplineSecondDerivatives(const vector<double>& x, const vector<double>& y) {
    size_t n = x.size();
    vector<double> m(n, 0.0);   // second derivatives
    vector<double> l(n, 0.0), mu(n, 0.0), z(n, 0.0);

    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    vector<double> h(n - 1, 0.0);
    for (size_t i = 0; i < n - 1; i++) {
        h[i] = x[i + 1] - x[i];
    }

    vector<double> alpha(n, 0.0);
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i])
                 - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    for (size_t i = 1; i < n - 1; i++) {
        l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    m[n - 1] = 0.0;

    // Back-substitution to compute second derivatives.
    for (std::vector<double>::difference_type j = static_cast<std::vector<double>::difference_type>(n) - 2; j >= 0; j--) {
        size_t jj = static_cast<size_t>(j);
        m[jj] = z[jj] - mu[jj] * m[jj + 1];
    }

    return m;
}

// ---------------------------------------------------------------------
// Function: cubicSplineInterpolation
//
// Purpose:
//   Given vectors x and y (data points) and a query point x0, this function 
//   computes the interpolated value y(x0) using natural cubic spline interpolation.
//
// Parameters:
//   x  - vector of independent variable values
//   y  - vector of dependent variable values corresponding to x
//   x0 - the query point where interpolation is required
//
// Returns:
//   The interpolated value y0 at x0.
//
// ---------------------------------------------------------------------
double cubicSplineInterpolation(const vector<double>& x, const vector<double>& y, double x0) {
    size_t n = x.size();
    if (n == 0) throw runtime_error("No data points provided in x.");
    if (n == 1) return y[0];

    // Linear extrapolation if x0 is outside the x-range.
    if (x0 <= x.front()) {
        double t = (x0 - x.front()) / (x[1] - x.front());
        return y.front() + t * (y[1] - y.front());
    }
    if (x0 >= x.back()) {
        double t = (x0 - x[n - 2]) / (x.back() - x[n - 2]);
        return y[n - 2] + t * (y.back() - y[n - 2]);
    }

    // Find the interval [x[i], x[i+1]] that contains x0.
    size_t i = 0;
    while (i < n - 1 && x0 > x[i + 1]) {
        i++;
    }

    double h = x[i + 1] - x[i];
    double A = (x[i + 1] - x0) / h;
    double B = (x0 - x[i]) / h;

    // Compute second derivatives using natural cubic spline.
    vector<double> m = computeSplineSecondDerivatives(x, y);

    // Cubic spline interpolation formula.
    double y0 = A * y[i] + B * y[i + 1]
              + ((A * A * A - A) * m[i] + (B * B * B - B) * m[i + 1]) * (h * h) / 6.0;

    return y0;
}

// ---------------------------------------------------------------------
// Main Program
//
// This program computes a local shallow-foreshore wave-height distribution based
// on the Composed Weibull distribution model as described in:
// 
//    "Shallow foreshore wave height statistics" by H. Groenendijk,
//    Master's Thesis, Delft University of Technology, 1998.
//    Available at: http://example.com/groenendijk1998_shallow_foreshore_wave_height_statistics.pdf
//
// The program calculates:
//   - Free surface variance m0 using the deep-water relation: m0 = (Hm0/4)^2
//   - Mean square wave height: Hrms = 3 * sqrt(m0)
//   - Dimensional transitional wave height: Htr = 0.12 * d / sqrt(m0)
//   - Dimensionless transitional parameter: H̃_tr = Htr / Hrms
//
// It then uses natural cubic spline interpolation on Groenendijk's Table 7.1, where
// each column (H̃₁, H̃₂, H̃₁/₃, H̃₁/₁₀, H̃₁/₅₀, H̃₁/₁₀₀, H̃₁/₁₀₀₀) is given as a ratio
// relative to Hrms. Final (dimensional) wave heights are computed as:
//      H_final = (nondimensional value) * Hrms
//
// The program builds a detailed report that includes:
//   - The input parameters (Hm0 and d)
//   - Calculated parameters (m0, Hrms, and H̃_tr)
//   - Dimensionless wave heights (ratios H/Hrms)
//   - Final wave heights in meters
//   - Various ratios among the wave heights
//
// The report is printed to the command line and saved to "report.txt".
// ---------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // ------------------------------
    // Input Section
    // ------------------------------
    double Hm0 = 0.0; // Local significant spectral wave height (m)
    double d   = 0.0; // Local depth (m)

    if (argc >= 3) {
        Hm0 = atof(argv[1]);
        d = atof(argv[2]);
    } else {
        cout << "Enter local significant spectral wave height Hm0 (m): ";
        cin >> Hm0;
        cout << "Enter local depth d (m): ";
        cin >> d;
    }

    if (Hm0 <= 0.0 || d <= 0.0) {
        cerr << "Error: Both Hm0 and d must be positive numbers." << endl;
        return 1;
    }

    // ------------------------------
    // Parameter Calculation Section
    // ------------------------------
    // Compute free-surface variance m0 (using deep-water approximation)
    double m0 = pow(Hm0 / 4.0, 2.0);

    // Compute mean square wave height Hrms = 3 * sqrt(m0)
    double Hrms = 3.0 * sqrt(m0);

    // Compute dimensional transitional wave height Htr = 0.12 * d / sqrt(m0)
    double Htr_dim = 0.12 * d / sqrt(m0);
    // Dimensionless transitional parameter relative to Hrms:
    double Htr_tilde = (Hrms > 0.0) ? (Htr_dim / Hrms) : 0.0;

    // ------------------------------
    // Table Data Section (Groenendijk's Table 7.1)
    //
    // tableX: dimensionless transitional wave height H̃_tr values (Htr/Hrms) from 0.05 to 3.50.
    // The table columns contain characteristic nondimensional wave heights relative to Hrms.
    // (For example, col1 gives H̃1 = H1/Hrms.)
    // ------------------------------
    vector<double> tableX = {
         0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
         0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
         1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50,
         1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00,
         2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50,
         2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00,
         3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50
    };

    // Corrected nondimensional data from Table 7.1 (each column is H/Hrms)
    // Column indices:
    // 0: H̃1, 1: H̃2, 2: H̃(1/3), 3: H̃(1/10), 4: H̃(1/50), 5: H̃(1/100), 6: H̃(1/1000)
    vector<double> col1 = {
      9.949, 5.916, 4.365, 3.518, 2.976, 2.595, 2.312, 2.092, 1.916, 1.772,
      1.651, 1.549, 1.462, 1.387, 1.322, 1.265, 1.216, 1.174, 1.137, 1.106,
      1.079, 1.056, 1.037, 1.021, 1.008, 0.998, 0.989, 0.983, 0.978, 0.974,
      0.972, 0.970, 0.970, 0.969, 0.970, 0.971, 0.972, 0.974, 0.975, 0.977,
      0.979, 0.981, 0.983, 0.985, 0.986, 0.988, 0.989, 0.991, 0.992, 0.993,
      0.994, 0.995, 0.996, 0.997, 0.997, 0.998, 0.998, 0.999, 0.999, 0.999,
      0.999, 0.999, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000
    };
    vector<double> col2 = {
      1.029, 1.029, 1.029, 1.029, 1.029, 1.029, 1.030, 1.030, 1.030, 1.030,
      1.031, 1.032, 1.033, 1.035, 1.037, 1.040, 1.043, 1.048, 1.053, 1.059,
      1.066, 1.075, 1.084, 1.094, 1.105, 1.117, 1.130, 1.144, 1.158, 1.172,
      1.187, 1.202, 1.218, 1.233, 1.249, 1.265, 1.281, 1.297, 1.312, 1.328,
      1.344, 1.359, 1.374, 1.390, 1.404, 1.419, 1.433, 1.448, 1.462, 1.475,
      1.489, 1.502, 1.515, 1.528, 1.540, 1.553, 1.565, 1.577, 1.589, 1.601,
      1.612, 1.623, 1.635, 1.646, 1.657, 1.668, 1.679, 1.689, 1.700, 1.711
    };
    vector<double> col3 = {
      1.249, 1.249, 1.249, 1.249, 1.249, 1.249, 1.249, 1.250, 1.250, 1.250,
      1.251, 1.252, 1.253, 1.255, 1.258, 1.262, 1.266, 1.271, 1.278, 1.285,
      1.294, 1.304, 1.314, 1.321, 1.327, 1.331, 1.335, 1.338, 1.341, 1.344,
      1.347, 1.350, 1.354, 1.357, 1.361, 1.364, 1.368, 1.372, 1.375, 1.379,
      1.382, 1.386, 1.389, 1.392, 1.395, 1.397, 1.400, 1.402, 1.404, 1.406,
      1.407, 1.409, 1.410, 1.411, 1.412, 1.413, 1.413, 1.414, 1.414, 1.414,
      1.415, 1.415, 1.415, 1.415, 1.416, 1.416, 1.416, 1.416, 1.416, 1.416
    };
    vector<double> col4 = {
      1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.438, 1.439,
      1.440, 1.441, 1.442, 1.445, 1.448, 1.452, 1.457, 1.463, 1.470, 1.479,
      1.489, 1.501, 1.514, 1.528, 1.544, 1.561, 1.578, 1.597, 1.617, 1.637,
      1.654, 1.669, 1.682, 1.694, 1.704, 1.714, 1.722, 1.730, 1.737, 1.744,
      1.750, 1.756, 1.761, 1.766, 1.770, 1.774, 1.778, 1.781, 1.784, 1.786,
      1.789, 1.791, 1.792, 1.794, 1.795, 1.796, 1.797, 1.797, 1.798, 1.798,
      1.799, 1.799, 1.799, 1.799, 1.799, 1.800, 1.800, 1.800, 1.800, 1.800
    };
    vector<double> col5 = {
      1.520, 1.520, 1.520, 1.520, 1.520, 1.520, 1.521, 1.520, 1.521, 1.521,
      1.522, 1.523, 1.525, 1.528, 1.531, 1.535, 1.540, 1.547, 1.555, 1.564,
      1.575, 1.587, 1.601, 1.616, 1.632, 1.650, 1.669, 1.689, 1.709, 1.731,
      1.753, 1.775, 1.798, 1.821, 1.844, 1.868, 1.891, 1.915, 1.926, 1.929,
      1.933, 1.936, 1.940, 1.944, 1.947, 1.951, 1.954, 1.957, 1.960, 1.962,
      1.965, 1.967, 1.969, 1.970, 1.971, 1.973, 1.974, 1.974, 1.975, 1.976,
      1.976, 1.977, 1.977, 1.977, 1.977, 1.978, 1.978, 1.978, 1.978, 1.978
    };
    vector<double> col6 = {
      1.592, 1.592, 1.592, 1.592, 1.593, 1.593, 1.593, 1.593, 1.593, 1.594,
      1.595, 1.596, 1.598, 1.600, 1.604, 1.608, 1.614, 1.621, 1.629, 1.638,
      1.650, 1.663, 1.677, 1.693, 1.710, 1.729, 1.748, 1.769, 1.791, 1.813,
      1.836, 1.860, 1.884, 1.908, 1.932, 1.957, 1.981, 2.006, 2.030, 2.055,
      2.079, 2.103, 2.105, 2.109, 2.113, 2.116, 2.120, 2.123, 2.126, 2.129,
      2.132, 2.134, 2.136, 2.137, 2.139, 2.140, 2.141, 2.142, 2.143, 2.144,
      2.144, 2.144, 2.145, 2.145, 2.145, 2.146, 2.146, 2.146, 2.146, 2.146
    };
    vector<double> col7 = {
      1.788, 1.788, 1.788, 1.788, 1.788, 1.788, 1.788, 1.789, 1.789, 1.790,
      1.791, 1.792, 1.794, 1.797, 1.801, 1.806, 1.812, 1.820, 1.829, 1.840,
      1.852, 1.867, 1.883, 1.901, 1.920, 1.941, 1.963, 1.987, 2.011, 2.036,
      2.062, 2.088, 2.115, 2.142, 2.170, 2.197, 2.225, 2.252, 2.280, 2.307,
      2.334, 2.361, 2.388, 2.414, 2.440, 2.465, 2.490, 2.515, 2.539, 2.563,
      2.586, 2.609, 2.616, 2.618, 2.620, 2.621, 2.623, 2.624, 2.625, 2.625,
      2.626, 2.626, 2.627, 2.627, 2.627, 2.628, 2.628, 2.628, 2.628, 2.628
    };

    vector<vector<double>> tableY_all = {col1, col2, col3, col4, col5, col6, col7};

    // ------------------------------
    // Interpolation Section
    // ------------------------------
    // Interpolate each column at x = H̃_tr using natural cubic spline interpolation.
    vector<double> dimensionlessOutputs(7, 0.0);
    for (size_t c = 0; c < 7; c++) {
        dimensionlessOutputs[c] = cubicSplineInterpolation(tableX, tableY_all[c], Htr_tilde);
    }

    // Compute final (dimensional) wave heights (in m) by multiplying nondimensional outputs by Hrms.
    vector<double> finalHeights(7, 0.0);
    for (size_t i = 0; i < 7; i++) {
        finalHeights[i] = dimensionlessOutputs[i] * Hrms;
    }

    // For clarity, assign names:
    double H1_ratio     = dimensionlessOutputs[0];  // H̃1 = H1/Hrms
    double H2_ratio     = dimensionlessOutputs[1];  // H̃2 = H2/Hrms
    double H13_ratio    = dimensionlessOutputs[2];  // H̃(1/3) = (H1/3)/Hrms
    double H110_ratio   = dimensionlessOutputs[3];  // H̃(1/10)
    double H150_ratio   = dimensionlessOutputs[4];  // H̃(1/50)
    double H1100_ratio  = dimensionlessOutputs[5];  // H̃(1/100)
    double H11000_ratio = dimensionlessOutputs[6];  // H̃(1/1000)

    double H1_m = finalHeights[0];
    double H2_m = finalHeights[1];
    double H1_3_m = finalHeights[2];
    double H1_10_m = finalHeights[3];
    double H1_50_m = finalHeights[4];
    double H1_100_m = finalHeights[5];
    double H1_1000_m = finalHeights[6];

    // Compute some diagnostic ratios:
    double ratio_110_13   = (H1_3_m != 0.0) ? (H1_10_m / H1_3_m) : 0.0;
    double ratio_150_13   = (H1_3_m != 0.0) ? (H1_50_m / H1_3_m) : 0.0;
    double ratio_1100_13  = (H1_3_m != 0.0) ? (H1_100_m / H1_3_m) : 0.0;
    double ratio_11000_13 = (H1_3_m != 0.0) ? (H1_1000_m / H1_3_m) : 0.0;
    double ratio_11000_110= (H1_10_m != 0.0) ? (H1_1000_m / H1_10_m) : 0.0;

    // ------------------------------
    // Report Generation Section
    // ------------------------------
    ostringstream reportStream;
    reportStream << fixed << setprecision(4);
    reportStream << "-----------------------------------------------------------\n";
    reportStream << "Input Parameters:\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << setw(40) << left << "Local significant spectral wave height (Hm0) [m]"
                 << " = " << right << setw(10) << Hm0 << "\n";
    reportStream << setw(40) << left << "Local depth (d) [m]"
                 << " = " << right << setw(10) << d << "\n\n";
    
    reportStream << "-----------------------------------------------------------\n";
    reportStream << "Calculated Parameters:\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << setw(40) << left << "Free surface variance (m^2) m0"
                 << " = " << right << setw(10) << m0 << "\n";
    reportStream << setw(40) << left << "Mean square wave height (m) Hrms"
                 << " = " << right << setw(10) << Hrms << "\n";
    reportStream << setw(40) << left << "Compound Weibull transition (adim) H~_tr"
                 << " = " << right << setw(10) << Htr_tilde << "\n\n";

    reportStream << setprecision(3);
    reportStream << "-----------------------------------------------------------\n";
    reportStream << "Dimensionless Wave Heights (H/Hrms):\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << setw(8)  << "H~_tr"
                 << setw(8)  << "H~_1"
                 << setw(8)  << "H~_2"
                 << setw(10) << "H~_(1/3)"
                 << setw(10) << "H~_(1/10)"
                 << setw(10) << "H~_(1/50)"
                 << setw(10) << "H~_(1/100)"
                 << setw(12) << "H~_(1/1000)" << "\n";
    reportStream << setw(8)  << Htr_tilde
                 << setw(8)  << H1_ratio
                 << setw(8)  << H2_ratio
                 << setw(10) << H13_ratio
                 << setw(10) << H110_ratio
                 << setw(10) << H150_ratio
                 << setw(10) << H1100_ratio
                 << setw(12) << H11000_ratio << "\n\n";

    reportStream << "-----------------------------------------------------------\n";
    reportStream << "Dimensional Wave Heights (m):\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << setw(8)  << "H1(m)"
                 << setw(8)  << "H2(m)"
                 << setw(9)  << "H1/3(m)"
                 << setw(10) << "H1/10(m)"
                 << setw(10) << "H1/50(m)"
                 << setw(10) << "H1/100(m)"
                 << setw(12) << "H1/1000(m)" << "\n";
    reportStream << setw(8)  << H1_m
                 << setw(8)  << H2_m
                 << setw(9)  << H1_3_m
                 << setw(10) << H1_10_m
                 << setw(10) << H1_50_m
                 << setw(10) << H1_100_m
                 << setw(12) << H1_1000_m << "\n\n";

    reportStream << "-----------------------------------------------------------\n";
    reportStream << "Ratios:\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << "  (H1/10)/(H1/3)   = " << ratio_110_13 << "\n";
    reportStream << "  (H1/50)/(H1/3)   = " << ratio_150_13 << "\n";
    reportStream << "  (H1/100)/(H1/3)  = " << ratio_1100_13 << "\n";
    reportStream << "  (H1/1000)/(H1/3) = " << ratio_11000_13 << "\n";
    reportStream << "  (H1/1000)/(H1/10)= " << ratio_11000_110 << "\n\n";
    reportStream << "-----------------------------------------------------------\n";
    reportStream << "End of Report\n";

    // ------------------------------
    // Output Section: Print and Write Report
    // ------------------------------
    cout << reportStream.str();

    ofstream ofs("report.txt");
    if (!ofs) {
        cerr << "Error: Could not open report.txt for writing.\n";
        return 1;
    }
    ofs << reportStream.str();
    ofs.close();

    return 0;
}
