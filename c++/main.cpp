#include <iostream>
#include <cmath>

using namespace std;

double function_1(double x) {//y= |_frac_{{1};{({x^{5}+x^{4}+1})}}-log({x})
    return 1 / (1 + (pow(x, 5) + pow(x, 4) + 1)) - log(x) / log(10);
}

double function_2(double x) {//y=x^{5}+0.1x^{4}+0.4abs({x})-1.2
    return pow(x, 5) + 0.1 * pow(x, 4) + 0.4 * abs(x) - 1.2;
}

vector<double> uniform_grid(double start, double end, int numPoints) {
    vector<double> grid;
    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        grid.push_back(start + i * step);
    }
    return grid;
}

// Функция для создания сетки Чебышева
vector<double> chebyshev_grid(double start, double end, int numPoints) {
    vector<double> grid;
    for (int i = 0; i < numPoints; ++i) {
        double point = 0.5 * ((start - end) * cos(M_PI * (2.0 * i + 1) / (2.0 * numPoints)) + (start + end));
        grid.push_back(point);
    }
    return grid;
}

double poly_newton_method_left_right(const std::vector<double> &x, const std::vector<double> &y, double xi) {
    int n = x.size();
    double result = 0;

    std::vector<std::vector<double>> f(n, std::vector<double>(n));

    // Initialize the first column of the divided difference table
    for (int i = 0; i < n; ++i) {
        f[i][0] = y[i];
    }

    // Construct the divided difference table
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n - i; ++j) {
            f[j][i] = (f[j + 1][i - 1] - f[j][i - 1]) / (x[j + i] - x[j]);
        }
    }

    // Evaluate the interpolating polynomial
    for (int i = 0; i < n; ++i) {
        double term = f[0][i];
        for (int j = 0; j < i; ++j) {
            term *= (xi - x[j]);
        }
        result += term;
    }

    return result;
}

double poly_newton_method_right_left(const std::vector<double> &x, const std::vector<double> &y, double xi) {
    int n = x.size();
    double result = 0;

    std::vector<std::vector<double>> f(n, std::vector<double>(n));

    // Initialize the last column of the divided difference table
    for (int i = 0; i < n; ++i) {
        f[i][n - 1] = y[i];
    }

    // Construct the divided difference table
    for (int i = n - 2; i >= 0; --i) {
        for (int j = n - 1; j > i; --j) {
            f[j][i] = (f[j][i + 1] - f[j - 1][i + 1]) / (x[j] - x[j - i - 1]);
        }
    }

    // Evaluate the interpolating polynomial
    for (int i = 0; i < n; ++i) {
        double term = f[i][0];
        for (int j = 0; j < i; ++j) {
            term *= (xi - x[n - 1 - j]);
        }
        result += term;
    }

    return result;
}

// Function to generate dataset for plotting
vector<pair<double, double>> generateDataset(double start, double end, int numPoints, const vector<double> &(*gridFunction)(double, double, int)) {
    vector<pair<double, double>> dataset;

    // Generate points using the function
    vector<double> xValues = gridFunction(start, end, numPoints);
    for (double x: xValues) {
        double y = function_1(x); // Change this to function_2 if needed
        dataset.push_back({x, y});
    }

    return dataset;
}

int main() {


}