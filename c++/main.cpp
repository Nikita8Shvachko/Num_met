#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

double function_1(double x) {//y= |_frac_{{1};{({x^{5}+x^{4}+1})}}-log({x})
    return 1 / (1 + (pow(x, 5) + pow(x, 4) + 1)) - log(x) / log(10);
}

double function_2(double x) {//y=x^{5}+0.1x^{4}+0.4abs({x})-1.2
    return pow(x, 5) + 0.1 * pow(x, 4) + 0.4 * abs(x) - 1.2;
}

// Функция для вычисления значений функции в узлах
vector<double> computeFunctionValues(const vector<double> &nodes, double (*func)(double)) {
    vector<double> values;

    for (double x: nodes) {
        values.push_back(func(x));
    }

    return values;
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

// Функция для вычисления значений полинома в узлах и серединах между узлами
vector<double> computePolynomialValues(const vector<double> &x_nodes, const vector<double> &y_values, int num_nodes) {
    vector<double> polynomial_values;

    // Проходим по узлам интерполяции
    for (int i = 0; i < num_nodes; ++i) {
        // Вычисляем значение полинома в текущем узле
        double y_interp = poly_newton_method_left_right(x_nodes, y_values, x_nodes[i]);
        polynomial_values.push_back(y_interp);

        // Если это не последний узел, вычисляем значение полинома в середине между текущим и следующим узлами
        if (i < num_nodes - 1) {
            double x_mid = 0.5 * (x_nodes[i] + x_nodes[i + 1]);
            y_interp = poly_newton_method_left_right(x_nodes, y_values, x_mid);
            polynomial_values.push_back(y_interp);
        }
    }

    return polynomial_values;
}

// Функция для выбора равномерно распределенных узлов в заданном интервале
vector<double> chooseNodes(double start, double end, int num_nodes) {
    vector<double> nodes;

    double step = (end - start) / (num_nodes - 1);
    for (int i = 0; i < num_nodes; ++i) {
        nodes.push_back(start + i * step);
    }

    return nodes;
}

// Функция для вычисления фактической ошибки в заданной точке
double computeActualError(double x, double y_true, double y_interp) {
    return abs(y_true - y_interp); // Просто возвращает модуль разности между значением функции и полинома
}

vector<double> computeInterpolationError(const vector<double> &x_nodes, const vector<double> &y_values, int num_nodes,
                                         double (*func)(double)) {
    vector<double> interpolation_error;

    // Проходим по узлам интерполяции
    for (int i = 0; i < num_nodes; ++i) {
        // Вычисляем значение полинома в текущем узле
        double y_interp = poly_newton_method_left_right(x_nodes, y_values, x_nodes[i]);
        double error = abs(y_values[i] - y_interp);
        interpolation_error.push_back(error);

        // Если это не последний узел, вычисляем значение полинома в середине между текущим и следующим узлами
        if (i < num_nodes - 1) {
            double x_mid = 0.5 * (x_nodes[i] + x_nodes[i + 1]);
            y_interp = poly_newton_method_left_right(x_nodes, y_values, x_mid);
            error = abs(func(x_mid) - y_interp);
            interpolation_error.push_back(error);
        }
    }

    return interpolation_error;
}

void generateAndSaveDataset(const string &filename, double start, double end, int numPoints,
                            const vector<double> &(*gridFunction)(double, double, int)) {
    ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        cerr << "Error: Unable to open file " << filename << " for writing." << endl;
        return;
    }

    // Generate points using the function
    vector<double> xValues = gridFunction(start, end, numPoints);
    for (double x: xValues) {
        double y = function_1(x); // Change this to function_2 if needed
        outputFile << x << " " << y << endl;
    }

    outputFile.close();
}

int main() {
// Интервал непрерывности для каждой функции
    double start_function_1 = 0.1, end_function_1 = 10;
    double start_function_2 = -10, end_function_2 = 10;

    // Открываем файлы для записи результатов
    ofstream outfile_function_1("results_function_1.txt");
    ofstream outfile_function_2("results_function_2.txt");

    if (!outfile_function_1.is_open() || !outfile_function_2.is_open()) {
        cerr << "Error: Unable to open files for writing." << endl;
        return 1;
    }

    // Перебираем числа узлов от 3 до 10
    for (int num_nodes = 3; num_nodes <= 10; ++num_nodes) {
        // Выбираем узлы интерполяции для функции 1
        vector<double> nodes_function_1 = chooseNodes(start_function_1, end_function_1, num_nodes);

        // Вычисляем значения функции в узлах для функции 1
        vector<double> values_function_1 = computeFunctionValues(nodes_function_1, function_1);

        // Вычисляем значения полинома и фактические ошибки для функции 1
        vector<double> polynomial_values_function_1 = computePolynomialValues(nodes_function_1, values_function_1,
                                                                              num_nodes);

        // Выводим результаты для функции 1 в файл
        outfile_function_1 << "Number of nodes: " << num_nodes << endl;
        outfile_function_1 << "Nodes for function 1: ";
        for (double node: nodes_function_1) {
            outfile_function_1 << node << " ";
        }
        outfile_function_1 << endl;

        outfile_function_1 << "Values of polynomial for function 1: ";
        for (double value: polynomial_values_function_1) {
            outfile_function_1 << value << " ";
        }
        outfile_function_1 << endl;

        // Перебираем узлы для функции 2, вычисляем значения функции и т.д.
        // То же самое нужно проделать для функции 2

        // Для примера я закомментирую код для функции 2

        // Выбираем узлы интерполяции для функции 2
        vector<double> nodes_function_2 = chooseNodes(start_function_2, end_function_2, num_nodes);

        // Вычисляем значения функции в узлах для функции 2
        vector<double> values_function_2 = computeFunctionValues(nodes_function_2, function_2);

        // Вычисляем значения полинома и фактические ошибки для функции 2
        vector<double> polynomial_values_function_2 = computePolynomialValues(nodes_function_2, values_function_2,
                                                                              num_nodes);

        // Выводим результаты для функции 2 в файл
        outfile_function_2 << "Number of nodes: " << num_nodes << endl;
        outfile_function_2 << "Nodes for function 2: ";
        for (double node: nodes_function_2) {
            outfile_function_2 << node << " ";
        }
        outfile_function_2 << endl;

        outfile_function_2 << "Values of polynomial for function 2: ";
        for (double value: polynomial_values_function_2) {
            outfile_function_2 << value << " ";
        }
        outfile_function_2 << endl;

        cout << endl;
    }

    // Закрываем файлы
    outfile_function_1.close();
    outfile_function_2.close();

    return 0;
}

