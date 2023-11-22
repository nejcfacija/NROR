#include <iostream>
#include <cmath>

double f(double x) {
    return std::pow(std::exp(1.0), 3*x) * std::pow(std::atan(x), 2);
}

double arctan_taylor(double x, int terms) {
    double result = 0;
    for (int n = 0; n < terms; ++n) {
        double term = std::pow(-1, n) * std::pow(x, 2*n + 1) / (2*n + 1);
        result += term;
    }
    return result;
}

double trapezoidal_rule(double a, double b, int n) {
    double h = (b - a) / n;
    double result = 0.5 * h * (f(a) + f(b));
    for (int i = 1; i < n; ++i) {
        result += h * f(a + i*h);
    }
    return result;
}

int main() {
    double a = 0;
    double b = M_PI / 4;
    int n = 4;

    double integral_approximation = trapezoidal_rule(a, b, n);
    std::cout << "Približna vrednost integrala: " << integral_approximation << std::endl;

    return 0;
}
