#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <sympy/sympy.h>

using namespace std;
using namespace sympy;

// Define the cyclotomic polynomials of order up to 10
vector<function<complex<double>(complex<double>)>> phi = {
     { return 1.0; },
     { return x - 1.0; },
     { return x + 1.0; },
     { return x*x + x + 1.0; },
     { return x*x - 1.0; },
     { return x*x - x + 1.0; },
     { return x*x + x - 1.0; },
     { return x*x*x + x*x + x + 1.0; },
     { return x*x*x - 1.0; },
     { return x*x*x - x + 1.0; },
     { return x*x*x + x - 1.0; }
};

// Define the inner product of two polynomials
complex<double> inner_product(function<complex<double>(complex<double>)> f, function<complex<double>(complex<double>)> g) {
    complex<double> result = 0.0;
    int n = 1000; // number of integration points
    double h = 1.0 / n; // integration step
    for (int i = 0; i < n; i++) {
        double x = i * h; // integration point
        result += f(x) * g(x) * h; // trapezoidal rule
    }
    return result;
}

// Define the function to approximate a polynomial root using analytical cyclotomic polynomials
complex<double> approximate_root(function<complex<double>(complex<double>)> f, int N) {
    complex<double> result = 0.0;
    for (int n = 1; n <= N; n++) {
        complex<double> a_n = inner_product(f, phi[n]); // coefficient of the nth cyclotomic polynomial
        result += a_n * phin; // evaluate the nth cyclotomic polynomial at zero
    }
    return result;
}

// Define the function to find one root of a polynomial of degree six using Cardano's formula
complex<double> cardano_root(double a, double b, double c, double d) {
    complex<double> p = a + 3.0 * d * d; // coefficient of the reduced cubic equation
    complex<double> q = b + 2.0 * d * (a + d * d); // constant term of the reduced cubic equation
    complex<double> r = -q / 2.0 + sqrt(q * q / 4.0 + p * p * p / 27.0); // first term of the Cardano's formula
    complex<double> s = -q / 2.0 - sqrt(q * q / 4.0 + p * p * p / 27.0); // second term of the Cardano's formula
    complex<double> t = pow(r, 1.0 / 3.0) + pow(s, 1.0 / 3.0); // sum of the cube roots
    complex<double> x = t - d; // one root of the original polynomial
    return x;
}

// Define the function to find the other roots of a polynomial of degree six using elementary symmetric polynomials
vector<complex<double>> other_roots(double a, double b, double c, complex<double> x1) {
    vector<complex<double>> result;
    complex<double> x2 = -x1 - a / 2.0 + sqrt(x1 * x1 + a * x1 + b / 2.0); // second root
    complex<double> x3 = -x1 - a / 2.0 - sqrt(x1 * x1 + a * x1 + b / 2.0); // third root
    complex<double> x4 = -x1 - x2 - x3; // fourth root
    complex<double> x5 = sqrt(-c / (x1 * x2 * x3 * x4)); // fifth root
    complex<double> x6 = -x5; // sixth root
    result.push_back(x2);
    result.push_back(x3);
    result.push_back(x4);
    result.push_back(x5);
    result.push_back(x6);
    return result;
}

// Define the main function to test the method
int main() {
    // Define the coefficients of a polynomial of degree six
    double a = 1.0;
    double b = -2.0;
    double c = -3.0;
    double d = 1.0;

    // Define the polynomial function
    function<complex<double>(complex<double>)> f = a, b, c {
        return x*x*x*x*x*x + a*x*x*x*x + b*x*x + c;
    };

    // Find one root using Cardano's formula
    complex<double> x1 = cardano_root(a, b, c, d);

    // Find the other roots using elementary symmetric polynomials
    vector<complex<double>> x2_6 = other_roots(a, b, c, x1);

    // Approximate the roots using analytical cyclotomic polynomials of order 10
    complex<double> y1 = approximate_root(f, 10);
    vector<complex<double>> y2_6;
    for (int i = 0; i < 5; i++) {
        y2_6.push_back(approximate_root([f, x2_6[i]](complex<double> x) {
            return f(x) / (x - x2_6[i]);
        }, 10));
    }

    // Compute the exact roots using sympy
    vector<complex<double>> z1_6;
    Symbol x("x");
    Expr p = pow(x, 6) + a * pow(x, 4) + b * pow(x, 2) + c;
    auto r = roots(p, x);
    for (auto& kv : r) {
        z1_6.push_back(complex<double>(N(kv.first).real(), N(kv.first).imag()));
    }

    // Print the results and compare the errors
    cout << "One root using Cardano's formula: " << x1 << endl;
    cout << "Other roots using elementary symmetric polynomials: " << endl;
    for (int i = 0; i < 5; i++) {
        cout << x2_6[i] << endl;
    }
    cout << "Approximated roots using analytical cyclotomic polynomials of order 10: " << endl;
    cout << y1 << endl;
    for (int i = 0; i < 5; i++) {
        cout << y2_6[i] << endl;
    }
    cout << "Exact roots using sympy: " << endl;
    for (int i = 0; i < 6; i++) {
        cout << z1_6[i] << endl;
    }
    cout << "Errors of the approximations: " << endl;
    cout << abs(y1 - z1_6[0]) << endl;
    for (int i = 0; i < 5; i++) {
        cout << abs(y2_6[i] - z1_6[i + 1]) << endl;
    }
    return 0;
}
