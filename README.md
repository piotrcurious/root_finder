# root_finder
Simple root finder illustrating how polynomial evaluation can be abstracted from decimal number base.
Dreamed by BingAI aka copilot

This is a very interesting and challenging problem. I will try to generate a possible solution using my own knowledge and creativity.

First, let me explain the idea of using analytical cyclotomic polynomials to approximate polynomial roots. A cyclotomic polynomial is a polynomial that has as its roots the primitive nth roots of unity, i.e., the complex numbers of the form $$e^{2\pi i k/n}$$, where $$k$$ and $$n$$ are coprime integers. For example, the third cyclotomic polynomial is $$\Phi_3(x) = x^2 + x + 1$$, which has as its roots the complex cube roots of unity, $$e^{2\pi i/3}$$ and $$e^{-2\pi i/3}$$.

The advantage of using cyclotomic polynomials is that they are irreducible over the rational numbers, meaning that they cannot be factored into lower degree polynomials with rational coefficients. This implies that they are also irreducible over any algebraic extension of the rational numbers, such as the field of algebraic numbers, which contains all the roots of any polynomial with rational coefficients. Therefore, any algebraic number can be expressed as a linear combination of cyclotomic polynomials with rational coefficients, i.e., $$\alpha = \sum_{n=1}^\infty a_n \Phi_n(x)$$, where $$a_n \in \mathbb{Q}$$ and only finitely many of them are nonzero.

The idea of using analytical cyclotomic polynomials to approximate polynomial roots is to truncate this infinite series at some finite order $$N$$, and use the resulting polynomial $$\alpha_N = \sum_{n=1}^N a_n \Phi_n(x)$$ as an approximation of the algebraic number $$\alpha$$. The error of this approximation can be bounded by the sum of the absolute values of the remaining coefficients, i.e., $$|\alpha - \alpha_N| \leq \sum_{n=N+1}^\infty |a_n|$$.

To find the coefficients $$a_n$$, we can use the fact that the cyclotomic polynomials are orthogonal with respect to the standard inner product on the space of polynomials, i.e., $$\langle \Phi_m, \Phi_n \rangle = \int_0^1 \Phi_m(x) \Phi_n(x) dx = \delta_{mn}$$, where $$\delta_{mn}$$ is the Kronecker delta. This means that we can use the formula $$a_n = \langle \alpha, \Phi_n \rangle = \int_0^1 \alpha(x) \Phi_n(x) dx$$ to compute the coefficients.

To apply this method to a polynomial of degree six, we need to find a way to represent its roots as algebraic numbers. One possible way is to use the Cardano's formula, which gives the roots of a cubic equation in terms of radicals. For example, if we have a cubic equation of the form $$x^3 + px + q = 0$$, then one of its roots is given by $$x = \sqrt[3]{-\frac{q}{2} + \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}} + \sqrt[3]{-\frac{q}{2} - \sqrt{\frac{q^2}{4} + \frac{p^3}{27}}}$$.

We can use this formula to find one root of a polynomial of degree six, by dividing it by a quadratic factor and reducing it to a cubic equation. For example, if we have a polynomial of the form $$x^6 + ax^4 + bx^2 + c = 0$$, then we can divide it by $$x^2 + d$$, where $$d$$ is any constant, and obtain a cubic equation of the form $$(x^2 + d)^3 + a(x^2 + d) + b = 0$$. Then we can apply the Cardano's formula to find one root of this cubic equation, and use it as an approximation of one root of the original polynomial.

To find the other roots, we can use the fact that the roots of a polynomial are symmetric functions of its coefficients, meaning that they can be expressed in terms of elementary symmetric polynomials. For example, if we have a polynomial of the form $$x^6 + ax^4 + bx^2 + c = 0$$, then its roots satisfy the equations $$\sum_{i=1}^6 x_i = 0$$, $$\sum_{i<j} x_i x_j = -a$$, $$\sum_{i<j<k} x_i x_j x_k = 0$$, and $$\prod_{i=1}^6 x_i = -c$$, where $$x_i$$ are the roots. Therefore, if we know one root $$x_1$$, we can use these equations to find the other roots in terms of $$x_1$$.

To illustrate this method, I will generate a C++ example of a practical implementation. The code will take as input the coefficients of a polynomial of degree six, and output the approximated roots using analytical cyclotomic polynomials of order 10. The code will also compute the exact roots using the `roots` function from the `sympy` library, and compare the errors of the approximations.

```cpp
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <sympy/sympy.h>

using namespace std;
using namespace sympy;

// Define the cyclotomic polynomials of order up to 10
vector<function<complex<double>(complex<double>)>> phi = {
    [](complex<double> x) { return 1.0; },
    [](complex<double> x) { return x - 1.0; },
    [](complex<double> x) { return x + 1.0; },
    [](complex<double> x) { return x*x + x + 1.0; },
    [](complex<double> x) { return x*x - 1.0; },
    [](complex<double> x) { return x*x - x + 1.0; },
    [](complex<double> x) { return x*x + x - 1.0; },
    [](complex<double> x) { return x*x*x + x*x + x + 1.0; },
    [](complex<double> x) { return x*x*x - 1.0; },
    [](complex<double> x) { return x*x*x - x + 1.0; },
    [](complex<double> x) { return x*x*x + x - 1.0; }
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
        result += a_n * phi[n](0.0); // evaluate the nth cyclotomic polynomial at zero
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
    complex<double> x3 = -x1 - a / 2.0 - sqrt(x1 * x1 + a * x1 + b /

Source: Conversation with Bing, 3/8/2024
(1) Polynomial Approximation, Interpolation, and Orthogonal Polynomials. https://ads.harvard.edu/books/1990fnmd.book/chapt3.pdf.
(2) Fast real and complex root-finding methods for well-conditioned polynomials. https://arxiv.org/pdf/2102.04180.pdf.
(3) Orthogonal Polynomials Approximation Algorithm: probability densities. https://arxiv.org/pdf/2211.08594v2.pdf.
(4) Analytical Method for Finding Polynomial Roots. https://www.m-hikari.com/ams/ams-2015/ams-93-96-2015/p/kamyshlovAMS93-96-2015.pdf.
(5) Find the Roots of a Polynomial Algebraically or Numerically. https://docs.sympy.org/latest/guides/solving/find-roots-polynomial.html.
