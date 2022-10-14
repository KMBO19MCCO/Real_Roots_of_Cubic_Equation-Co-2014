#include <iostream>
#include <vector>
#include <cmath>
#include "excerpt.h"

using namespace std;
typedef float fp_t;
//const double PI = acos(-1.0);


template<typename fp_t>
void solve_cubic(std::vector<fp_t> coefficients, std::vector<fp_t> &roots) {
    fp_t a, b, c, d, phi, A, B, C;
    //Coefficients
    a = coefficients[3];
    b = coefficients[2];
    c = coefficients[1];
    d = coefficients[0];
//    cout << endl << a << "*x^3 + " << b << "*x^2 + " << c << "*x + " << d << endl;
    fp_t p = (3 * a * c - pow(b, 2)) / (3 * pow(a, 2));
    fp_t q = (2 * pow(b, 3) - 9 * a * b * c + 27 * pow(a, 2) * d) / (27 * pow(a, 3));

    A = 2 * sqrt(-p / 3);

    fp_t acos_arg = (3 * q) / (A * p);
    phi = acos(acos_arg) / 3; // TODO fix nan
    fp_t twothirdspi = 2 * std::numbers::pi / 3;


    if (p < 0) {
        A = 2 * sqrt(-p / 3);
        B = -b / (3 * a);
        C = 3 * q / (A * p);

        if (abs(C) <= 1) {
            // phi = acos(C)/3;
            roots[0] = A * cos(phi + 0 * twothirdspi) + B;
            roots[1] = A * cos(phi + 1 * twothirdspi) + B;
            roots[2] = A * cos(phi + 2 * twothirdspi) + B;
        } else if ((C < -1) || (C > 1)) {
//            cout << "Complex roots" << endl;
        }

//        cout << "p = " << p << " " << "q = " << q << endl;
//        cout << "Roots: " << roots[0] << "; " << roots[1] << "; " << roots[2] << endl;
    } else if (p > 0) {
//        cout << "One real root ";
        A = 2 * sqrt(p / 3);
        fp_t asinh_arg = (3 * q) / (A * p);
        phi = asinh(asinh_arg);

        roots[0] = -A * sinh(phi / 3);
//        cout << "p = " << p << " " << "q = " << q << endl;
//        cout << "Root: " << roots[0] << endl;
    }
}

template<typename fp_t>
auto testPolynomial(unsigned int roots_count) {
    fp_t deviation;
    vector<fp_t> roots_computed(roots_count);
    vector<fp_t> roots(roots_count), coefficients(roots_count + 1);
    generate_polynomial<fp_t>(roots_count, 0, roots_count, 0, 1e-5, -1, 1, roots, coefficients);
    solve_cubic(coefficients, roots_computed);
    auto result = compare_roots<fp_t>(roots_computed.size(), roots.size(), roots_computed, roots, deviation);
    switch (result) {
        case PR_2_INFINITE_ROOTS:
            cout << "INFINITE ROOTS";
            break;
        case PR_AT_LEAST_ONE_ROOT_IS_FAKE:
            cout << "AT LEAST ONE ROOT IS FAKE";
            break;
        case PR_AT_LEAST_ONE_ROOT_LOST:
            cout << "AT LEAST ONE ROOT LOST";
            break;
        default:
            break;
    }
    return deviation;
}

int main() {
//    vector<float> koef1 = {1, -10, 7, 1};
//    vector<float> koef2 = {1, 0, 4, 16};
    //vector<float> koef3 = { 1, -2,-1, 1 };
    //vector<float> koef4 = { 2, 4, 6, -5 };
//    vector<float> r1(3);
//    vector<float> r2(3);
    //vector<float> r3(3);
    //vector<float> r4(3);
    //{ 3, -2, 3, -2 }

    //3 different roots
//    solve_cubic(koef1, r1);
    //1 root of multiplicity 2;
//    solve_cubic(koef2, r2);
    //complex roots
    //solve_cubic(koef4, r4);
    float deviation, max_deviation = 0;
    for (auto i = 0; i < 100'000; ++i) {
        deviation = testPolynomial<float>(3);
        if (deviation > max_deviation) {
            max_deviation = deviation;
        }
    }
    cout << "Max deviation: " << max_deviation << endl;
}
