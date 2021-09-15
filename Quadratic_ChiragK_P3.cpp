/**
 * @author Chirag Kaushik
 * @date August 31, 2021
 * 
 * Various methods of testing floating point precision.
 */

#include <iostream>
#include <cmath>

#define ENDL "\n"


float logistic1(float x) {
    return std::expf(x) / (std::expf(x) - 1);
}

float logistic2(float x) {
    return 1 / (1 - std::expf(-x));
}

/**
 * Runs a quadratic formula test for precision, and uses multiple methods to calculate the roots.
 */
void quad() {
    float a,b,c;
    std::cout << "a: ";
    std::cin >> a;
    std::cout << "b: ";
    std::cin >> b;
    std::cout << "c: ";
    std::cin >> c;
    float r1 = (-b - std::sqrtf(b * b - 4.f * a * c)) / (2.f * a);
    float r2 = (-b + std::sqrtf(b * b - 4.f * a * c)) / (2.f * a);
    std::cout << "root 1: " << r1 << ENDL << "root 2: " << r2 << ENDL;
    float x2 = c / (a * r1);
    std::cout << "root 2 (using root 1): " << x2 << ENDL;

    float plug1 = a * r1 * r1 + b * r1 + c;
    float plug2 = a * r2 * r2 + b * r2 + c;
    float plug3 = a * x2 * x2 + b * x2 + c;

    std::cout << "plugged back in root 1: " << plug1 << ENDL;
    std::cout << "plugged back in root 2: " << plug2 << ENDL;
    std::cout << "plugged back in root 2 (using root 1): " << plug3 << ENDL;
}

// WIP
void logistic() {
    float x;
    std::cout << "x: ";
    std::cin >> x;

    float f1 = logistic1(x);
    float f2 = logistic2(x);

    std::cout << "f1: " << f1 << ENDL;
    std::cout << "f2: " << f2 << ENDL;
}

/**
 * Runs either a logistic or quadratic test.
 */
int main() 
{
    // logistic();
    quad();
} // int main()