/**
 * @author Chirag Kaushik
 * @date August 31, 2021
 * 
 */

#include <iostream>
#include <cmath>

#define ENDL "\n"


/**
 * 
 */
int main() 
{
    double a,b,c;
    std::cout << "a: ";
    std::cin >> a;
    std::cout << "b: ";
    std::cin >> b;
    std::cout << "c: ";
    std::cin >> c;
    double r1 = (-b - std::sqrt(b * b - 4 * a * c)) / (2 * a);
    double r2 = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
    std::cout << "root 1: " << r1 << ENDL << "root 2: " << r2 << ENDL;
    double x2 = c / (a * r1);
    std::cout << "root 2 (using root 1): " << x2 << ENDL;

    double plug1 = a * r1 * r1 + b * r1 + c;
    double plug2 = a * r2 * r2 + b * r2 + c;
    double plug3 = a * x2 * x2 + b * x2 + c;

    std::cout << "plugged back in root 1: " << plug1 << ENDL;
    std::cout << "plugged back in root 2: " << plug2 << ENDL;
    std::cout << "plugged back in root 2 (using root 1): " << plug3 << ENDL;
    
} // int main()