/**
 * This file contains two methods of finding the zeroes of a function, the 
 * bisection and the Newton-Ralphson method. It also has various test functions
 * to examine how the methods behavior under different conditions.
 * 
 * @author Chirag Kaushik
 * @version October 12, 2021
 */

#include <iostream>
#include <cmath>

#define ENDL "\n"

/**
 * Calculates x^2
 * 
 * @param x the x value
 * @return x^2
 */
double fn1(double x) 
{
    return x * x;
}

/**
 * Calculates x^2 + 1
 * 
 * @param x the x value
 * @return x^2 + 1
 */
double fn2(double x) 
{
    return 1 + x * x;
}

/**
 * Calculates e^(-x^2)
 * 
 * @param x the x value
 * @return e^(-x^2)
 */
double fn3(double x) 
{
    return std::exp(-x * x);
}

/**
 * Calculates xe^(-x^2)
 * 
 * @param x the x value
 * @return xe^(-x^2)
 */
double fn4(double x) 
{
    return x * std::exp(-x * x);
}

/**
 * Calculates xe^(x^2)
 * 
 * @param x the x value
 * @return xe^(x^2)
 */
double fn5(double x) 
{
    return x * std::exp(x * x);
}

/**
 * Calculates xe^(-x)
 * 
 * @param x the x value
 * @return xe^(-x)
 */
double fn6(double x) 
{
    return x * std::exp(-x);
}

/**
 * Calculates xe^(x)
 * 
 * @param x the x value
 * @return xe^(x)
 */
double fn7(double x) 
{
    return x * std::exp(x);
}

/**
 * Calculates √|x|
 * 
 * @param x the x value
 * @return √|x|
 */
double fn8(double x) 
{
    return std::sqrt(std::abs(x));
}

/**
 * Returns the sign of a number
 * 
 * @param num the number
 * @return the sign of the number. 1 for positive, -1 for negative, and 0 for 0.
 */
double sign(double num) 
{
    if(num > 0) return 1;
    if(num < 0) return -1;
    return 0;
}

/**
 * Runs the bisection algorithm on a function to find a zero. Finds a zero 
 * between the points l and r, inclusive, provided that the points l and r
 * have opposite signs, and there are no asymptotes or discontinuities between
 * the two points. The method runs for a defined number of iterations and 
 * accepts a parameter for the tolerance of the solution.
 * 
 * @param l the x coordinate of the left bound
 * @param r the x cordinate of the right bound
 * @param fn the function to find zeroes of
 * @param iterations the number of iterations to run the algorithm for, each
 *                   iteration halves the range
 * @param tolerance how close the found zero must be to zero. A sufficient
 *                  number of iterations can find an exact zero, due to floating
 *                  point precision.
 * @return the x coordinate of the zero. if an error is encountered, returns nan
 */
double bisection(double l, double r, double (*fn)(double), double iterations, 
                                                           double tolerance) 
{
    double value;
    if(l > r)
    {
        std::cout << "Range is backwards" << ENDL;
        return nan("");
    }
    if(fn(l) == 0) 
    {
        return l;
    }
    if(fn(r) == 0) 
    {
        return r;
    }
    if(sign(fn(l)) == sign(fn(r))) 
    {
        std::cout << "The left and right bound must have opposite sign" << ENDL;
        return nan("");
    }
    for(int i = 0;i < iterations;i ++) 
    {
        double mid = (l + r) / 2;
        if(fn(mid) == INFINITY || fn(l) == INFINITY || fn(r) == INFINITY) 
        {
            std::cout << "Reached infinite asymptote" << ENDL;
            return nan("");
        }
        value = fn(mid);
        // std::cout << l << " " << r << " " << mid << " " << value << ENDL;
        if(std::abs(value) <= tolerance) 
        {
            return mid;
        }
        if(sign(mid) == sign(fn(l))) 
        {
            l = mid;
        } else 
        {
            r = mid;
        }
    }
    std::cout << "Did not find zero, exited at " << ((l+r)/2) << ENDL; 
    return nan("");
}

/**
 * Numerically approximates the derivative of a function using the five point
 * stencil, with a specified step size.
 * 
 * @param fn the function to find the derivative of
 * @param t the coordinate to find the derivative at
 * @param h the step size
 * @return the derivative 
 */
double calcFivePointDeriv(double (*fn)(double), double t, double h) 
{
    return (-fn(t + 2 * h) + 8 * fn(t + h) - 8 * fn(t - h) + fn(t - 2 * h)) / 
           (12 * h);
}

/**
 * Numerically approximates the second derivative of a function using the five 
 * point stencil, with a specified step size. It uses the five point stencil to
 * find the first derivative of 4 adjacent points, then uses those to find the
 * second derivative.
 * 
 * @param fn the function to find the derivative of
 * @param t the coordinate to find the derivative at
 * @param h the step size
 * @return the second derivative 
 */
double calcFivePointDoubleDeriv(double (*fn)(double), double t, double h) 
{
    double p1 = calcFivePointDeriv(fn, t + 2 * h, h);
    double p2 = calcFivePointDeriv(fn, t + h, h);
    double p3 = calcFivePointDeriv(fn, t - h, h);
    double p4 = calcFivePointDeriv(fn, t - 2 * h, h);
    return (-p1 + 8 * p2 - 8 * p3 + p4) / (12 * h);
}

/**
 * Uses the Newton-Raphson method to find the zero of a function in a specific
 * range, based on the derivative of the function calculated using the five
 * point stencil. It takes in a starting point as well as a range to find the
 * zero over. If the method exits the range, it will throw an error. It will 
 * also throw an error if it finds a stationary point where the derivative is 
 * zero, or if the algorithm does not find a root. 
 * 
 * @param t the start point
 * @param l the left bound
 * @param r the right bound
 * @param fn the function to find a zero of
 * @param stepSize the step size of the five point stencil
 * @param iterations the max number of iterations to run the algorithm for
 * @param tolerance how close the found zero must be to zero. A sufficient
 *                  number of iterations can find an exact zero, due to floating
 *                  point precision.
 * @return the x coordinate of the zero. if an error is encountered, returns nan
 */
double newton(double t, double l, double r, double (*fn)(double), 
              double stepSize, double iterations, double tolerance) 
{
    for(int i = 0;i < iterations;i ++) 
    {
        if(std::abs(fn(t)) <= tolerance) 
        {
            return t;
        }  
        double deriv = calcFivePointDeriv(fn, t, stepSize);
        // std::cout << t << " " << deriv << ENDL;
        if(std::abs(deriv) < 1e-11) 
        {
            std::cout << "stationary point or asymptote at " << t << ENDL;
            return nan("");
        }
        t = t - fn(t) / deriv;
        if(t < l || t > r) 
        {
            std::cout << "Exited bounds at " << t << ENDL;
            return nan("");
        }
    }
    std::cout << "Newton algorithm did not find root, exited at " << t << ENDL;
    return nan("");
}


double (*currentFn)(double);
double currentStepSize;
/**
 * Calculates the derivative of a function at a specific time. Takes in only the
 * time as a parameter, so it can be wrapped as a lambda into the methods for
 * finding zeroes. The function and step size are set outside of the function.
 * 
 * @param t the time to take the derivative
 * @return the derivative
 */
double minDeriv(double t) {
    return calcFivePointDeriv(currentFn, t, currentStepSize);
}

/**
 * Finds a relative minima of a function. Takes in an initial starting point
 * estimate, and a range over which the algorithm will run. First, it finds
 * a zero of the first derivative of the input function by approximating with
 * the five point stencil. Tries to use both the bisection and newton algorithms
 * to find the zero, and errors if no zero is find. Then, it checks the sign of
 * the double derivative of the function to see if a local maxima, minima, or
 * saddle point was found. If a minima is not found, it errors.
 * 
 * @param t the initial time estimate
 * @param l the left bound to search over
 * @param r the right bound to search over
 * @param fn to function to find minimum of
 * @param stepSize the step size of the five point stencil
 * @param iterations number of iterations to run the zero algorithms for
 * @param tolerance how close to zero the zero must be
 * @return the time value of the minima, or nan if there is an error.
 */
double findMinimum(double t, double l, double r, double (*fn)(double),
                   double stepSize, double iterations, double tolerance) 
{
    currentFn = fn;
    currentStepSize = stepSize;
    double zero = newton(t, l, r, minDeriv, stepSize, iterations, tolerance);
    if(std::isnan(zero))
    {
        std::cout << "Could not find zero with newton, using bisection" << ENDL;
        zero = bisection(l, r, minDeriv, iterations, tolerance);
        if(std::isnan(zero))
        {
            std::cout << "Could not find zero with either algorithm" << ENDL;
            return nan("");
        }
    }
    double secondDeriv = calcFivePointDoubleDeriv(fn, t, stepSize);
    // std::cout << secondDeriv << ENDL;
    if(secondDeriv == 0)
    {
        std::cout << "Found saddle point at " << zero << ENDL;
        return nan("");
    }
    else if(secondDeriv < 0)
    {
        std::cout << "Found local max at " << zero << ENDL; 
        return nan("");
    }
    else
    {
        return zero;
    }
}

/**
 * Runs the bisection and newton algorithms on a variety of functions to 
 * examine their behavior.
 */
int main() 
{
    // std::cout << "BISECTION" << ENDL << "=========" << ENDL << ENDL;
    // std::cout << "x^2" << ENDL;
    // std::cout << bisection(-3, 2.1, fn1, 100000, 0) << ENDL << ENDL;
    // std::cout << "x^2 + 1" << ENDL;
    // std::cout << bisection(-3, 2.1, fn2, 100000, 0) << ENDL << ENDL;
    // std::cout << "e^(-x^2)" << ENDL;
    // std::cout << bisection(-3, 2.1, fn3, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(-x^2)" << ENDL;
    // std::cout << bisection(-3, 2.1, fn4, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(x^2)" << ENDL;
    // std::cout << bisection(-3, 2.1, fn5, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(-x)" << ENDL;
    // std::cout << bisection(-3, 2.1, fn6, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(x)" << ENDL;
    // std::cout << bisection(-3, 2.1, fn7, 100000, 0) << ENDL << ENDL;
    // std::cout << "√|x|" << ENDL;
    // std::cout << bisection(-3, 2.1, fn8, 100000, 0) << ENDL << ENDL;

    // std::cout << "NEWTON" << ENDL << "=======" << ENDL << ENDL;
    // std::cout << "x^2" << ENDL;
    // std::cout << newton(0.2, -10, 10, fn1, 0.01, 10000, 1e-10) << ENDL << ENDL;
    // std::cout << "x^2 + 1" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn2, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "e^(-x^2)" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn3, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(-x^2)" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn4, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(x^2)" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn5, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(-x)" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn6, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "xe^(x)" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn7, 0.001, 100000, 0) << ENDL << ENDL;
    // std::cout << "√|x|" << ENDL;
    // std::cout << newton(0.3, -10, 10, fn8, 0.001, 100000, 0) << ENDL << ENDL;

    std::cout << "MINIMUM" << ENDL << "=======" << ENDL << ENDL;
    std::cout << "x^2" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn1, 0.001, 100000, 0) << ENDL << ENDL;
    std::cout << "x^2 + 1" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn2, 0.001, 100000, 1e-10) << ENDL << ENDL;
    std::cout << "e^(-x^2)" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn3, 0.001, 100000, 1e-10) << ENDL << ENDL;
    std::cout << "xe^(-x^2)" << ENDL;
    std::cout << findMinimum(-1, -10, 10, fn4, 0.001, 100000, 1e-10) << ENDL << ENDL;
    std::cout << "xe^(x^2)" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn5, 0.001, 100000, 0) << ENDL << ENDL;
    std::cout << "xe^(-x)" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn6, 0.001, 100000, 1e-10) << ENDL << ENDL;
    std::cout << "xe^(x)" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn7, 0.001, 100000, 1e-10) << ENDL << ENDL;
    std::cout << "√|x|" << ENDL;
    std::cout << findMinimum(0.3, -10, 10, fn8, 0.001, 100000, 1e-10) << ENDL << ENDL;
}