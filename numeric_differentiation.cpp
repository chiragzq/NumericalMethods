/**
 * This file contains tests for various methods of approximating the 
 * derivative of a function over a set a points, and calculating the
 * RMS error of the approximation compared to the closed form 
 * expression for the derivative.
 * 
 * @author Chirag Kaushik
 * @version September 17, 2021
 */

#include <iostream>
#include <fstream>
#include <cmath>

#define ENDL "\n"
#define NUM_POINTS 1000

/**
 * Point struct that stores a time and y value.
 */
struct point
{
    double t;
    double y;
};

/**
 * Calculates the value of e^(-(t^2)).
 * 
 * @param t the time value
 * @return the value of the function
 */
double calculateF(double t) {
    return std::exp(-t * t);
}

/**
 * Calculates the closed form derivative of the previous function 
 * = -2t*e^(-(t^2)).
 * 
 * @param t the time value
 * @return the derivative of the function
 */
double calculateF2(double t) {
    return -2 * t * std::exp(-t * t);
}

/**
 * Prints a list of points to a file by printing the time and y 
 * value separated by a space.
 * 
 * @param name the name of the file
 * @param values the array of points
 * @param numPoints the number of points to print out
 */
void printToFile(std::string name, point* values, int numPoints) {
    std::ofstream fout(name);
    for (int i = 0;i < numPoints;i ++) 
    {
        fout << values[i] << ENDL;
    }
}

/**
 * Calculates the RMS error between two sets of derivatives.
 * 
 * @param realDerivs the derivative values from the closed form expression
 * @param approxDerivs the approximate derivatives
 * @param offset since some calculations for the derivatives use adjacent 
 *               points, the endpoints are cut off, so the offset represents
 *               where to start the comprisons in the second array
 * @param len the number of points in the approxDerivs array
 * @return the RMS error between the real and approximate derivs
 */
double calcrms(point* realDerivs, point* approxDerivs, int offset, int len) {
    double errorsum = 0;
    for(int i = 0;i < len;i ++) {
        double diff = std::abs(realDerivs[offset + i].y * 
                               realDerivs[offset + i].y
                               - approxDerivs[i].y * approxDerivs[i].y);
        errorsum += diff;
    }
    return std::sqrt(errorsum / len);
}

/**
 * The default c++ print function overloads the bit left shift operator, but
 * there is no equivalent of toString in java that a custom class can override
 * to get a nice formatted string output. In c++ the actual printing function
 * must be override with implementation for the custom class. This function 
 * achieves that.
 * 
 * @param os the output stream being printed to
 * @param p the point to print
 * @return the output stream, for method chaining
 */
std::ostream& operator<<(std::ostream& os, const point& p) {
    os << p.t << ' ' << p.y;
    return os;
}

/**
 * Calculates numPoints evenly spaced out points of a function over the range
 * [lower, upper).
 * 
 * @param lower the inclusive lower bound of the range
 * @param upper the exclusive upper bound of the range
 * @param numPoints the number of points to generate
 * @param fn the function to calculate values with
 * @return the resulting array of points
 */
point* calculateValues(double lower, 
                       double upper, int numPoints, double (*fn)(double)) 
{
    point* ret = new point[numPoints];
    double step = (upper - lower) / (numPoints );
    for (int i = 0;i < numPoints;i ++) 
    {
        double xCoord = lower + step * i;
        double y = (*fn)(xCoord);
        ret[i] = {xCoord, y};
    }
    return ret;
}

/**
 * Calculates the approximate derivative by taking the slope between adjacent 
 * points. Returns a tuple containing three arrays which have the same 
 * derivative values, by assign the time value to the left time value, right
 * time value, and midpoint of the time values, respectively.
 * 
 * @param points the list of points to calculate derivatives of
 * @param numPoints the number of points in the array
 * @return a tuple containing the three derivative arrays
 */
std::tuple<point*, point*, point*> simpleSlopeDerivs(point* points, 
                                                     int numPoints) {
    point* leftDerivs = new point[numPoints - 1];
    point* rightDerivs = new point[numPoints - 1];
    point* midDerivs = new point[numPoints - 1];
    for(int i = 0;i < numPoints - 1;i ++) {
        double deriv = (points[i + 1].y - points[i].y) / 
                       (points[i + 1].t - points[i].t);
        leftDerivs[i] = {points[i].t, deriv};
        rightDerivs[i] = {points[i + 1].t, deriv};
        midDerivs[i] = {(points[i].t + points[i + 1].t) / 2, deriv};
    }
    std::tuple<point*, point*, point*> ret(leftDerivs, rightDerivs, midDerivs);
    return ret;
}

/**
 * Approximates the derivative for a point by using the slope between the
 * previous and next point in the array.
 * 
 * @param points the array of points
 * @param numPoints the number of points in the array
 * @return an array of derivatives
 */
point* threePointDeriv(point* points, int numPoints) {
    point* ret = new point[numPoints - 2];
    for(int i = 1;i < numPoints - 1;i ++) {
        double deriv = (points[i + 1].y - points[i - 1].y) / 
                       (points[i + 1].t - points[i - 1].t);
        ret[i - 1] = {(points[i - 1].t + points[i + 1].t) / 2, deriv};
    }
    return ret;
}

/**
 * Calculates the parabolic fit between three points, and returns the 
 * coefficients of the quadratic and linear term. Used mathematica to calculate
 * a closed form expression for the coefficients.
 * 
 * @param p1 the first point
 * @param p2 the second point
 * @param p3 the third point
 * @return a pair with two values, which are the quadratic and linear term.
 */
std::pair<double, double> calcParabolaCoeffs(point p1, point p2, point p3) {
    double t1 = p1.t;
    double y1 = p1.y;
    double t2 = p2.t;
    double y2 = p2.y;
    double t3 = p3.t;
    double y3 = p3.y;
    double a = (t3 * (-y1 + y2) + t2 * (y1 - y3) + t1 * (-y2 + y3)) /
               ((t1 - t2) * (t1 - t3) * (t2 - t3));
    double b = (t3 * t3 * (y1 - y2) + t1 * t1 * (y2 - y3) + 
                                                t2 * t2 * (-y1 + y3)) /
               ((t1 - t2) * (t1 - t3) * (t2 - t3));
    return std::make_pair(a, b);
}

/**
 * Approximates the derivative for a point by fitting the point and its two 
 * neighbors to a parabolic function, then taking the closed form derivative
 * of that parabola. For endpoints, uses the 3 points on the edge of the array.
 * @param points the array of points 
 * @param numPoints the number of points in the array
 * @return an array of derivatives
 */
point* calcParabolaDerivs(point* points, int numPoints) {
    point* ret = new point[numPoints];
    for(int i = 0;i < numPoints;i ++) {
        int centerIndex = i;
        if(centerIndex == 0)
            centerIndex = 1;
        else if(centerIndex == numPoints - 1)
            centerIndex = numPoints - 2;
        
        std::pair<double, double> parabolaCoeffs = 
            calcParabolaCoeffs(points[centerIndex - 1], 
                               points[centerIndex], 
                               points[centerIndex + 1]);
        double a = parabolaCoeffs.first;
        double b = parabolaCoeffs.second;
        ret[i] = {points[i].t, 2 * a * points[i].t + b};
    }
    return ret;
}

/**
 * Approximates derivatives using the five point stencil, which takes 
 * the two neighbors on each side of a point for the approximation.
 * 
 * @param points an array of points
 * @param numPoints the number of points in the array
 * @return an array of derivatives
 */
point* calcFivePointDeriv(point* points, int numPoints) {
    point* ret = new point[numPoints - 4];
    for(int i = 2;i < numPoints - 2;i ++) {
        double derivative = (-points[i + 2].y + 
                             8 * points[i + 1].y - 
                             8 * points[i - 1].y + 
                             points[i - 2].y) / 
                             (12 * (points[i].t - points[i+1].t));
        ret[i - 2] = {points[i].t, derivative};
    }
    return ret;
}


/**
 * Calculates simple slope derivatives and prints out the RMS error.
 */
void runSimpleSlopeDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    std::tuple<point*, point*, point*> derivs = simpleSlopeDerivs(values, 
                                                                  NUM_POINTS);
    point* leftDerivs = std::get<0>(derivs);
    point* rightDerivs = std::get<1>(derivs);
    point* midDerivs = std::get<2>(derivs);

    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    point* realMidDerivs = calculateValues(-9.9, 10.1, NUM_POINTS, calculateF2);
    printToFile("pointdata.out", values, NUM_POINTS);
    printToFile("leftderivs.out", leftDerivs, NUM_POINTS - 1);
    printToFile("rightderivs.out", rightDerivs, NUM_POINTS - 1);
    printToFile("midderivs.out", midDerivs, NUM_POINTS - 2);
    printToFile("realmidderivs.out", realMidDerivs, NUM_POINTS - 2);
    printToFile("realderivs.out", realDerivs, NUM_POINTS);

    std::cout << calcrms(realDerivs, leftDerivs, 0, NUM_POINTS - 1) << ENDL;
    std::cout << calcrms(realDerivs, rightDerivs, 1, NUM_POINTS - 1) << ENDL;
    std::cout << calcrms(realMidDerivs, midDerivs, 0, NUM_POINTS - 2) << ENDL;
}

/**
 * Calculates the three point derivative and prints the RMS error
 */ 
void runThreePointDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    point* threePtDerivs = threePointDeriv(values, NUM_POINTS);
    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    printToFile("pointdata.out", values, NUM_POINTS);
    printToFile("realderivs.out", realDerivs, NUM_POINTS - 2);
    printToFile("threepointderivs.out", threePtDerivs, NUM_POINTS - 2);
    std::cout << calcrms(realDerivs, threePtDerivs, 1, NUM_POINTS - 2) << ENDL;
}

/**
 * Calculates parabolic approximation derivative and prints the RMS error
 */
void runParabolaDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    point* parabolaDerivs = calcParabolaDerivs(values, NUM_POINTS);
    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    printToFile("pointdata.out", values, NUM_POINTS);
    printToFile("realderivs.out", realDerivs, NUM_POINTS);
    printToFile("paraboladerivs.out", parabolaDerivs, NUM_POINTS);
    std::cout << calcrms(realDerivs, parabolaDerivs, 0, NUM_POINTS) << ENDL;
}

/**
 * Calculates the 5 point stencil derivative and prints the RMS error
 */
void runFivePointStencil() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    point* fivePtDerivs = calcFivePointDeriv(values, NUM_POINTS);
    printToFile("fivepointderivs.out", fivePtDerivs, NUM_POINTS - 4);
    std::cout << calcrms(realDerivs, fivePtDerivs, 2, NUM_POINTS - 4) << ENDL;
}

/**
 * In the main method one of the derivative procedures can be run.
 */ 
int main() {
    // runSimpleSlopeDeriv();
    // runThreePointDeriv();
    // runParabolaDeriv();
    runFivePointStencil();
}