#include <iostream>
#include <fstream>
#include <cmath>

#define ENDL "\n"
#define x first
#define y second

#define NUM_POINTS 100

typedef std::pair<double, double> point;


double calculateF(double t) {
    return std::exp(-t * t);
}

double calculateF2(double t) {
    return -2 * t * std::exp(-t * t);
}

std::ostream& operator<<(std::ostream& _os, const point& _p) {
    _os << _p.first << ' ' << _p.second;
    return _os;
}

point* calculateValues(double lower, double upper, int numPoints, double (*fn)(double)) 
{
    point* ret = new point[numPoints];
    double step = (upper - lower) / numPoints;
    for (int i = 0;i < numPoints;i ++) 
    {
        double xCoord = lower + step * i;
        double y = (*fn)(xCoord);
        ret[i] = point(xCoord, y);
    }
    return ret;
}

std::tuple<point*, point*, point*> calculateDerivs(point* points, int numPoints) {
    point* leftDerivs = new point[numPoints - 1];
    point* rightDerivs = new point[numPoints - 1];
    point* midDerivs = new point[numPoints - 1];
    for(int i = 0;i < numPoints - 1;i ++) {
        double deriv = (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
        leftDerivs[i] = point(points[i].x, deriv);
        rightDerivs[i] = point(points[i + 1].x, deriv);
        midDerivs[i] = point((points[i].x + points[i + 1].x) / 2, deriv);
    }
    std::tuple<point*, point*, point*> ret(leftDerivs, rightDerivs, midDerivs);
    return ret;
}

// point* calculateDerivativesRight(point* points, int numPoints) {
//     point* ret = new point[numPoints - 1];
//     for(int i = 0;i < numPoints - 1;i ++) {
//         double deriv = (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
//         ret[i] = point(points[i + 1].x, deriv);
//     }
//     return ret; 
// }

// point* calculateDerivativesMidpoint(point* points, int numPoints) {
//     point* ret = new point[numPoints - 2];
//     for(int i = 1;i < numPoints - 1;i ++) {
//         double deriv = (points[i + 1].y - points[i - 1].y) / (points[i + 1].x - points[i - 1].x);
//         ret[i - 1] = point(points[i].x, deriv);
//     }
//     return ret;
// }

void printToFile(std::string name, point* values, point* leftDerivs, point* rightDerivs, point* midDerivs, point* realDerivs, int numPoints) {
    std::ofstream fout(name);
    for (int i = 0;i < numPoints;i ++) 
    {
        fout << values[i].x << " " << values[i].y;
        if(i != numPoints - 1)
            fout << " " <<  leftDerivs[i].y;
        else
            fout << " -";

        if(i != 0)
            fout << " " << rightDerivs[i - 1].y;
        else
            fout << " -";
        if(i != 0 && i != numPoints - 1)
            fout << " " << midDerivs[i - 1].y;
        else 
            fout << " -";

        fout << " " << realDerivs[i].y << ENDL;
    }
}

double calcrms(point* realDerivs, point* approxDerivs, int offset, int len) {
    double errorsum = 0;
    for(int i = 0;i < len;i ++) {
        // std::cout << realDerivs[offset + i].y << " " << approxDerivs[i].y << ENDL;
        double diff = std::abs(realDerivs[offset + i].y * realDerivs[offset + i].y - approxDerivs[i].y * approxDerivs[i].y);
        errorsum += diff;
    }
    return std::sqrt(errorsum / len);

}


int main() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    std::tuple<point*, point*, point*> derivs = calculateDerivs(values, NUM_POINTS);
    point* leftDerivs = std::get<0>(derivs);
    point* rightDerivs = std::get<1>(derivs);
    point* midDerivs = std::get<2>(derivs);

    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    printToFile("pointdata1.out", values, leftDerivs, rightDerivs, midDerivs, realDerivs, NUM_POINTS);

    std::cout << calcrms(realDerivs, leftDerivs, 0, NUM_POINTS - 1) << ENDL;
    std::cout << calcrms(realDerivs, rightDerivs, 1, NUM_POINTS - 1) << ENDL;
    std::cout << calcrms(realDerivs, midDerivs, 1, NUM_POINTS - 2) << ENDL;


}