#include <iostream>
#include <fstream>
#include <cmath>

#define ENDL "\n"
#define t first
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
    double step = (upper - lower) / (numPoints );
    for (int i = 0;i < numPoints;i ++) 
    {
        double xCoord = lower + step * i;
        double y = (*fn)(xCoord);
        ret[i] = point(xCoord, y);
    }
    return ret;
}

std::tuple<point*, point*, point*> simpleSlopeDerivs(point* points, int numPoints) {
    point* leftDerivs = new point[numPoints - 1];
    point* rightDerivs = new point[numPoints - 1];
    point* midDerivs = new point[numPoints - 1];
    for(int i = 0;i < numPoints - 1;i ++) {
        double deriv = (points[i + 1].y - points[i].y) / (points[i + 1].t - points[i].t);
        leftDerivs[i] = point(points[i].t, deriv);
        rightDerivs[i] = point(points[i + 1].t, deriv);
        midDerivs[i] = point((points[i].t + points[i + 1].t) / 2, deriv);
    }
    std::tuple<point*, point*, point*> ret(leftDerivs, rightDerivs, midDerivs);
    return ret;
}

point* threePointDeriv(point* points, int numPoints) {
    point* ret = new point[numPoints - 2];
    for(int i = 1;i < numPoints - 1;i ++) {
        double deriv = (points[i + 1].y - points[i - 1].y) / (points[i + 1].t - points[i - 1].t);
        ret[i - 1] = point((points[i - 1].t + points[i + 1].t) / 2, deriv);
    }
    return ret;
}

void printToFile(std::string name, point* values, int numPoints) {
    std::ofstream fout(name);
    for (int i = 0;i < numPoints;i ++) 
    {
        fout << values[i] << ENDL;
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

std::pair<double, double> calcParabolaCoeffs(point p1, point p2, point p3) {
    double t1 = p1.t;
    double y1 = p1.y;
    double t2 = p2.t;
    double y2 = p2.y;
    double t3 = p3.t;
    double y3 = p3.y;
    double a = (t3 * (-y1 + y2) + t2 * (y1 - y3) + t1 * (-y2 + y3)) /
               ((t1 - t2) * (t1 - t3) * (t2 - t3));
    double b = (t3 * t3 * (y1 - y2) + t1 * t1 * (y2 - y3) + t2 * t2 * (-y1 + y3)) /
               ((t1 - t2) * (t1 - t3) * (t2 - t3));
    return std::make_pair(a, b);
}

point* calcParabolaDerivs(point* points, int numPoints) {
    point* ret = new point[numPoints];
    for(int i = 0;i < numPoints;i ++) {
        int centerIndex = i;
        if(centerIndex == 0)
            centerIndex = 1;
        else if(centerIndex == numPoints - 1)
            centerIndex = numPoints - 2;
        
        std::pair<double, double> parabolaCoeffs = calcParabolaCoeffs(points[centerIndex - 1], points[centerIndex], points[centerIndex + 1]);
        double a = parabolaCoeffs.first;
        double b = parabolaCoeffs.second;

        ret[i] = point(points[i].t, 2 * a * points[i].t + b);
    }
    return ret;
}


void runSimpleSlopeDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    std::tuple<point*, point*, point*> derivs = simpleSlopeDerivs(values, NUM_POINTS);
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

void runThreePointDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    point* threePointDerivs = threePointDeriv(values, NUM_POINTS);
    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    printToFile("pointdata.out", values, NUM_POINTS);
    printToFile("realderivs.out", realDerivs, NUM_POINTS - 2);
    printToFile("threepointderivs.out", threePointDerivs, NUM_POINTS - 2);
    std::cout << calcrms(realDerivs, threePointDerivs, 1, NUM_POINTS - 2) << ENDL;
}

void runParabolaDeriv() {
    point* values = calculateValues(-10, 10, NUM_POINTS, calculateF);
    point* parabolaDerivs = calcParabolaDerivs(values, NUM_POINTS);
    point* realDerivs = calculateValues(-10, 10, NUM_POINTS, calculateF2);
    printToFile("pointdata.out", values, NUM_POINTS);
    printToFile("realderivs.out", realDerivs, NUM_POINTS);
    printToFile("paraboladerivs.out", parabolaDerivs, NUM_POINTS);
    std::cout << calcrms(realDerivs, parabolaDerivs, 0, NUM_POINTS) << ENDL;
}


int main() {
    // runParabolaDeriv();
    // runSimpleSlopeDeriv();
    runParabolaDeriv();
}