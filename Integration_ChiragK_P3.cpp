#include <iostream>
#include <fstream>
#include <cmath>

typedef std::pair<double, double> point;

#define x first
#define y second

double fn1(double x) {
    return 1 + x * x;
}

double fn2(double x) {
    return x * std::exp(-x * x);
}

double fn3(double x) {
    return x * std::exp(-x);
}

double fn4(double x) {
    return std::sin(x);
}

double fn1_int(double x) {
    return x + x*x*x / 3;
}

double fn2_int(double x) {
    return -0.5 * std::exp(-x * x);
}

double fn3_int(double x) {
    return -std::exp(-x) * (x + 1);
}

double fn4_int(double x) {
    return -std::cos(x);
}

double trapezoid(double (*fn)(double), double a, double b, double n) {
    double tot = 0;
    double delta = (b - a) / n;
    for(int i = 0;i < n;i ++) {
        // calculate trapezoidal integral
        double l = a + delta * i;
        double r = l + delta;
        tot += (fn(l) + fn(r)) * delta / 2;
    }
    return tot;
}

point* calcIntegrationCurveTrapezoid(double (*fn)(double), double a, double b, int n) {
    point* points = new point[n+1];
    double delta = (b - a) / n;
    double integral = 0;
    points[0] = std::make_pair(a, 0);
    for(int i = 0;i < n;i ++) {
        double l = a + delta * i;
        double r = l + delta;
        integral += (fn(l) + fn(r)) * delta / 2;
        points[i+1] = std::make_pair(r, integral);
    }
    return points;
}

point* calcIntegrationCurveSimpson(double (*fn)(double), double a, double b, int n) {
    point* points = new point[n+1];
    double delta = (b - a) / n;
    double integral = 0;
    points[0] = std::make_pair(a, 0);
    for(int i = 0;i < n;i ++) {
        double l = a + delta * i;
        double r = l + delta;
        //simpsons rule
        integral += (fn(l) + 4 * fn((l + r) / 2) + fn(r)) * delta / 6;
        points[i+1] = std::make_pair(r, integral);
    }
    return points;
}

void printToFile(std::string filename, point** pts, int n, int numsets) {
    std::ofstream fout(filename);
    for(int i = 0;i < n;i ++) {
        fout << pts[0][i].x;
        for(int j = 0;j < numsets;j ++) {
            fout << "\t" << pts[j][i].y;
        }
        fout << std::endl;
    }
    fout.close();
}

/**
 * Calculates the RMS error between two sets of derivatives.
 * 
 * @param realDerivs the derivative values from the closed form expression
 * @param approxDerivs the approximate derivatives
 * @param len the number of points in the approxDerivs array
 * @return the RMS error between the real and approximate derivs
 */
double calcrms(double (*realInt)(double), point* approxInt, int len) 
{
    double errorsum = 0;
    for (int i = 0; i < len; i++) 
    {
        double realVal = realInt(approxInt[i].x) - realInt(approxInt[0].x);
        double diff = std::abs(realVal * realVal
                               - approxInt[i].y * approxInt[i].y);
        errorsum += diff;
    }
    return std::sqrt(errorsum / len);
}

double percentError(double actual, double approx) {
    return (1 - std::abs(approx / actual)) * 100;
}

int main() {
    int n = 100;
    point* fn1Int_trap = calcIntegrationCurveTrapezoid(fn1, -2, 2, n);
    point* fn2Int_trap = calcIntegrationCurveTrapezoid(fn2, -2, 2, n);
    point* fn3Int_trap = calcIntegrationCurveTrapezoid(fn3, -2, 2, n);
    point* fn4Int_trap = calcIntegrationCurveTrapezoid(fn4, -2, 2, n);

    point* fn1Int_simp = calcIntegrationCurveSimpson(fn1, -2, 2, n);
    point* fn2Int_simp = calcIntegrationCurveSimpson(fn2, -2, 2, n);
    point* fn3Int_simp = calcIntegrationCurveSimpson(fn3, -2, 2, n);
    point* fn4Int_simp = calcIntegrationCurveSimpson(fn4, -2, 2, n);

    // point* fn1Int_real = calcIntegrationCurveTrapezoid(fn1_int, -2, 2, n);
    // point* fn2Int_real = calcIntegrationCurveTrapezoid(fn2_int, -2, 2, n);
    // point* fn3Int_real = calcIntegrationCurveTrapezoid(fn3_int, -2, 2, n);
    // point* fn4Int_real = calcIntegrationCurveTrapezoid(fn4_int, -2, 2, n); 

    point* data_trap[] = {fn1Int_trap, fn2Int_trap, fn3Int_trap, fn4Int_trap};
    point* data_simp[] = {fn1Int_simp, fn2Int_simp, fn3Int_simp, fn4Int_simp};
    printToFile("data/integral_data_trap.txt", data_trap, n, 4);
    printToFile("data/integral_data_simp.txt", data_simp, n, 4);

    std::cout << "rms trap fn1: " << calcrms(fn1_int, fn1Int_trap, n) << std::endl;
    std::cout << "rms trap fn2: " << calcrms(fn2_int, fn2Int_trap, n) << std::endl;
    std::cout << "rms trap fn3: " << calcrms(fn3_int, fn3Int_trap, n) << std::endl;
    std::cout << "rms trap fn4: " << calcrms(fn4_int, fn4Int_trap, n) << std::endl << std::endl;

    std::cout << "\%e trap fn1: " << percentError(fn1_int(fn1Int_trap[n].x) - fn1_int(fn1Int_trap[0].x), fn1Int_trap[n].y) << std::endl;
    std::cout << "\%e trap fn2: " << percentError(fn2_int(fn2Int_trap[n].x) - fn2_int(fn2Int_trap[0].x), fn2Int_trap[n].y) << std::endl;
    std::cout << "\%e trap fn3: " << percentError(fn3_int(fn3Int_trap[n].x) - fn3_int(fn3Int_trap[0].x), fn3Int_trap[n].y) << std::endl;
    std::cout << "\%e trap fn4: " << percentError(fn4_int(fn4Int_trap[n].x) - fn4_int(fn4Int_trap[0].x), fn4Int_trap[n].y) << std::endl << std::endl;

    std::cout << "rms simp fn1: " << calcrms(fn1_int, fn1Int_simp, n) << std::endl;
    std::cout << "rms simp fn2: " << calcrms(fn2_int, fn2Int_simp, n) << std::endl;
    std::cout << "rms simp fn3: " << calcrms(fn3_int, fn3Int_simp, n) << std::endl;
    std::cout << "rms simp fn4: " << calcrms(fn4_int, fn4Int_simp, n) << std::endl << std::endl;
}