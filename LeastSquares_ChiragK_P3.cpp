#include <iostream>
#include <fstream>
#include <cmath>

// read in an array of points based on filename
void readPointsFromFile(std::string filename, int n, double* x, double* y) {
	std::ifstream file;
	file.open(filename);
	if(!file.is_open()){
		std::cout << "Error: file not found" << std::endl;
		return;
	}
	for(int i = 0; i < n; i++)
    {
		file >> x[i] >> y[i];
	}
	file.close();
}

double getTotalError(int n, int numQ, double* x, double* y, double* q, double (*model)(double, double*)) {
	double totalError = 0;
	for(int i = 0; i < n; i++) {
		// std::cout << "x: " << x[i] << " y: " << y[i] << "mode: " << model(x[i], q) << std::endl;
		double error = model(x[i], q) - y[i];
		totalError += error * error;
	}
	return 0.5 * totalError;
}

double modelFunctionLinear(double x, double* q) {
    return q[0] + q[1] * x;
}

double modelFunctionGaussian(double x, double* q) {
	return q[0] * q[0] * std::exp(-(x - q[1]) * (x - q[1]) / (q[2] * q[2])) + q[3] * q[3];
}

double modelFunctionSin(double x, double* q) {
	return q[0] * std::sin(q[1] * x + q[2]) + q[3];
}

double fivePointPartial(double (*model)(double, double*), double x, double* q, int j, double h) {
	double qj = q[j];
	q[j] = qj + 2 * h;
	double a = model(x, q);
	q[j] = qj + h;
	double b = model(x, q);
	q[j] = qj - h;
	double c = model(x, q);
	q[j] = qj - 2 * h;
	double d = model(x, q);
	q[j] = qj;
	return (-a + 8 * b - 8 * c + d) / (12 * h);
}

/**
 * 
 * @param n the number of data points
 * @param numQ the number of parameters in the model function
 * @param x array of x values
 * @param y array of y values
 * @param q array of q parameters
 * @param lambda the learning rate
 * @param iterations number of iterations to run least squares for
 * @param model the model function with q parameters
 * @param partial the partial derivative function with q parameters
 */
void train(int n, int numQ, double* x, double* y, double* q, double lambda, double iterations, 
             double (*model)(double, double*), double tolerance) {
	// std::cout << n << std::endl;
	double MOMENTUM = 0.9;
	std::cout << "it " << "starting" << " error: " << getTotalError(n, numQ, x, y, q, model) << std::endl;
	double* lastDeltas = new double[n];
	double prevError = 0;
	for(int i = 0; i < n; i++) {
		lastDeltas[i] = 0;
	}
	for(int k = 0;k < iterations;k ++) {
		for(int j = 0;j < numQ;j ++) {
			double sum = 0;
			for(int i = 0;i < n;i ++) {
				// sum += (y[i] - model(x[i], q)) * partial(x[i], q, j);
				sum += (y[i] - model(x[i], q)) * fivePointPartial(model, x[i], q, j, 0.001);
			}
			double deltaQ = sum * lambda + MOMENTUM * lastDeltas[j];
			q[j] += deltaQ;
			lastDeltas[j] = deltaQ;
		}
		double currentError = getTotalError(n, numQ, x, y, q, model);
		if(k % 100000 == 0) {
			std::cout << "it " << k << " error: " << currentError << std::endl;
			std::cout << "q0: " << q[0] << " q1: " << q[1] << " q2: " << q[2] << " q3: " << q[3] << std::endl;
			// std::cout << abs(currentError - prevError) << std::endl;
		}
		if(abs(currentError - prevError) < tolerance) {
			std::cout << "total it, " << k << " error: " << currentError << std::endl;
			std::cout << "final q0: " << q[0] << " q1: " << q[1] << " q2: " << q[2] << " q3: " << q[3] << std::endl;
			return;
		}
		prevError = currentError;
	}
}

void gaussiantest() {
	int numPoints = 42;
	int numQ = 4;
	double lr = 1e-10;
	int iterations = 100000000;

    double* x = new double[numPoints];
	double* y = new double[numPoints];
	readPointsFromFile("./GeigerHistoOriginal.txt", numPoints, x, y);
	double* q = new double[numQ];

	for(int i = 0;i < numQ;i ++) {
		q[i] = 1;
	}
	
	train(numPoints, numQ, x, y, q, lr, iterations, modelFunctionGaussian, 0.1);
	std::cout << "q0: " << q[0] << " q1: " << q[1] << " q2: " << q[2] << " q3: " << q[3] << std::endl;
}


void sinHe1Test() {
	int numPoints = 15;
	int numQ = 4;
	double lr = 1e-5;
	int iterations = 100000000;

    double* x = new double[numPoints];
	double* y = new double[numPoints];
	readPointsFromFile("./UCrBHe1.txt", numPoints, x, y);
	double* q = new double[numQ];
	
	for(int i = 0;i < numQ;i ++) {
		q[i] = 1;
	}
	
	train(numPoints, numQ, x, y, q, lr, iterations, modelFunctionSin, 1e-10);
	// std::cout << "q0: " << q[0] << " q1: " << q[1] << " q2: " << q[2] << " q3: " << q[3] << std::endl;
}

void sinMg2Test() {
	int numPoints = 15;
	int numQ = 4;
	double lr = 1e-6;
	int iterations = 100000000;

    double* x = new double[numPoints];
	double* y = new double[numPoints];
	readPointsFromFile("./UCrBMg2.txt", numPoints, x, y);
	double* q = new double[numQ];
	
	for(int i = 0;i < numQ;i ++) {
		q[i] = 1;
	}
	
	train(numPoints, numQ, x, y, q, lr, iterations, modelFunctionSin, 1e-10);
	std::cout << "q0: " << q[0] << " q1: " << q[1] << " q2: " << q[2] << " q3: " << q[3] << std::endl;
}

int main() {
	sinMg2Test();
}



