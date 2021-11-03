#include <iostream>
#include <fstream>

// read in an array of points based on filename
void readPointsFromFile(std::string filename, int &n, double* x, double* y) {
	std::ifstream file;
	file.open(filename);
	if(!file.is_open()){
		std::cout << "Error: file not found" << std::endl;
		return;
	}
	file >> n;
	for(int i = 0; i < n; i++)
    {
		file >> x[i] >> y[i];
	}
	file.close();
}

double train(int n, int nq, double* x, double* y, double* q, double lambda, double iterations, 
             double (*fn)(double, double*), double (*partial)(double, double*, int)) {
	for(int i = 0;i < iterations;i ++) {
		for(int j = 0;j < nq;j ++) {
			double sum = 0;
			for(int k = 0;k < n;k ++) {
				sum += fn(q[j], x + k * 2) - y[k];
			}
		}
	}
}

double modelFunction(double x, double* q) {
    return q[0] + q[1] * x;
}

double modelFunctionPartial(double x, double* q, int i) {
    if(i == 0)
        return 1;
    else
        return x;
}


int main() {
    
}



