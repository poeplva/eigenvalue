#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace std;

class Matrix {
	int m, n;                                                  // m is the number of rows, and n is the number of columns
	double** arr;
public:
	// constructors (no default constructor)
	Matrix (int m_1, int n_1) : m(m_1), n(n_1) {               // constructor that creates m by n matrix
		arr = new double*[m];
		for (int i = 0; i < m; i++) arr[i] = new double[n];
	}
	Matrix (int n_1) : m(n_1), n(n_1) {                        // constructor that creates n by n (square) matrix
		arr = new double*[m];
		for (int i = 0; i < m; i++) arr[i] = new double[n];
	}

	// operator overloads
	// overload addition
	friend Matrix operator+ (const Matrix&, const Matrix&);
	// overload subtraction
	friend Matrix operator- (const Matrix&, const Matrix&);
	// overload matrix multiplication
	friend Matrix operator* (const Matrix&, const Matrix&);
	// overload scalar multiplication
	friend Matrix operator* (const double&, const Matrix&);

	// return jth column as an array
	double* column(int j) const;

	// take the transpose of given matrix, will be useful when creating the householder transformation
	Matrix transpose() const;

	// bottom right (m-1) by (n-1) part of the matrix
	Matrix bottom_right() const;

	// set values for the matrix
	void set_array (double** arr_1) {
		arr = arr_1;
	}

	// return matrix row size
	int get_m() const {
		return m;
	}
	// return matrix column size
	int get_n() const {
		return n;
	}
	// return the absolute maximal entry of the matrix
	double get_max() const {
		double max = arr[0][0];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (abs(arr[i][j]) > abs(max)) {
					max = arr[i][j];
				}
			}
		}
		return max;
	}
	// return the matrix array
	double** get_arr() const {
		return arr;
	}

	// print the matrix to given file handler
	void print(ostream& o) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				o << arr[i][j];
				if (j != n - 1) o << " ";
			}
			o << endl;
		}
	}

};

double* Matrix::column(int j) const {
	double* col = new double[m];

	for (int i = 0; i < m; i++) {
		col[i] = arr[i][j];
	}
	return col;
}

Matrix Matrix::transpose() const {
	Matrix trans(n, m);

	double** trarr = new double*[n];
	for (int i = 0; i < n; i++) {
		trarr[i] = this->column(i);
	}
	trans.set_array(trarr);
	return trans;
}

Matrix Matrix::bottom_right() const {
	Matrix C(m - 1, n - 1);
	double** C_arr = new double*[m - 1];
	for (int i = 0; i < m - 1; i++) {
		C_arr[i] = new double[n - 1];
		for (int j = 0; j < n - 1; j++) {
			C_arr[i][j] = arr[i + 1][j + 1];
		}
	}
	C.set_array(C_arr);
	return C;
}

Matrix operator+ (const Matrix& A, const Matrix& B) {
	if (A.n != B.n || A.m != B.m) {
		throw invalid_argument("matrices to add must have the same size");
	}
	int m = A.m;
	int n = A.n;

	Matrix C(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C.arr[i][j] = A.arr[i][j] + B.arr[i][j];
		}
	}
	return C;
}

Matrix operator* (const Matrix& A, const Matrix& B) {
	if (A.n != B.m) {
		throw invalid_argument("matrices to mutliply should have compatible sizes");
	}

	Matrix C(A.m, B.n);
	for (int i = 0; i < A.m; i++) {
		for (int j = 0; j < B.n; j++) {
			C.arr[i][j] = 0;
			for (int k = 0; k < A.n; k++) {
				C.arr[i][j] += A.arr[i][k] * B.arr[k][j];
			}
		}
	}
	return C;
}

Matrix operator* (const double& b, const Matrix& A) {
	Matrix C(A.m, A.n);
	for (int i = 0; i < A.m; i++) {
		for (int j = 0; j < A.n; j++) {
			C.arr[i][j] = b * A.arr[i][j];
		}
	}
	return C;
}

Matrix operator- (const Matrix& A, const Matrix& B) {
	return A + (-1) * B;
}

// return the n by n identity matrix
Matrix id (int n) {
	Matrix A(n);
	double** arr = new double*[n];
	for (int i = 0; i < n; i++) {
		arr[i] = new double[n];
		for (int j = 0; j < n; j++) {
			arr[i][j] = (i == j) ? 1 : 0;
		}
	}
	A.set_array(arr);
	return A;
}

// implement normalized power iteration algorithm, return the eigenvalue and modify eigenvec
double normalized_power_iteration (const Matrix A, Matrix& eigenvec, double tolerance) {
	if (A.get_m() != A.get_n() || A.get_m() != eigenvec.get_m() || eigenvec.get_n() != 1) {
		throw invalid_argument("there is an error in the dimensions of the matrix and/or the vector");
	}

	int n = A.get_n();

	double** eigenarr = new double*[n];
	for (int i = 0; i < n; i++) {                              // make your first guess e_1
		eigenarr[i] = new double[1];
		eigenarr[i][0] = (i == 0) ? 1 : 0;
	}
	eigenvec.set_array(eigenarr);

	double eigenval_pre, eigenval;
	eigenval = 1;
	do {                                                       // calculate the eigenvalue using power iteration method, continue until the difference between consecutive approximations is within tolerance
		eigenval_pre = eigenval;
		eigenvec = A * eigenvec;
		eigenval = eigenvec.get_max();
		eigenvec = abs(1/eigenval) * eigenvec;
	} while (abs(eigenval_pre - eigenval) >= tolerance);

	return eigenval;
}

// hotelling's deflation algorithm
Matrix hotelling_deflation (Matrix A, Matrix vec, double eigenval) {
	int n = vec.get_m();
	Matrix I = id(n);
	Matrix trans = vec.transpose();
	double norm_square = (vec.transpose() * vec).get_arr()[0][0];
	return (A - (eigenval / norm_square) * (vec * vec.transpose()));
}

Matrix read_square_matrix_from_file (string filename) {
	ifstream fh;
	fh.open(filename);

	string line;
	int n;
	for (n = 0; getline(fh, line); n++);
	fh.clear();
	fh.seekg(0);

	double** arr = new double*[n];
	for (int i = 0; i < n; i++) {
		arr[i] = new double[n];
		for (int j = 0; j < n; j++) {
			fh >> arr[i][j];
		}
	}
	fh.close();

	Matrix A(n);
	A.set_array(arr);
	return A;
}


int main(int argc, char* argv[]) {
	int nec_no_args = 3;
	if (argc != nec_no_args + 1) {
		cout << "You should pass exactly " << nec_no_args << " arguments." << endl;
		return 1;
	}

	string A_file = argv[1];                                   // name of the input file
	double tolerance = stod(argv[2]);                          // tolerance level for the error in the eigenvalue
	string x_file = argv[3];                                   // name of the output file

	Matrix A = read_square_matrix_from_file(A_file);
	int n = A.get_n();
	Matrix eigenvec1(n, 1);
	double eigenval1 = normalized_power_iteration(A, eigenvec1, tolerance);

	Matrix B = hotelling_deflation(A, eigenvec1, eigenval1);
	Matrix eigenvec2(n, 1);
	double eigenval2 = normalized_power_iteration(B, eigenvec2, tolerance);

	cout << "Eigenvalue#1: " << eigenval1 << endl;
	eigenvec1.print(cout);
	cout << "Eigenvalue#2: " << eigenval2 << endl;

	ofstream fh;
	fh.open(x_file);
	fh << "Eigenvalue#1: " << eigenval1 << endl;
	eigenvec1.print(fh);
	fh << "Eigenvalue#2: " << eigenval2 << endl;
	fh.close();
}
