#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
using namespace std;

class Matrix {
public:

	void set_m(int _m) { m = _m;   }
	void set_n(int _n) { n = _n;   }
	int get_m()		   { return m; }
	int get_n()        { return n; }

	//fill the matrix through user input
	void fill_matrix();
	//fill the matrix randomly
	void fill_matrix_rand();

	//default and non default csonstructors
	Matrix(int _m, int _n);
	//copy constructor
	Matrix(const Matrix& rhs);
	//destructor
	~Matrix();
	//copy assignment operator
	Matrix& operator=(const Matrix& rhs);
	//move constructor
	Matrix(Matrix&& rhs) noexcept;
	//move assignment operator
	Matrix& operator=(Matrix&& rhs);

	//printing
	friend ostream& operator<<(ostream& outs, const Matrix& obj);

	//Matrix addition:
	Matrix operator+(const Matrix& rhs);

	//Scalar multiplication (2 functions for commutativity)
	Matrix operator*(int k);
	friend Matrix operator*(int k, Matrix obj);

	//Matrix Multiplication
	Matrix operator* (const Matrix& rhs);

	//Matrix negation
	Matrix operator-();

	//Matrix subtraction
	Matrix operator-(Matrix rhs);

	//Equality
	bool operator==(const Matrix& rhs);

	//swaps row of idx r1 with row of idx r2 in a matrix 
	void swap_rows(int r1, int r2);

	//divides every num in the row by the passed div_val
	void divide_row(int row_idx, double div_val);

	//row  i - row j
	void subtract_rows(int i, int j);
	//row i = row i - f*rowj (f is a factor)
	void subtract_rows(int i, int j, double factor);

	//row i + row j
	void add_rows(int i, int j);
	// Add row j to row i multiplied by a factor
	void add_rows(int i, int j, double factor);

	//Transpose
	Matrix transpose();

	//identity
	void identity(int m);

	//Make REF
	Matrix to_ref() const;

	//Make RREF
	Matrix to_rref() const;

	//Find the determinant: Expansion along the first row
	double find_det() const;

	//Find the inverse
	Matrix inverse() const;

	//Extract the matrix after expanding along the first row
	Matrix extract_matrix(int col) const;


protected:
	int m, n;
	double** matrix ;
};

//default and non default constructors
Matrix::Matrix(int _m = 0, int _n = 0) {
	m = _m;
	n = _n;
	matrix = new double* [m];
	for (int i = 0; i < m; i++) {
		matrix[i] = new double[n];
	}
}
//copy constructor
Matrix::Matrix(const Matrix& rhs) {
	m = rhs.m;
	n = rhs.n;
	matrix = new double* [m];

	for (int i = 0; i < m; i++) {
		matrix[i] = new double[n];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = rhs.matrix[i][j];
		}
	}
}
//destructor
Matrix::~Matrix() {
	for (int i = 0; i < m; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}
//copy assignment operator
Matrix& Matrix::operator=(const Matrix& rhs) {

	//check for self assignment
	if (this == &rhs) {
		return *this;
	}
	//destroy and reconstruct if needed
	if (m != rhs.m || n != rhs.n) {

		for (int i = 0; i < m; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;

		m = rhs.m;
		n = rhs.n;
		matrix = new double* [m];
		for (int i = 0; i < m; i++) {
			matrix[i] = new double[n];
		}
	}
	//deep copy
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = rhs.matrix[i][j];
		}
	}
	return *this;
}
//move constructor
Matrix::Matrix(Matrix&& rhs) noexcept : m(rhs.m), n(rhs.n), matrix(rhs.matrix) {
	rhs.m = 0;
	rhs.n = 0;
	rhs.matrix = nullptr;
}
//move assignment operator
Matrix& Matrix::operator=(Matrix&& rhs) {

	if (this != &rhs) {

		for (int i = 0; i < m; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;

		m = rhs.m;
		n = rhs.n;
		matrix = rhs.matrix;

		rhs.m = 0;
		rhs.n = 0;
		rhs.matrix = nullptr;
	}

	return *this;

}

//filling in
void Matrix::fill_matrix() {
	cout << "Enter the numbers going left to right and top down" << endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cin >> matrix[i][j];
		}
	}
}
void Matrix::fill_matrix_rand() {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = rand() % 10 + 1;
		}
	}
}
//printing
ostream& operator<<(ostream& outs, const Matrix& obj) {
	outs << std::fixed << std::setprecision(3); // Set fixed-point notation and precision

	for (int i = 0; i < obj.m; i++) {
		for (int j = 0; j < obj.n; j++) {
			double value = obj.matrix[i][j];

			if (std::fabs(value - round(value)) < 0.001) {
				outs << std::setw(8) << static_cast<int>(round(value)) << " ";
			}
			else {
				outs << std::setw(8) << value << " "; // Adjust the width as needed
			}
		}
		outs << '\n';
	}
	outs << '\n';

	return outs;
}
//Matrix addition:
Matrix Matrix::operator+(const Matrix& rhs) {

	//check if valid
	if (m != rhs.m || n != rhs.n) {
		cout << "Dimensions aren't compatible, operation cannot be done." << endl;
		return *this;
	}

	Matrix answer(this->m, this->n);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			answer.matrix[i][j] = matrix[i][j] + rhs.matrix[i][j];
		}
	}

	return answer;
}
//Scalar multiplication (2 functions for commutativity)
Matrix Matrix::operator*(int k) {
	Matrix answer(this->m, this->n);

	for (int i = 0; i < this->m; i++) {
		for (int j = 0; j < this->n; j++) {
			if (this->matrix[i][j] == 0) {
				answer.matrix[i][j] = 0;
			}
			else {
				answer.matrix[i][j] = this->matrix[i][j] * k;
			}
		}
	}

	return answer;

}
Matrix operator*(int k, Matrix obj) {

	Matrix answer(obj.m, obj.n);

	for (int i = 0; i < obj.m; i++) {
		for (int j = 0; j < obj.n; j++) {
			if (obj.matrix[i][j] == 0) {
				answer.matrix[i][j] = 0;
			}
			else {
				answer.matrix[i][j] = obj.matrix[i][j] * k;
			}
		}
	}
	return answer;
}
//Matrix Multiplication
Matrix Matrix::operator* (const Matrix& rhs) {

	//check for validity
	if (n != rhs.m) {
		cout << "Dimensions aren't compatible, operation cannot be done." << endl;
		return *this;
	}

	int ans_m = m;
	int ans_n = rhs.n;

	Matrix answer(ans_m, ans_n);

	for (int i = 0; i < ans_m; i++) {
		for (int j = 0; j < ans_n; j++) {
			double sum = 0;
			for (int r = 0; r < this->m; r++) {
				sum += this->matrix[i][r] * rhs.matrix[r][j];
			}
			answer.matrix[i][j] = sum;
		}
	}

	return answer;

}
//Matrix negation
Matrix Matrix::operator-() {
	Matrix answer;
	answer = *this * -1;
	return answer;
}
//Matrix subtraction
Matrix Matrix::operator-(Matrix rhs) {
	Matrix answer;
	answer = *this + -rhs;
	return answer;
}
//Equality
bool Matrix::operator==(const Matrix& rhs) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (matrix[i][j] != rhs.matrix[i][j]) {
				return false;
			}
		}
	}
	return true;
}
//swaps row of idx r1 with row of idx r2 in a matrix 
void Matrix::swap_rows(int r1, int r2) {
	if (r1 == r2) {
		return; // no need to swap if the rows are the same
	}

	double* temp_row = matrix[r1];
	matrix[r1] = matrix[r2];
	matrix[r2] = temp_row;
}
//divides every num in the row by the passed div_val
void Matrix::divide_row(int row_idx, double div_val) {
	for (int c = 0; c < n; c++) {
		matrix[row_idx][c] /= div_val;
	}
}
//row  i - row j
void Matrix::subtract_rows(int i, int j) {
	for (int c = 0; c < n; c++) {
		matrix[i][c] -= matrix[j][c];
	}
}
//row i = row i - f*rowj (f is a factor)
void Matrix::subtract_rows(int i, int j, double factor) {
	for (int c = 0; c < n; c++) {
		matrix[i][c] -= factor * matrix[j][c];
	}
}
//row i + row j
void Matrix::add_rows(int i, int j) {
	for (int c = 0; c < n; c++) {
		matrix[i][c] += matrix[j][c];
	}
}
// Add row j to row i multiplied by a factor
void Matrix::add_rows(int i, int j, double factor) {
	for (int c = 0; c < n; c++) {
		matrix[i][c] += factor * matrix[j][c];
	}
}
//Transpose
Matrix Matrix::transpose() {
	Matrix answer(n, m);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			answer.matrix[j][i] = matrix[i][j];
		}
	}
	return answer;
}
//identity
void Matrix::identity(int m) {
	Matrix identity_matrix(m, m);
	*this = identity_matrix;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			if (i == j) {
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
	}
}
//Make REF
Matrix Matrix::to_ref() const {
	Matrix ref_matrix(*this); // Create a copy of the original matrix

	int pivot_col = 0;

	for (int r = 0; r < ref_matrix.m; r++) {
		if (pivot_col >= ref_matrix.n) {
			break; // If you reach the final column, everything before that is 0
			// The pivot is in the final column, and so we're done with the row
		}

		int i = r;
		while (ref_matrix.matrix[i][pivot_col] == 0) {
			i++;
			// If the last element in the column is zero, the whole column is 0s,
			// so check the next column
			if (i == ref_matrix.m) {
				i = r;
				pivot_col++;
				if (ref_matrix.n == pivot_col) {
					break;
				}
			}
		}

		if (i != r) {
			ref_matrix.swap_rows(i, r);
		}

		// Make the pivot 1
		if (ref_matrix.matrix[r][pivot_col] != 1) {
			ref_matrix.divide_row(r, ref_matrix.matrix[r][pivot_col]);
		}

		// Make everything under the pivot 0
		for (int i = r + 1; i < ref_matrix.m; i++) {
			double factor = ref_matrix.matrix[i][pivot_col];
			ref_matrix.subtract_rows(i, r, factor);
		}
		pivot_col++;
	}

	return ref_matrix;
}
//Make RREF
Matrix Matrix::to_rref() const {
	Matrix ref_matrix(*this); // Create a copy of the original matrix

	int pivot_col = 0;

	for (int r = 0; r < ref_matrix.m; r++) {
		if (pivot_col >= ref_matrix.n) {
			break; // If you reach the final column, everything before that is 0
			// The pivot is in the final column, and so we're done with the row
		}

		int i = r;
		while (ref_matrix.matrix[i][pivot_col] == 0) {
			i++;
			// If the last element in the column is zero, the whole column is 0s,
			// so check the next column
			if (i == ref_matrix.m) {
				i = r;
				pivot_col++;
				if (ref_matrix.n == pivot_col) {
					break;
				}
			}
		}

		if (i != r) {
			ref_matrix.swap_rows(i, r);
		}

		// Make the pivot 1
		if (ref_matrix.matrix[r][pivot_col] != 1) {
			ref_matrix.divide_row(r, ref_matrix.matrix[r][pivot_col]);
		}

		// Make everything above and under the pivot 0
		for (int i = 0; i < ref_matrix.m; i++) {
			if (i != r) {
				double factor = ref_matrix.matrix[i][pivot_col];
				ref_matrix.subtract_rows(i, r, factor);
			}
		}
		pivot_col++;
	}

	return ref_matrix;
}
//Find the determinant: Expansion along the first row
double Matrix::find_det() const {

	//verify validity
	if (m != n) {
		cout << "Error: Determinant is only defined for square matrices\n";
		return 0;
	}

	double det = 0;
	if (m == 2) {
		det += matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]; //ad - bc
	}
	else {
		Matrix small_mat(m - 1, m - 1);
		for (int i = 0; i < n; i++) {
			//extract the smaller matrix
			small_mat = this->extract_matrix(i);
			det += pow(-1, i) * matrix[0][i] * small_mat.find_det();
		}
	}

	return det;
}
//Find the inverse
Matrix Matrix::inverse() const {

	if (this->find_det() == 0) {
		cout << "Error: Inverse does not exist\n";
		return *this;
	}

	if (m != n) {
		cout << "Error: Inverse is only defined for square matrices\n";
		return *this;
	}

	Matrix inverse;
	inverse.identity(m);

	Matrix ref_matrix(*this); // Create a copy of the original matrix

	int pivot_col = 0;

	for (int r = 0; r < ref_matrix.m; r++) {
		if (pivot_col >= ref_matrix.n) {
			break; // If you reach the final column, everything before that is 0
			// The pivot is in the final column, and so we're done with the row
		}

		int i = r;
		while (ref_matrix.matrix[i][pivot_col] == 0) {
			i++;
			// If the last element in the column is zero, the whole column is 0s,
			// so check the next column
			if (i == ref_matrix.m) {
				i = r;
				pivot_col++;
				if (ref_matrix.n == pivot_col) {
					break;
				}
			}
		}

		if (i != r) {
			ref_matrix.swap_rows(i, r);
			inverse.swap_rows(i, r);
		}

		// Make the pivot 1
		if (ref_matrix.matrix[r][pivot_col] != 1) {
			double divisor = ref_matrix.matrix[r][pivot_col];
			ref_matrix.divide_row(r, divisor);
			inverse.divide_row(r, divisor);
		}

		// Make everything under the pivot 0
		for (int i = 0; i < ref_matrix.m; i++) {
			if (i != r) {
				double factor = ref_matrix.matrix[i][pivot_col];
				ref_matrix.subtract_rows(i, r, factor);
				inverse.subtract_rows(i, r, factor);
			}
		}
		pivot_col++;
	}

	return inverse;
}
//Extract the matrix after expanding along the first row
Matrix Matrix::extract_matrix(int col) const {

	Matrix answer(m - 1, n - 1);

	for (int i = 1; i < m; i++) {
		for (int j = 0, dest_col = 0; j < n; j++) {
			if (j != col) {
				answer.matrix[i - 1][dest_col++] = matrix[i][j];
			}
		}
	}

	return answer;
}

int main() {

	srand(time(0));

	/*Matrix m1(2,2);
	Matrix m2(2,2);
	Matrix m3;
	Matrix m4;
	Matrix m5;
	Matrix m6;
	Matrix m7;

	m1.fill_matrix();
	m2.fill_matrix();
	m3 = m1 + m2;
	m4 = m1 * m2;
	m5 = 3 * m1;
	m6 = m1 * 3;
	m7 = -m2;

	cout << "Matrix 1 : \n" << m1;
	cout << "Matrix 2 : \n" << m2;
	cout << "Matrix 1 + 2 : \n" << m3;
	cout << "Matrix 1 * 2 : \n" << m4;
	cout << "3 * Matrix 1 : \n" << m5;
	cout << "Matrix 1 * 3 : \n" << m6;
	cout << "- Matrix 2 : \n" << m7;*/

	int num_rows, num_cols;
	char dummy;
	cout << "Enter the dimension of the matrix as m x n : ";
	cin >> num_rows >> dummy >> num_cols;

	Matrix original(num_rows, num_cols), toREF, toRREF, inverse;
	original.fill_matrix();

	//toREF = original.to_ref();
	toRREF = original.to_rref();
	toREF = original.to_ref();
	inverse = original.inverse();

	cout << "Original:\n" << original;
	cout << "After REF:\n" << toREF;
	cout << "After RREF:\n" << toRREF;
	cout << "Inverse:\n" << inverse;
	cout << "Determinant of original: " << original.find_det() << endl;
	cout << "Determinant of REF: " << toREF.find_det() << endl;
	cout << "Determinant of RREF: " << toRREF.find_det() << endl;


	return 0;
}