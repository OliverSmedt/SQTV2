//This code implements a complex matrix class
#ifndef COMPLEX_MATRIX_H
#define COMPLEX_MATRIX_H
#include<iostream>
#include<complex>
#include<iomanip>
#include<sstream>
#include<regex>

class matrix
{
    // Friends
    friend std::ostream& operator<<(std::ostream& os, const matrix& mat);
    friend matrix hstack(matrix mat1, matrix mat2);
    friend matrix vstack(matrix mat1, matrix mat2);
private:
    std::complex<double>* matrix_data{ nullptr };

    int rows{ 0 };
    int columns{ 0 };
    int size{ 0 };
public:
    // Default constructor
    matrix() = default;
    // Parameterized constructor
    matrix(int m, int n) :
        rows{ m }, columns{ n }
    {
        if (n < 0 || m < 0) { std::cout << "dimensions must be positive" << std::endl; exit(1); }
        else
        {
            size = m * n;
            matrix_data = new std::complex<double>[size];
        }
    }

    //deep copy
    matrix(matrix& original)
    {
        rows = original.rows;
        columns = original.columns;
        size = original.size;

        matrix_data = new std::complex<double>[size];
        for (int i = 0; i < (original.rows * original.columns); i++) {
            matrix_data[i] = original.matrix_data[i];
        }
    }
    // Move constructor
    matrix(matrix&& original) noexcept
    {
        columns = original.columns;
        rows = original.rows;
        size = original.size;
        matrix_data = original.matrix_data;
        original.columns = 0;
        original.rows = 0;
        original.matrix_data = nullptr;
    }

    // Destructor
    ~matrix() { delete[] matrix_data; }// std::cout << "matrix deleted" << std::endl;

    //unasigns a matrix
    void clear() {
        delete[] matrix_data;
        matrix_data = nullptr;
        rows = 0;
        columns = 0;
        size = 0;
    }

    //get the probabilities
    matrix get_modulus()
    {
        matrix new_matrix{ rows, columns };
        for (int index{ 0 }; index < size; index++) {
            double a = 0;
            new_matrix.matrix_data[index] = std::pow(std::abs(matrix_data[index]), 2);
        }
        return new_matrix;
    }

    // Access functions
    int get_rows() const { return rows; } // Return number of rows
    int get_cols() const { return columns; } // Return number of columns
    int get_size() const { return size; }
    int index(int m, int n) const // Return position in array of element (m,n)
    {
        if (m > 0 && m <= rows && n > 0 && n <= columns) return (n - 1) + (m - 1) * columns;
        else { std::cout << "Error: out of range" << std::endl; exit(1); }
    }
    std::complex<double>& operator()(int m, int n) const { return matrix_data[index(m, n)]; } //This function can be an Lvalue and hence can be used to assign elements of matrix

    //Copy  Assignment operator
    matrix& operator=(matrix& original);

    //Move Assignment operator
    matrix& operator=(matrix&& original) noexcept;

    //Addition, subtraction
    matrix operator+(const matrix& mat);
    matrix operator-(const matrix& mat);
    matrix& operator+=(const matrix& mat);
    matrix& operator-=(const matrix& mat);

    //multiplication overloads 
    matrix operator*(const matrix& mat);
    matrix operator*(const std::complex<double>& factor);
    matrix operator*(const double& factor);
    matrix operator*(const int& factor);
    matrix& operator*=(const matrix& mat);
    matrix& operator*=(const std::complex<double>& factor);
    matrix& operator*=(const double& factor);
    matrix& operator*=(const int& factor);

    //division overloads
    matrix operator/(const std::complex<double>& factor);
    matrix operator/(const double& factor);
    matrix operator/(const int& factor);
    matrix& operator/=(const std::complex<double>& factor);
    matrix& operator/=(const double& factor);
    matrix& operator/=(const int& factor);

};

// Member functions defined outside class
matrix& matrix::operator=(matrix& original) {
    if (&original == this) return *this; //No self assignemnt
    delete[] matrix_data; matrix_data = nullptr;
    columns = original.get_cols();
    rows = original.get_rows();
    size = columns * rows;

    if (size > 0) { //load in ne matrix_data
        matrix_data = new std::complex<double>[size];
        for (int i = 0; i < size; i++) {
            matrix_data[i] = original.matrix_data[i];
        }
    }
    return *this;
}
matrix& matrix::operator=(matrix&& original) noexcept
{
    if (&original == this) return *this; //No self assignment
    std::swap(columns, original.columns); //swap data
    std::swap(rows, original.rows);
    std::swap(size, original.size);
    std::swap(matrix_data, original.matrix_data);

    original.columns = 0; //unassign original matrix
    original.rows = 0;
    original.size = 0;
    original.matrix_data = nullptr;
    return *this;
}

//addition, subtraction
matrix matrix::operator+(const matrix& mat) {
    if (this->rows != mat.rows || this->columns != mat.columns) {
        std::cout << "Error: non-matching matrix sizes" << std::endl;
        matrix unasigned_matrix;
        return unasigned_matrix;
    }
    else {
        matrix added{ this->rows,this->columns };
        for (int i = 0; i < mat.rows * mat.columns; i++) {
            added.matrix_data[i] = this->matrix_data[i] + mat.matrix_data[i]; // simply add the contents of the individual sets of matrix_data
        }
        return added;
    }
}
matrix matrix::operator-(const matrix& mat) {
    if (this->rows != mat.rows || this->columns != mat.columns) {
        std::cout << "Error: non-matching matrix sizes" << std::endl;
        matrix unasigned_matrix;
        return unasigned_matrix;
    }
    else {
        matrix subtracted{ this->rows,this->columns };


        for (int i = 0; i < mat.size; i++) {
            subtracted.matrix_data[i] = this->matrix_data[i] - mat.matrix_data[i];
        }

        return subtracted;
    }
}
matrix& matrix::operator+=(const matrix& mat) {
    if (this->rows != mat.rows || this->columns != mat.columns) {
        std::cout << "Error: non-matching matrix sizes" << std::endl;
        this->clear();
        return *this;
    }
    else {
        for (int i = 0; i < mat.rows * mat.columns; i++) {
            this->matrix_data[i] += mat.matrix_data[i];
        }
        return *this;
    }
}
matrix& matrix::operator-=(const matrix& mat) {
    if (this->rows != mat.rows || this->columns != mat.columns) {
        std::cout << "Error: non-matching matrix sizes" << std::endl;
        this->clear();
        return *this;
    }
    else {
        for (int i = 0; i < mat.rows * mat.columns; i++) {
            this->matrix_data[i] -= mat.matrix_data[i];
        }
        return *this;
    }
}

//multiplication
matrix matrix::operator*(const matrix& mat) {
    if (this->columns != mat.rows) {
        std::cout << "Error: invalid matrix sizes" << std::endl;
        matrix unasigned_matrix;
        return unasigned_matrix;

    }
    else {
        matrix multi{ this->rows,mat.columns };
        for (int i = 1; i <= this->rows; i++) {
            for (int j = 1; j <= mat.columns; j++) {
                for (int k = 1; k <= mat.rows; k++) {
                    multi(i, j) += this->operator()(i, k) * mat(k, j);
                }

            }
        }
        return multi;
    }
}
matrix matrix::operator*(const std::complex<double>& factor) {
    matrix multi = *this;
    for (int index{ 0 }; index < size; index++) {
        multi.matrix_data[index] *= factor;
    }
    return multi;
}
matrix matrix::operator*(const double& factor) {
    matrix multi = *this;
    for (int index{ 0 }; index < size; index++) {
        multi.matrix_data[index] *= factor;
    }
    return multi;
}
matrix matrix::operator*(const int& factor) {
    matrix multi = *this;
    for (int index{ 0 }; index < size; index++) {
        multi.matrix_data[index] *= factor;
    }
    return multi;
}
matrix& matrix::operator*=(const matrix& mat) {
    if (this->columns != mat.rows) {
        std::cout << "Error: invalid matrix sizes" << std::endl;
        matrix unasigned_matrix;
        *this = unasigned_matrix;
        return *this;

    }
    else {
        matrix multi{ this->rows,mat.columns };
        //std::complex<double> sum{ 0, 0 };
        for (int i = 1; i <= this->rows; i++) {
            for (int j = 1; j <= mat.columns; j++) {
                for (int k = 1; k <= mat.rows; k++) {
                    multi(i, j) += this->operator()(i, k) * mat(k, j); //simple application of usual matrix multiplication formula
                }
                //multi(i, j) = sum;
                //sum.imag = 0;
                //sum.real = 0;
            }
        }
        *this = multi;
        return *this;
    }
}
matrix& matrix::operator*=(const std::complex<double>& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] *= factor;
    }
    return *this;
}
matrix& matrix::operator*=(const double& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] *= factor;
    }
    return *this;
}
matrix& matrix::operator*=(const int& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] *= factor;
    }
    return *this;
}

//division by scalar
matrix matrix::operator/(const std::complex<double>& factor) {
    matrix div = *this;
    for (int index{ 0 }; index < size; index++) {
        div.matrix_data[index] /= factor;
    }
    return div;
}
matrix matrix::operator/(const double& factor) {
    matrix div = *this;
    for (int index{ 0 }; index < size; index++) {
        div.matrix_data[index] /= factor;
    }
    return div;
}
matrix matrix::operator/(const int& factor) {
    matrix div = *this;
    for (int index{ 0 }; index < size; index++) {
        div.matrix_data[index] /= factor;
    }
    return div;
}
matrix& matrix::operator/=(const std::complex<double>& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] /= factor;
    }
    return *this;
}
matrix& matrix::operator/=(const double& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] /= factor;
    }
    return *this;
}
matrix& matrix::operator/=(const int& factor) {
    for (int index{ 0 }; index < size; index++) {
        matrix_data[index] /= factor;
    }
    return *this;
}

// Overload insertion operator
std::ostream& operator<<(std::ostream& os, const matrix& mat) //this code does not put each row on a new line to keep it compatible with the extraction operator
{
    if (mat.size == 0) return os;

    for (int i = 1; i <= mat.rows; i++) {

        for (int j = 1; j <= mat.columns; j++) {
            os << std::setprecision(3) << mat(i, j) << " ";
            if (j == mat.columns);
            else;
        }
        os << '\n';
    }

    return os;
}

matrix hstack(matrix mat1, matrix mat2) {
    if (mat1.get_size() == 0) return mat2; //This is so you can append to an unassigned matrix
    else if (mat2.get_size() == 0) return mat1;
    matrix result{ mat1.rows,mat1.columns + mat2.columns };
    for (int i = 1; i <= result.rows; i++) {
        for (int j = 1; j <= result.columns; j++) {
            if (j <= mat1.columns) result(i, j) = mat1(i, j);
            else result(i, j) = mat2(i, j - mat1.columns);
        }
    }
    return result;
}

matrix vstack(matrix mat1, matrix mat2) {
    if (mat1.get_size() == 0) return mat2;
    else if (mat2.get_size() == 0) return mat1;
    matrix result{ mat1.rows + mat2.rows,mat1.columns };
    for (int i = 1; i <= result.rows; i++) {
        for (int j = 1; j <= result.columns; j++) {
            if (i <= mat1.rows) result(i, j) = mat1(i, j);
            else result(i, j) = mat2(i - mat1.rows, j);

        }
    }
    return result;
}

matrix tensor_product(matrix mat1, matrix mat2) {
    if (mat1.get_size() == 0) return mat2; //return other matrix if one of them is zero
    else if (mat2.get_size() == 0) return mat1; //this can be useful for iterative statements (e.g. tensor product of parallel gates)
    matrix result{};
    matrix temp_row{};
    for (int i = 1; i <= mat1.get_rows(); i++) {
        for (int j = 1; j <= mat1.get_cols(); j++) {
            temp_row = hstack(temp_row, mat2 * mat1(i, j));
        }
        result = vstack(result, temp_row);
        temp_row.clear();
    }

    return result;
}
#endif