// matrix_operations.h
#include <vector>
#include <algorithm>


class Matrix{
private:
    Matrix mult_slow(const Matrix& other) const;
    void switch_rows(int i, int j);
    void switch_cols(int i, int j);
    std::vector<double> data;

public:
    Matrix(int rows, int cols);
    int rows, cols;

    void display() const;
    double& at(int i, int j);
    const double& at(int i, int j) const;
    Matrix add(const Matrix& other) const;
    Matrix sub(const Matrix& other) const;
    Matrix mult(const Matrix& other) const;
    Matrix transpose() const;
    std::vector<Matrix> PLU_decomp() const;
    Matrix inverse() const;
    double determinant() const;
    Matrix power(int p);
};


