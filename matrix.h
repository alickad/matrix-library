// matrix_operations.h
#include <vector>
#include <algorithm>


class Matrix{
private:
    

public:
    Matrix(int rows, int cols);

    std::vector<std::vector<double>> data;
    int rows, cols;
    Matrix mult_slow(const Matrix& other) const;
    void display() const;
    Matrix add(const Matrix& other) const;
    Matrix sub(const Matrix& other) const;
    Matrix mult(const Matrix& other) const;
    Matrix transpose() const;
    Matrix P_inverse() const;
    void switch_rows(int i, int j);
    void switch_cols(int i, int j);
    std::vector<Matrix> PLU_decomp() const;
    Matrix inverse() const;
};


