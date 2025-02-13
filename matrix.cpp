#include "matrix.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
using namespace std;

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(r * c, 0.0) {}

double& Matrix::at(int i, int j){
    if (i >= rows || i < 0 || j < 0 || j >= cols)
        throw std::runtime_error("Indexing in matrix out of range");
    return data[i*cols + j];
}
const double& Matrix::at(int i, int j) const {
    if (i >= rows || i < 0 || j < 0 || j >= cols)
        throw std::runtime_error("Indexing in matrix out of range");
    return data[i * cols + j];  
}

Matrix Matrix::add(const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Matrix dimensions must match for addition.");

    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.at(i,j) = at(i,j) + other.at(i,j);

    return result;
}

Matrix Matrix::sub(const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Matrix dimensions must match for subtraction.");

    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.at(i,j) = at(i,j) - other.at(i,j);

    return result;
}

Matrix Matrix::mult_slow(const Matrix& other) const{
    if (cols != other.rows)
        throw std::runtime_error("Matrix dimensions must match for slow multiplication.");

    Matrix result(rows, other.cols);
    for (int m = 0; m<rows; m++){
        for (int n = 0; n<other.cols; n++){
            result.at(m,n) = 0;
            for (int k = 0; k<cols; k++){
                result.at(m,n) += at(m,k) * other.at(k,n);
            }
        }
    }

    return result;
}

Matrix Matrix::mult(const Matrix& other) const{ //A, B
    int M0 = rows;
    int K0 = cols;
    int N0 = other.cols;
    if (cols != other.rows)
        throw std::runtime_error("Matrix dimensions must match for multiplication.");

    if (M0 < 30 && K0<30 && N0<30){
        return mult_slow(other);          
    }

    int M = M0 + M0%2;
    int N = N0 + N0%2;
    int K = K0 + K0%2;

    Matrix A11(M/2, K/2);
    Matrix A12(M/2, K/2);
    Matrix A21(M/2, K/2);
    Matrix A22(M/2, K/2);

    for (int m = 0; m<M0; m++){
        for (int k = 0; k<K0; k++){
            if (m<M/2 && k<K/2){
                A11.at(m,k) = at(m,k);
            }
            else if (m<M/2){
                A12.at(m,k-K/2) = at(m, k);
            }
            else if (k<K/2){
                A21.at(m-M/2, k) = at(m, k);
            }
            else{
                A22.at(m-M/2, k-K/2) = at(m, k);
            }
        }
    }

    Matrix B11(K/2, N/2);
    Matrix B12(K/2, N/2);
    Matrix B21(K/2, N/2);
    Matrix B22(K/2, N/2);
// TODO check, if initializes to zero
    for (int k = 0; k<K0; k++){
        for (int n = 0; n<N0; n++){
            if (k<K/2 && n<N/2){
                B11.at(k, n) = other.at(k, n);
            }
            else if (k<K/2){
                B12.at(k, n-N/2) = other.at(k, n);
            }
            else if (n<N/2){
                B21.at(k-K/2, n) = other.at(k, n);
            }
            else{
                B22.at(k-K/2, n-N/2) = other.at(k, n);
            }
        }
    }
    
    Matrix p1 = (A11.add(A22)).mult(B11.add(B22)); //1
    Matrix p2 = (A21.add(A22)).mult(B11); //6
    Matrix p3 = A11.mult(B12.sub(B22)); //5
    Matrix p4 = A22.mult(B21.sub(B11)); //2
    Matrix p5 = (A11.add(A12)).mult(B22); //3
    Matrix p6 = (A11.sub(A21)).mult(B11.add(B12)); //7
    Matrix p7 = (A12.sub(A22)).mult(B21.add(B22)); //4

    Matrix C(M0, N0);
    for (int m = 0; m<M0; m++){
        for (int n = 0; n<N0; n++){
            if (m < M/2 && n<N/2){
                C.at(m, n) = p1.at(m, n) + p4.at(m, n) - p5.at(m, n) + p7.at(m, n);
            }
            else if (n < N/2){
                C.at(m, n) = p2.at(m-M/2, n) + p4.at(m-M/2, n);
            }
            else if (m < M/2){
                C.at(m, n) = p3.at(m, n-N/2) + p5.at(m, n-N/2);
            }
            else{
                C.at(m, n) = p1.at(m-M/2, n-N/2) + p3.at(m-M/2, n-N/2) - p2.at(m-M/2, n-N/2) - p6.at(m-M/2, n-N/2);
            }
        }
    }

    return C;
}

Matrix Matrix::transpose() const{
    Matrix result(cols, rows);
    for (int y = 0; y<rows; y++){
        for (int x = 0; x < cols; x++){
            result.at(x, y) = at(y, x);
        }
    }
    return result;
}

void Matrix::display() const{
    for (int y = 0; y < rows; y++){
        for (int x = 0; x < cols; x++){
            cout << at(y, x);
            if (x < cols-1) cout << " ";
            else cout << '\n';
        }
    }
}

void Matrix::switch_rows(int i, int j){
    for (int y = 0; y < cols; y++){
        swap(at(i, y), at(j, y));
    }
}
void Matrix::switch_cols(int i, int j){
    for (int x = 0; x < rows; x++){
        swap(at(x, i), at(x, j));
    }
}

vector<Matrix> Matrix::PLU_decomp() const{
    Matrix P(rows,rows);
    for (int i = 0; i < rows; i++)
        P.at(i, i) = 1;
    Matrix L(rows,rows);
    for (int i = 0; i < rows; i++)
        L.at(i, i) = 1;
    Matrix U(rows,cols);
    U.data = data;
    int row_now = 0;
    for (int x = 0; x < cols; x++){
        int max_row = row_now;
        if (row_now >= rows) break;
        for (int y = row_now; y < rows; y++){
            if (abs(U.at(y, x)) > abs(U.at(max_row, x))){
                max_row = y;
            }
        }
        if (U.at(max_row, x) == 0) continue;
        U.switch_rows(row_now, max_row);
        L.switch_rows(row_now, max_row);
        L.switch_cols(row_now, max_row);
        P.switch_cols(row_now, max_row);
        for (int y = row_now + 1; y < rows; y++){
            double pomer = U.at(y, x) / U.at(row_now, x);
            U.at(y, x) = 0;
            for (int xx = x+1; xx<cols; xx++){
                U.at(y, xx) -= U.at(row_now, xx) * pomer;
            }
            L.at(y,x) = pomer;
        }
        row_now++;
    }

    return {P,L,U};
}

Matrix Matrix::inverse() const{
    if (rows != cols)
        throw std::runtime_error("Matrix must be square for inversion");
    Matrix A(rows,cols);
    A.data = data;
    Matrix invA(rows, cols);
    for (int y = 0; y < rows; y++){
        invA.at(y,y) = 1;
    }
    for (int x = 0; x < cols; x++){
        int maxi = x;
        for (int y = x; y < rows; y++){
            if (abs(A.at(y,x)) > abs(A.at(maxi,x)))
                maxi = y;
        }
        if (A.at(maxi, x) == 0)
            throw std::runtime_error("Matrix is not invertible.");
        A.switch_rows(x, maxi);
        invA.switch_rows(x, maxi);
        // divide the row by pivot value
        for (int xx = x+1; xx<cols; xx++){
            A.at(x, xx) /= A.at(x, x);
        }
        for (int xx = 0; xx<cols; xx++){
            invA.at(x,xx) /= A.at(x,x);
        }
        A.at(x,x) = 1;
        // make the column except pivot zero 
        for (int y = 0; y < rows; y++){
            if (y == x) continue;
            double pomer = A.at(y,x);
            for (int xx = 0; xx < cols; xx++){
                A.at(y,xx) -= A.at(x,xx) * pomer;
                invA.at(y,xx) -= invA.at(x,xx) * pomer;
            }
            A.at(y,x) = 0;
        }        
    }

    return invA;
}

double Matrix::determinant() const{
    if (rows != cols)
        throw runtime_error("Determinant is defined only on square matrices");
    vector<Matrix> PLU = PLU_decomp();
    Matrix P = PLU[0];
    Matrix L = PLU[1];
    Matrix U = PLU[2];

    double det = 1;

    // first determinant of P (using number of cycles of permutation)
    vector<int> perm(rows,0);
    for (int y = 0; y < rows; y++){
        for (int x = 0; x < cols; x++){
            if (P.at(y,x)){
                perm[x] = y;
                break;
            }
        }
    }
    vector<bool> visited(rows, false);
    for (int i = 0; i < rows; i++){
        if (visited[i]) continue;
        det *= -1;
        visited[i] = true;
        int j = perm[i];
        while (!visited[j]){
            visited[j] = true;
            j = perm[j];
        }
    }
    if (rows % 2 == 1) det *= -1;

    // multiply by determinant of L
    for (int i = 0; i<rows; i++){
        det *= L.at(i,i);
    }
    // multiply by determinant of U
    for (int i = 0; i<rows; i++){
        det *= U.at(i,i);
    }
    
    return det;
}

Matrix Matrix::power(int p){
    if (rows != cols) 
        throw runtime_error("The matrix must be square to be multiplied by itself");
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++){
        result.at(i,i) = 1;
    }
    Matrix A(rows, cols);
    A.data = data;
    if (p < 0){
        A = A.inverse();
        p = -p;
    }
    for (int i = 0; (1 << i) <= p; i++){
        if ((p >> i) & 1)
            result = result.mult(A);
        A = A.mult(A);
    }
    return result;
}
