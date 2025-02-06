#include "matrix.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
using namespace std;

Matrix::Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

Matrix Matrix::add(const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Matrix dimensions must match for addition.");

    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] + other.data[i][j];

    return result;
}

Matrix Matrix::sub(const Matrix& other) const{
    if (rows != other.rows || cols != other.cols)
        throw std::runtime_error("Matrix dimensions must match for subtraction.");

    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result.data[i][j] = data[i][j] - other.data[i][j];

    return result;
}

Matrix Matrix::mult_slow(const Matrix& other) const{
    if (cols != other.rows)
        throw std::runtime_error("Matrix dimensions must match for slow multiplication.");

    Matrix result(rows, other.cols);
    for (int m = 0; m<rows; m++){
        for (int n = 0; n<other.cols; n++){
            result.data[m][n] = 0;
            for (int k = 0; k<cols; k++){
                result.data[m][n] += data[m][k] * other.data[k][n];
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
        return mult_slow(other);            ///////////////////////////////////////////
    }

    int M = M0 + M0%2;
    int N = N0 + N0%2;
    int K = K0 + K0%2;

    Matrix A11(M/2, K/2);
    Matrix A12(M/2, K/2);
    Matrix A21(M/2, K/2);
    Matrix A22(M/2, K/2);
//TODO check if initializes to zero

    for (int m = 0; m<M0; m++){
        for (int k = 0; k<K0; k++){
            if (m<M/2 && k<K/2){
                A11.data[m][k] = data[m][k];
            }
            else if (m<M/2){
                A12.data[m][k-K/2] = data[m][k];
            }
            else if (k<K/2){
                A21.data[m-M/2][k] = data[m][k];
            }
            else{
                A22.data[m-M/2][k-K/2] = data[m][k];
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
                B11.data[k][n] = other.data[k][n];
            }
            else if (k<K/2){
                B12.data[k][n-N/2] = other.data[k][n];
            }
            else if (n<N/2){
                B21.data[k-K/2][n] = other.data[k][n];
            }
            else{
                B22.data[k-K/2][n-N/2] = other.data[k][n];
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
                C.data[m][n] = p1.data[m][n] + p4.data[m][n] - p5.data[m][n] + p7.data[m][n];
            }
            else if (n < N/2){
                C.data[m][n] = p2.data[m-M/2][n] + p4.data[m-M/2][n];
            }
            else if (m < M/2){
                C.data[m][n] = p3.data[m][n-N/2] + p5.data[m][n-N/2];
            }
            else{
                C.data[m][n] = p1.data[m-M/2][n-N/2] + p3.data[m-M/2][n-N/2] - p2.data[m-M/2][n-N/2] - p6.data[m-M/2][n-N/2];
            }
        }
    }

    return C;
}

Matrix Matrix::transpose() const{
    Matrix result(cols, rows);
    for (int y = 0; y<rows; y++){
        for (int x = 0; x < cols; x++){
            result.data[x][y] = data[y][x];
        }
    }
    return result;
}

Matrix Matrix::P_inverse() const{
    Matrix Pinv(rows, cols);
    for (int y = 0; y < rows; y++){
        for (int x = 0; x < cols; x++){
            if (data[y][x] == 1){
                Pinv.data[x][y] = 1;
            }
        }
    }
    return Pinv;
}

void Matrix::display() const{
    for (int y = 0; y < rows; y++){
        for (int x = 0; x < cols; x++){
            cout << data[y][x];
            if (x < cols-1) cout << " ";
            else cout << '\n';
        }
    }
}

void Matrix::switch_rows(int i, int j){
    for (int y = 0; y < cols; y++){
        swap(data[i][y], data[j][y]);
    }
}
void Matrix::switch_cols(int i, int j){
    for (int x = 0; x < rows; x++){
        swap(data[x][i], data[x][j]);
    }
}

vector<Matrix> Matrix::PLU_decomp() const{
    Matrix P(rows,rows);
    for (int i = 0; i < rows; i++)
        P.data[i][i] = 1;
    Matrix L(rows,rows);
    for (int i = 0; i < rows; i++)
        L.data[i][i] = 1;
    Matrix U(rows,cols);
    U.data = data;
    int row_now = 0;
    for (int x = 0; x < cols; x++){
        int max_row = row_now;
        for (int y = row_now; y < rows; y++){
            if (abs(U.data[y][x]) > abs(U.data[max_row][x])){
                max_row = y;
            }
        }
        if (U.data[max_row][x] == 0) continue;
        U.switch_rows(row_now, max_row);
        L.switch_rows(row_now, max_row);
        L.switch_cols(row_now, max_row);
        P.switch_cols(row_now, max_row);
        for (int y = row_now + 1; y < rows; y++){
            double pomer = U.data[y][x] / U.data[row_now][x];
            U.data[y][x] = 0;
            for (int xx = x+1; xx<cols; xx++){
                U.data[y][xx] -= U.data[row_now][xx] * pomer;
            }
            L.data[y][x] = pomer;
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
        invA.data[y][y] = 1;
    }
    for (int x = 0; x < cols; x++){
        int maxi = x;
        for (int y = x; y < rows; y++){
            if (abs(A.data[y][x]) > abs(A.data[maxi][x]))
                maxi = y;
        }
        if (A.data[maxi][x] == 0)
            throw std::runtime_error("Matrix is not invertible.");
        A.switch_rows(x, maxi);
        invA.switch_rows(x, maxi);
        // divide the row by pivot value
        for (int xx = x+1; xx<cols; xx++){
            A.data[x][xx] /= A.data[x][x];
        }
        for (int xx = 0; xx<cols; xx++){
            invA.data[x][xx] /= A.data[x][x];
        }
        A.data[x][x] = 1;
        // make the column except pivot zero 
        for (int y = 0; y < rows; y++){
            if (y == x) continue;
            double pomer = A.data[y][x];
            for (int xx = 0; xx < cols; xx++){
                A.data[y][xx] -= A.data[x][xx] * pomer;
                invA.data[y][xx] -= invA.data[x][xx] * pomer;
            }
            A.data[y][x] = 0;
        }        
    }

    return invA;
}