# Matrix-library
A simple library with matrix operations in C++.

## Setup
1. Clone this repository or download the files:
```
git clone <your-repo-url>
```
2. Ensure you have a C++ compiler (e.g., g++ for GCC).
3. Compile the project using:
```
g++ -o matrix_program matrix.cpp main.cpp
./matrix_program
```

Now to use the library, simply include "matrix.h" in your file.

## Matrix class
The class contains data about your matrix, as well as functions you can use on it.
### data
| rows  | number of rows of matrix             |
| cols  | number of columns of matrix          |
| data  | the values in a 2D vector of doubles |
### functions
|add(Matrix other)  |	Returns sum of two matrices                     |
|sub(Matrix other)  |	Returns difference of two matrices              |
|mult(Matrix other)	|   Returns product of two matrices                 |
|transpose()	    |   Returns the transposed matrix                   |
|PLU_decomp()	    |   Performs PLU decomposition                      |
|determinant()      |	Computes the determinant using PLU decomposition|
|inverse()	        |   Computes the inverse if the matrix is invertible|

### code example
```c++
#include "matrix.h"
#include <iostream>

int main() {
    Matrix A(3, 3); // specify number of rows and columns
    Matrix B(3, 3);
    
    // Example: Setting values
    A.data = {{2, 1, 3}, {4, 1, 6}, {7, 8, 9}};
    B.data = {{3, 4, 1}, {2, 5, 7}, {3, 2, 3}};
    
    // Compute determinant
    std::cout << "Determinant: " << A.determinant() << std::endl;

    // Find the product of A and B
    std::cout << "Product: " << A.mult(B) << std::endl;

    return 0;
}

```

## Used algorithms and their asymptotic complexity
**add, sub** simply add or subtract numbers from matrices one by one. For Matrix $n \times n$ this gives $O(n^2)$ complexity.
**transpose** just transpose in $O(n^2)$.
**mult** for smaller matrices, $O(n^3)$$. For big matrices we use Strassen`s algorithm, which gives us $$~O(n^2.81)$$.
**PLU_decomp** calculates 3 matrices, P, L and U. Their product is the given matrix. P is a permutation matrix, L is a lower triangular matrix and U is a rectangular matrix with all non-zero values on or above diagonal. We compute this using Gaussian elimination with partial pivoting. That gives us $O(n^3)$ time complexity for a $n \times n$ matrix.
**determinant** first calculate the PLU decomposition. Then calculate the determinants of the 3 matrices and multiply them. That gives us $O(n^3)$ time complexity.
**inverse** using Gauss-Jordan elimination with partial pivoting in $O(n^3)$.

