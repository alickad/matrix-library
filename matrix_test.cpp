#include <iostream>
#include <vector>
#include "matrix.h"
#include <random>
typedef long long ll;
using namespace std;

int main(){
    int m,k,n;
    cin >> m >> k >> n;    
    Matrix a(m, k);
    Matrix b(k, n);

    for (int mm = 0; mm<m; mm++){
        for (int kk = 0; kk<k; kk++){
            a.data[mm][kk] = double(rand() % 10);
        }
    }
    for (int kk = 0; kk<k; kk++){
        for (int nn = 0; nn<n; nn++){
            b.data[kk][nn] = double(rand() % 10);
        }
    }
    Matrix c = a.mult(b);
    Matrix d = a.mult_slow(b);
    for (int i = 0; i<c.rows; i++){
        for (int j = 0; j<c.cols; j++){
            cout << c.data[i][j] << " ";
            if (c.data[i][j] != d.data[i][j]) cout << "!!!!!!!!!!!!!!!!!!!!";
        }
        cout << '\n';
    }

    cout << '\n';
    Matrix x = a.inverse();
    x.display();
}