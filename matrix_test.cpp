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
    
    // test determinant
    a.display();
    cout << a.determinant() << '\n';



    
}