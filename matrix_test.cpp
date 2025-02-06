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
    
    // test PLU decomp
    vector<Matrix> PLU_test = a.PLU_decomp();
    a.display();
    cout << '\n';
    PLU_test[0].display();
    cout << '\n';
    PLU_test[1].display();
    cout << '\n';
    PLU_test[2].display();
    cout << '\n' << "a ich sucin je \n";
    ((PLU_test[0].mult(PLU_test[1])).mult(PLU_test[2])).display();



    
}