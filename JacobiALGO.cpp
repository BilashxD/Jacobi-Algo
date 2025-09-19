#include <stdio.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

class JacobiEigenSolver{
private: 
    const double epsilon;
    const int max_iterations;

public: 
     JacobiEigenSolver(double eps = 1e-12, int max_iter = 1000)
     :epsilon(eps),
     max_iterations(max_iter){}

     void findMaxOfDiagonal(const vector<vector<double>>& A, int& p, int& q, double& max_val){
      max_val = 0.0;
      int n = A.size();
      for(int i=0; i<n; ++i)
      for(int j=i+1; j<n; ++j)
      if(abs(A[i][j])>max_val) {

        max_val = abs(A[i][j]);
        p=i;
        q=j;
      }

     }

};

// 2 start from here (down)




//3 start from here




//4 start from here




//5 start from here