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

// Compute rotation values (cosine c and sine s) 
// for Jacobi eigenvalue method
void computeRotation(const vector<vector<double>>& A, int p, int q, double& c, double& s) {
    if (abs(A[p][q]) < epsilon) {  // If value is almost zero → no rotation needed
        c = 1.0; 
        s = 0.0; 
        return; 
    }
    // tau helps decide rotation direction
    double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
    double t = (tau >= 0) ? 1.0 / (tau + sqrt(1.0 + tau * tau)) : -1.0 / (-tau + sqrt(1.0 + tau * tau));
    c = 1.0 / sqrt(1.0 + t * t);   // cos value
    s = t * c;                     // sin value
}

// Apply the rotation on matrix A and update eigenvectors in V
void applyRotation(vector<vector<double>>& A, vector<vector<double>>& V, int p, int q, double c, double s) {
    int n = A.size();

    // Update diagonal and (p,q) elements
    double app = A[p][p], aqq = A[q][q], apq = A[p][q];
    A[p][p] = c*c*app - 2.0*c*s*apq + s*s*aqq;
    A[q][q] = s*s*app + 2.0*c*s*apq + c*c*aqq;
    A[p][q] = A[q][p] = 0.0;   // make off-diagonal zero

    // Update other entries of A
    for (int j = 0; j < n; ++j) {
        if (j != p && j != q) {
            double apj = A[p][j], aqj = A[q][j];
            A[p][j] = A[j][p] = c*apj - s*aqj;
            A[q][j] = A[j][q] = s*apj + c*aqj;
        }
    }

    // Update eigenvector matrix V
    for (int i = 0; i < n; ++i) {
        double vip = V[i][p], viq = V[i][q];
        V[i][p] = c*vip - s*viq;
        V[i][q] = s*vip + c*viq;
    }
}














