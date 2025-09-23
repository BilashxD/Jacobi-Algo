#include <stdio.h>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

// aDDING A COMMENT LINE

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



// do here number 3



void sortEigenvalues(vector<double>& eigenvalues, vector<vector<double>>& eigenvectors) {
        int n = eigenvalues.size();
        vector<pair<double,int>> pairs;
        for (int i = 0; i < n; ++i) pairs.emplace_back(eigenvalues[i], i);
        sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b){ return a.first > b.first; });
        vector<double> sorted_eigenvalues(n);
        vector<vector<double>> sorted_eigenvectors(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            sorted_eigenvalues[i] = pairs[i].first;
            int idx = pairs[i].second;
            for (int j = 0; j < n; ++j) sorted_eigenvectors[j][i] = eigenvectors[j][idx];
        }
        eigenvalues = sorted_eigenvalues;
        eigenvectors = sorted_eigenvectors;
    }

    void solve(const vector<vector<double>>& A, vector<double>& eigenvalues, vector<vector<double>>& eigenvectors) {
        int n = A.size();
        vector<vector<double>> A_work = A;
        eigenvectors = vector<vector<double>>(n, vector<double>(n,0.0));
        for (int i=0;i<n;i++) eigenvectors[i][i] = 1.0;

        int iteration = 0; double max_off_diag; int p,q;
        do {
            findMaxOffDiagonal(A_work, p, q, max_off_diag);
            if (abs(max_off_diag) > epsilon) {
                double c,s;
                computeRotation(A_work,p,q,c,s);
                applyRotation(A_work,eigenvectors,p,q,c,s);
            }
            iteration++;
        } while(abs(max_off_diag) > epsilon && iteration < max_iterations);

        eigenvalues.resize(n);
        for(int i=0;i<n;i++) eigenvalues[i] = A_work[i][i];
        sortEigenvalues(eigenvalues,eigenvectors);
    }



// Function to print eigenvectors where each column of the matrix is treated as an eigenvector
static void printEigenvectorsAsColumns(const vector<vector<double>>& V){
    int n = V.size();   // Get the size of the matrix (assuming square matrix: n x n)
    cout << "Eigenvectors (columns):\n";

    // Loop through each column (eigenvector)
    for(int col = 0; col < n; ++col){
        cout << "v[" << col << "] = [ ";
        
        // Loop through each row to print the current column values
        for(int row = 0; row < n; ++row){
            cout << setprecision(6) << V[row][col];   // Print element with 6-digit precision
            if(row < n-1) cout << ", ";              // Add a comma unless it's the last element
        }
        cout << " ]\n";  // End of current eigenvector
    }
    cout << "\n";  // Extra newline for readability
}

// Function to print a matrix with a given name
static void printMatrix(const vector<vector<double>>& M, const string& name){
    cout << name << ":\n";  // Print the name of the matrix

    // Loop through each row of the matrix
    for(auto &row : M){
        cout << "[ ";
        
        // Loop through each element in the row
        for(size_t j=0; j<row.size(); j++){
            cout << setprecision(6) << row[j];   // Print element with 6-digit precision
            if(j < row.size()-1) cout << ", ";  // Add a comma unless it's the last element
        }
        cout << " ]\n";  // End of row
    }
    cout << "\n";  // Extra newline after matrix
}

// Function to print a vector with a given name
static void printVector(const vector<double>& v, const string& name){
    cout << name << ": [";

    // Loop through each element of the vector
    for(size_t i=0; i<v.size(); i++){
        cout << setprecision(6) << v[i];   // Print element with 6-digit precision
        if(i < v.size()-1) cout << ", ";   // Add a comma unless it's the last element
    }
    cout << "]\n\n";  // Close vector and add extra newline
}






bool validateEigenSolution(const vector<vector<double>>& A, const vector<double>& eigenvalues, const vector<vector<double>>& V, double tol = 1e-8) {
    int n = A.size();
    vector<vector<double>> AV(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            for(int k=0;k<n;k++)
                AV[i][j] += A[i][k]*V[k][j];

    vector<vector<double>> VLambda(n, vector<double>(n,0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            VLambda[i][j] = V[i][j]*eigenvalues[j];

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            if(abs(AV[i][j]-VLambda[i][j])>tol) return false;

    return true;
}

bool validateOrthogonality(const vector<vector<double>>& V, double tol = 1e-8) {
    int n = V.size();
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            double dot = 0.0;
            for(int k=0;k<n;k++) dot += V[k][i]*V[k][j];
            if(abs(dot) > tol) return false;
        }
    }
    return true;
}

bool validateEigenvectors(const vector<vector<double>>& A, const vector<double>& eigenvalues, const vector<vector<double>>& V, double tol=1e-8) {
    int n = A.size();
    for(int j=0;j<n;j++){
        for(int i=0;i<n;i++){
            double Av_i = 0.0;
            for(int k=0;k<n;k++) Av_i += A[i][k]*V[k][j];
            if(abs(Av_i - eigenvalues[j]*V[i][j]) > tol) return false;
        }
    }
    return true;
}

bool isSymmetric(const vector<vector<double>>& A, double tol = 1e-8) {
    int n = A.size();
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            if(abs(A[i][j]-A[j][i])>tol) return false;
    return true;
}

















