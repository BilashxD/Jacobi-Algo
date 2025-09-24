#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

// This class implements the Jacobi method to compute eigenvalues and eigenvectors of a symmetric matrix
class JacobiEigenSolver {
private:
    const double epsilon;       // Small tolerance value to determine when off-diagonal elements are effectively zero
    const int max_iterations;   // Maximum number of iterations to prevent infinite loops

public:
    // Constructor allows custom tolerance and maximum iteration values
    JacobiEigenSolver(double eps = 1e-12, int max_iter = 1000)
        : epsilon(eps), max_iterations(max_iter) {}

    // This function finds the largest off-diagonal element in the matrix
    // p and q will store the row and column of that element
    void findMaxOffDiagonal(const vector<vector<double>>& A, int& p, int& q, double& max_val) {
        max_val = 0.0;          // Start with zero for comparison
        int n = A.size();
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (abs(A[i][j]) > max_val) {  // Compare absolute values
                    max_val = abs(A[i][j]);
                    p = i;
                    q = j;
                }
    }

    // Compute the cosine and sine values for the Jacobi rotation
    void computeRotation(const vector<vector<double>>& A, int p, int q, double& c, double& s) {
        if (abs(A[p][q]) < epsilon) { c = 1.0; s = 0.0; return; } // If element is already very small, skip rotation
        double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
        double t = (tau >= 0) ? 1.0 / (tau + sqrt(1.0 + tau*tau)) : -1.0 / (-tau + sqrt(1.0 + tau*tau));
        c = 1.0 / sqrt(1.0 + t*t);  // cosine of rotation
        s = t * c;                  // sine of rotation
    }

    // Apply the rotation to zero out the off-diagonal element at (p,q)
    void applyRotation(vector<vector<double>>& A, vector<vector<double>>& V, int p, int q, double c, double s) {
        int n = A.size();

        // Update the diagonal elements
        double app = A[p][p], aqq = A[q][q], apq = A[p][q];
        A[p][p] = c*c*app - 2.0*c*s*apq + s*s*aqq;
        A[q][q] = s*s*app + 2.0*c*s*apq + c*c*aqq;
        A[p][q] = A[q][p] = 0.0; // Off-diagonal element becomes zero

        // Update the rest of the matrix
        for (int j = 0; j < n; ++j)
            if (j != p && j != q) {
                double apj = A[p][j], aqj = A[q][j];
                A[p][j] = A[j][p] = c*apj - s*aqj;
                A[q][j] = A[j][q] = s*apj + c*aqj;
            }

        // Update the eigenvector matrix
        for (int i = 0; i < n; ++i) {
            double vip = V[i][p], viq = V[i][q];
            V[i][p] = c*vip - s*viq;
            V[i][q] = s*vip + c*viq;
        }
    }

    // Sort eigenvalues in descending order and rearrange eigenvectors accordingly
    void sortEigenvalues(vector<double>& eigenvalues, vector<vector<double>>& eigenvectors) {
        int n = eigenvalues.size();
        vector<pair<double,int>> pairs;
        for (int i = 0; i < n; ++i) pairs.emplace_back(eigenvalues[i], i);

        // Sort based on eigenvalues
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

    // Main solver function that computes eigenvalues and eigenvectors
    void solve(const vector<vector<double>>& A, vector<double>& eigenvalues, vector<vector<double>>& eigenvectors) {
        int n = A.size();
        vector<vector<double>> A_work = A; // Work on a copy of A
        eigenvectors = vector<vector<double>>(n, vector<double>(n,0.0));
        for (int i=0;i<n;i++) eigenvectors[i][i] = 1.0; // Start with identity matrix

        int iteration = 0; double max_off_diag; int p,q;
        do {
            findMaxOffDiagonal(A_work, p, q, max_off_diag); // Find biggest element to eliminate
            if (abs(max_off_diag) > epsilon) {               // If it's still significant
                double c,s;
                computeRotation(A_work,p,q,c,s);           // Compute rotation
                applyRotation(A_work,eigenvectors,p,q,c,s);// Apply it
            }
            iteration++;
        } while(abs(max_off_diag) > epsilon && iteration < max_iterations); // Keep going until convergence

        // Extract diagonal elements as eigenvalues
        eigenvalues.resize(n);
        for(int i=0;i<n;i++) eigenvalues[i] = A_work[i][i];

        sortEigenvalues(eigenvalues,eigenvectors); // Sort eigenvalues largest first
    }

    // Helper function to print a matrix nicely
    static void printMatrix(const vector<vector<double>>& M, const string& name){
        cout << name << ":\n";
        for(auto &row : M){
            cout << "[ ";
            for(size_t j=0;j<row.size();j++){
                cout << setprecision(6) << row[j];
                if(j < row.size()-1) cout << ", ";
            }
            cout << " ]\n";
        }
        cout << "\n";
    }

    // Helper function to print a vector nicely
    static void printVector(const vector<double>& v, const string& name){
        cout << name << ": [";
        for(size_t i=0;i<v.size();i++){
            cout << setprecision(6) << v[i];
            if(i<v.size()-1) cout << ", ";
        }
        cout << "]\n\n";
    }

    // Print eigenvectors column by column
    static void printEigenvectorsAsColumns(const vector<vector<double>>& V){
        int n = V.size();
        cout << "Eigenvectors (columns):\n";
        for(int col = 0; col < n; ++col){
            cout << "v[" << col << "] = [ ";
            for(int row = 0; row < n; ++row){
                cout << setprecision(6) << V[row][col];
                if(row < n-1) cout << ", ";
            }
            cout << " ]\n";
        }
        cout << "\n";
    }
};

// Validation functions check if results are correct
bool validateEigenSolution(const vector<vector<double>>& A, const vector<double>& eigenvalues, const vector<vector<double>>& V, double tol = 1e-8) {
    // Check if A * V = V * Lambda
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
    // Check if columns of V are orthogonal (dot product ~ 0)
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
    // Check if each eigenvector satisfies A*v = lambda*v
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
    // Check if A is symmetric (required for Jacobi)
    int n = A.size();
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
            if(abs(A[i][j]-A[j][i])>tol) return false;
    return true;
}

// Main program
int main(){
    cout << "=================== Jacobi Eigenvalue Solver ===================\n\n";

    int n;
    cout << "Enter the size of the square matrix: ";
    cin >> n;

    vector<vector<double>> A(n, vector<double>(n,0.0));
    cout << "\nEnter all elements row by row:\n";
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout << "Element [" << i+1 << "," << j+1 << "]: ";
            cin >> A[i][j];
        }
    }

    cout << "\nInput matrix:\n";
    JacobiEigenSolver::printMatrix(A,"Matrix A");

    if(!isSymmetric(A)){
        cout << "Error: The matrix is not symmetric. Jacobi method requires a symmetric matrix.\n";
        return 0;
    } else {
        cout << "Matrix is symmetric. Proceeding...\n\n";
    }

    JacobiEigenSolver solver(1e-10,1000); // Initialize solver with tolerance and max iterations
    vector<double> eigenvalues;
    vector<vector<double>> eigenvectors;
    solver.solve(A,eigenvalues,eigenvectors); // Compute eigenvalues and eigenvectors

    cout << "Computed eigenvalues and eigenvectors:\n";
    JacobiEigenSolver::printVector(eigenvalues,"Eigenvalues");
    JacobiEigenSolver::printEigenvectorsAsColumns(eigenvectors);

    cout << "Validation Results:\n";
    cout << "  Eigen decomposition reconstruction: " << (validateEigenSolution(A,eigenvalues,eigenvectors)?"Passed":"Failed") << "\n";
    cout << "  Eigenvector orthogonality: " << (validateOrthogonality(eigenvectors)?"Passed":"Failed") << "\n";
    cout << "  Eigenvectors satisfy Av = lambda*v: " << (validateEigenvectors(A,eigenvalues,eigenvectors) ? "Passed" : "Failed") << "\n";

    cout << "\n========================= Program Finished =========================\n";

    return 0;
}



















