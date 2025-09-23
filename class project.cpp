


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