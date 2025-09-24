# Jacobi Eigenvalue Solver

A robust C++ implementation of the classical Jacobi eigenvalue algorithm for computing eigenvalues and eigenvectors of symmetric matrices.

## 🔍 Overview

This program implements the classical Jacobi method to find all eigenvalues and eigenvectors of a real symmetric matrix. The algorithm works by iteratively applying rotation matrices to eliminate off-diagonal elements, eventually diagonalizing the matrix while accumulating the eigenvectors.

## ✨ Features

- **Complete Eigenvalue Decomposition**: Computes all eigenvalues and eigenvectors
- **Automatic Sorting**: Eigenvalues sorted in descending order with corresponding eigenvectors
- **Built-in Validation**: Comprehensive validation checks for accuracy
- **Interactive Interface**: User-friendly input and formatted output
- **Configurable Parameters**: Customizable tolerance and iteration limits
- **No Dependencies**: Pure C++ with STL only

## 🧮 Algorithm Details

The Jacobi method is an iterative algorithm that:

1. **Identification**: Finds the largest off-diagonal element
2. **Rotation Calculation**: Computes the optimal rotation angle to eliminate it
3. **Matrix Transformation**: Applies the rotation to the matrix and accumulates transformations
4. **Iteration**: Repeats until all off-diagonal elements are negligible
5. **Extraction**: Returns diagonal elements as eigenvalues and accumulated rotations as eigenvectors

## 📋 Requirements

- C++ compiler with C++11 support or higher
- Standard Template Library (STL)
- No external dependencies required

## 🔨 Compilation

```bash
g++ -o jacobi_eigen jacobi_eigen.cpp -std=c++11
```

For optimized build:
```bash
g++ -O3 -o jacobi_eigen jacobi_eigen.cpp -std=c++11
```

## 🚀 Usage

1. **Run the program**:
   ```bash
   ./jacobi_eigen
   ```

2. **Enter matrix size** when prompted

3. **Input matrix elements** row by row (symmetric matrices only)

4. **View results**: The program displays:
   - Original input matrix
   - Computed eigenvalues (sorted in descending order)
   - Corresponding eigenvectors
   - Validation results

## 💡 Example

**Input Matrix:**
```
[2, 1]
[1, 2]
```

**Output:**
```
Eigenvalues: [3.000000, 1.000000]
Eigenvectors: 
v[0] = [0.707107, 0.707107]
v[1] = [-0.707107, 0.707107]
```

## ✅ Validation

The program performs three comprehensive validation checks:

1. **Eigendecomposition Verification**: Confirms `A × V = V × Λ`
2. **Orthogonality Check**: Ensures eigenvectors are orthogonal
3. **Individual Eigenpair Validation**: Verifies `A·v = λ·v` for each pair

## 🔧 Key Functions

| Function | Description |
|----------|-------------|
| `JacobiEigenSolver` | Main solver class with configurable parameters |
| `solve()` | Computes eigenvalues and eigenvectors |
| `findMaxOffDiagonal()` | Locates largest off-diagonal element |
| `computeRotation()` | Calculates Jacobi rotation parameters |
| `applyRotation()` | Applies rotation to matrix and eigenvectors |

## ⚙️ Customization

Modify solver parameters in `main()`:

```cpp
// Custom tolerance and maximum iterations
JacobiEigenSolver solver(1e-10, 1000);
```

Default parameters:
- **Tolerance**: `1e-9`
- **Max Iterations**: `500`

## ⚠️ Limitations

- **Symmetric Matrices Only**: Algorithm requires symmetric input matrices
- **Computational Complexity**: O(n³) - suitable for small to medium matrices (n < 1000)
- **Linear Convergence**: Reliable but slower than some modern methods for large matrices

## 🎯 Applications

- **Principal Component Analysis (PCA)**
- **Vibration Analysis** in mechanical systems
- **Quantum Mechanics** calculations
- **Image Processing** and computer vision
- **Machine Learning** feature extraction
- **Structural Analysis** and modal analysis
- Any application requiring **spectral decomposition** of symmetric matrices

## 📊 Performance Notes

- Numerically stable for well-conditioned matrices
- Eigenvectors are automatically normalized and orthogonal
- Results include built-in validation
- Ideal for educational purposes and scientific computing

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📄 License

This project is open source. Please check the repository for license details.

## 📚 References

- Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.)
- Press, W. H., et al. (2007). *Numerical Recipes: The Art of Scientific Computing*

---

**Perfect for educational purposes, research, and small-scale scientific computing applications!**


