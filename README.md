# Jacobi Eigenvalue Solver

A C++ implementation of the **Jacobi method** to compute **eigenvalues** and **eigenvectors** of a **symmetric matrix**.

---

## Features

- Computes eigenvalues and eigenvectors using the Jacobi rotation method.
- Automatically sorts eigenvalues in descending order and rearranges corresponding eigenvectors.
- Validates:
  - Reconstruction \(A = V \Lambda V^T\)
  - Eigenvector orthogonality
  - Eigenvectors satisfy \(Av = \lambda v\)
- User-friendly input/output for matrices and vectors.

---

## Usage

1. Clone the repository:

```bash
git clone https://github.com/yourusername/jacobi-eigen-solver.git
cd jacobi-eigen-solver
