# PNMS-Real-Schur-Decomposition

## Description
This project implements the computation of the **real Schur decomposition** of a real square matrix using the **implicit QR algorithm**, developed in MATLAB/Octave.  
It was completed as part of a university course in numerical linear algebra.

Given a real matrix $A \in \mathbb{R}^{n \times n}$, the algorithm computes an orthogonal matrix $Q$ such that

$$Q^T A Q = T$$

where $T$ is in **real Schur form**, i.e., a quasi-upper triangular matrix with $1 \times 1$ blocks corresponding to real eigenvalues and $2 \times 2$ blocks corresponding to complex conjugate eigenvalue pairs.

## Mathematical Background
The implementation is based on the **implicitly shifted QR iteration**, which applies successive orthogonal similarity transformations to drive subdiagonal entries to zero while preserving eigenvalues.  

The implementation reduces $A$ to Hessenberg form $H$ in order to improve on the number of operations required to perform a QR iteration, then uses implicit shift strategies to improve convergence while maintaining numerical stability.

A relative conditio is used to detect deflation which serves to split the matrix into smaller block, further improving convergence.

$$|H(i+1, i)| \leq tol*(|H(i, i)| + |H(i+1, i+1)|), \quad i = 1, ..., n$$

## Implementation Details
- Language: MATLAB
- Input: Real square matrices
- Output: Real Schur form and associated orthogonal factors
- Emphasis on algorithmic clarity and numerical correctness

## References
Golub, G. H., & Van Loan, C. F.  
*Matrix Computations*, 4th Edition, Johns Hopkins University Press.
