# RationalLLL

A rigorous, standalone Python implementation of the Lenstra–Lenstra–Lovász (LLL) lattice reduction algorithm. 

This implementation is designed specifically to avoid the catastrophic cancellation and precision loss inherent in standard IEEE 754 floating-point arithmetic. By enforcing exact rational arithmetic via Python's `fractions` module and strictly decoupling the Gram-Schmidt Orthogonalization (GSO) state management, **RationalLLL** guarantees absolute computational accuracy even in high-dimensional or highly skewed lattices.

## Features

* **Exact Rational Arithmetic:** All GSO projections and lengths are maintained in exact fractions, preventing state drift.
* **Decoupled GSO Manager:** The orthogonalization state is isolated from the main integer-basis mutations, dynamically and deterministically recomputing local states during Lovász condition failures.
* **Zero Dependencies:** Built entirely on the Python Standard Library. 
* **Defensive Edge-Case Handling:** Audited against highly skewed matrices and linear dependence risks.
* **Practical Cryptographic Application:** Includes a built-in Lagarias-Odlyzko subset-sum (knapsack) solver.

## Mathematical Methodology

Lattice algorithms are highly susceptible to precision errors. In the LLL algorithm, standard floating-point arithmetic quickly suffers from round-off accumulation during the projection phase. 

This architecture tracks two distinct mathematical states:
1. **The Lattice Basis (B):** Resides strictly in the integer domain. All size reductions and swaps are applied directly to these integer coordinates.
2. **The GSO State (B and μ):** Evaluated strictly in the rational domain. When b_k = b_k - q * b_j, we execute an O(1) update to the coefficients. When the Lovász condition fails, the local orthogonality is fractured, and the GSO state is deterministically recomputed from the fractured index onward.

### The Lovász Condition
The algorithm tests whether swapping two adjacent basis vectors improves the quality of the basis using the adjustable parameter delta (default is 3/4):

||b_k*||^2  >=  (delta - μ_{k, k-1}^2) * ||b_{k-1}*||^2

## Installation & Usage

Clone the repository and run the standalone script. No external dependencies (like NumPy or SciPy) are required.

```bash
git clone https://github.com/Raj123-0/RationalLLL.git
cd RationalLLL
python lll.py
