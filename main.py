import math
from fractions import Fraction
from typing import List, Tuple

Vector = List[int]
RationalVector = List[Fraction]
Matrix = List[Vector]
RationalMatrix = List[RationalVector]


def rational_round(q: Fraction) -> int:
    return int(math.floor(q + Fraction(1, 2)))


class VectorMath:
    @staticmethod
    def dot_product(v1: RationalVector, v2: RationalVector) -> Fraction:
        return sum((x * y for x, y in zip(v1, v2)), Fraction(0))

    @staticmethod
    def norm_sq(v: RationalVector) -> Fraction:
        return VectorMath.dot_product(v, v)


class GramSchmidtManager:
    def __init__(self, basis: Matrix):
        self.n = len(basis)
        self.m = len(basis[0]) if self.n > 0 else 0
        
        self.mu: RationalMatrix = [[Fraction(0) for _ in range(self.n)] for _ in range(self.n)]
        
        self.b_star: RationalMatrix = [[Fraction(0) for _ in range(self.m)] for _ in range(self.n)]
        self.b_star_norm_sq: List[Fraction] = [Fraction(0)] * self.n
        
        self.recompute_from(0, basis)

    def recompute_from(self, start_idx: int, basis: Matrix) -> None:
        for i in range(start_idx, self.n):
            b_i_star = [Fraction(x) for x in basis[i]]
            
            for j in range(i):
                num = sum(Fraction(basis[i][k]) * self.b_star[j][k] for k in range(self.m))
                den = self.b_star_norm_sq[j]
                
                self.mu[i][j] = num / den if den != Fraction(0) else Fraction(0)
                
                for k in range(self.m):
                    b_i_star[k] -= self.mu[i][j] * self.b_star[j][k]
            
            self.b_star[i] = b_i_star
            self.b_star_norm_sq[i] = VectorMath.norm_sq(b_i_star)
            self.mu[i][i] = Fraction(1)


class LLLAlgorithm:
    def __init__(self, basis: Matrix, delta: Fraction = Fraction(3, 4)):
        self.basis = [list(row) for row in basis]
        self.delta = delta
        self.gso = GramSchmidtManager(self.basis)

    def reduce(self) -> Matrix:
        k = 1
        n = len(self.basis)

        while k < n:
            for j in range(k - 1, -1, -1):
                mu_kj = self.gso.mu[k][j]
                if abs(mu_kj) > Fraction(1, 2):
                    q = rational_round(mu_kj)
                    
                    for i in range(len(self.basis[k])):
                        self.basis[k][i] -= q * self.basis[j][i]
                    
                    for i in range(j + 1):
                        self.gso.mu[k][i] -= Fraction(q) * self.gso.mu[j][i]

            lhs = self.gso.b_star_norm_sq[k]
            rhs = (self.delta - self.gso.mu[k][k-1]**2) * self.gso.b_star_norm_sq[k-1]

            if lhs < rhs:
                self.basis[k], self.basis[k-1] = self.basis[k-1], self.basis[k]
                
                self.gso.recompute_from(k - 1, self.basis)
                
                k = max(1, k - 1)
            else:
                k += 1

        return self.basis


def test_lll_edge_cases():
    print("Running defensive mathematical audits...")

    skewed_basis = [
        [1000000, 1],
        [999999, 1]
    ]
    reducer = LLLAlgorithm(skewed_basis)
    reduced = reducer.reduce()
    assert sum(abs(x) for row in reduced for x in row) < 10, "Skewed reduction failed."

    neg_basis = [
        [-5, -7, 2],
        [1, -1, 3],
        [8, -2, -6]
    ]
    reducer_neg = LLLAlgorithm(neg_basis)
    reduced_neg = reducer_neg.reduce()
    gso_check = GramSchmidtManager(reduced_neg)
    delta = Fraction(3, 4)
    for k in range(1, len(reduced_neg)):
        lhs = gso_check.b_star_norm_sq[k]
        rhs = (delta - gso_check.mu[k][k-1]**2) * gso_check.b_star_norm_sq[k-1]
        assert lhs >= rhs, f"Lovász condition failed at index {k}"

    print("All defensive tests passed cleanly.\n")


def solve_subset_sum(weights: List[int], target: int) -> List[int]:
    n = len(weights)
    N = int(math.ceil(math.sqrt(n))) * 1000 
    
    basis = []
    
    for i in range(n):
        row = [0] * n + [weights[i] * N]
        row[i] = 1
        basis.append(row)
        
    target_row = [0] * n + [-target * N]
    basis.append(target_row)

    print(f"Constructed {(n+1)}x{(n+1)} Lattice Basis.")
    
    lll = LLLAlgorithm(basis)
    reduced_basis = lll.reduce()

    for row in reduced_basis:
        if row[-1] == 0:
            subset_indicators = row[:n]
            if all(x in (0, 1) for x in subset_indicators):
                return subset_indicators
            elif all(x in (0, -1) for x in subset_indicators):
                return [-x for x in subset_indicators]

    return []


if __name__ == "__main__":
    test_lll_edge_cases()

    public_weights = [2, 3, 7, 14, 30, 57, 120, 251]
    target_sum = 137 
    
    print(f"Attempting to break subset-sum encryption...")
    print(f"Public Weights: {public_weights}")
    print(f"Target Sum: {target_sum}")

    solution_vector = solve_subset_sum(public_weights, target_sum)

    if solution_vector:
        print("\n[+] SUCCESS: Short vector found mapping to exact sum.")
        print(f"Solution Indicator Vector: {solution_vector}")
        elements = [public_weights[i] for i, val in enumerate(solution_vector) if val == 1]
        print(f"Extracted Elements: {elements} (Sum = {sum(elements)})")
    else:
        print("\n[-] FAILED: LLL did not converge on the target vector. Density may be too high.")
