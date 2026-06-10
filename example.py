from lll import LLLAlgorithm
from fractions import Fraction

print("Initializing LLL Reduction...")

# Define a basis as a list of integer arrays
basis = [
    [1, -1, 3],
    [100005, 1, 2],
    [8, -2, -6]
]

# Initialize the reducer (delta defaults to 3/4)
reducer = LLLAlgorithm(basis, delta=Fraction(3, 4))
reduced_basis = reducer.reduce()

print("\nReduced Basis:")
for row in reduced_basis:
    print(row)
