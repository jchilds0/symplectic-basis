import symplectic_basis
import snappy

if __name__ == "__main__":
    M = snappy.Manifold("4_1")
    print(symplectic_basis.symplectic_basis(M))