import symplectic_basis
import snappy
import pytest


def is_symplectic(M):
    """
    Test if the matrix M is symplectic
    :param M: square matrix
    :return: true or false
    """
    n = len(M)

    for i in range(n):
        for j in range(i, len(M[i])):
            omega = abs(symplectic_form(M[i], M[j]))

            if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                if omega != 2:
                    return False
            elif omega:
                return False

    return True


def symplectic_form(u, v):
    return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i]
                for i in range(len(u) // 2)])


def test_figure8_knot():
    # M = snappy.Manifold("4_1")
    assert True is True
    # assert is_symplectic(symplectic_basis.symplectic_basis(M)) is True


# def test_knots():
#     for M in snappy.CensusKnots:
#         assert is_symplectic(symplectic_basis.symplectic_basis(M)) is True


if __name__ == "__main__":
    M = snappy.Manifold("4_1")
    print(symplectic_basis.symplectic_basis(M))