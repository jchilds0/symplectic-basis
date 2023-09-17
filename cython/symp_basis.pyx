"""

Cython file

"""

cdef extern from "triangulation.h":
    ctypedef struct c_Triangulation "Triangulation":
        pass

cdef extern from "SnapPea.h":
    void peripheral_curves(c_Triangulation *manifold)

cdef extern from "unix_file_io.h":
    c_Triangulation *read_triangulation_from_string(char *string)

cdef extern from "symplectic_basis.h":
    int **get_symplectic_basis(c_Triangulation *manifold, int *num_rows, int *num_cols, int log)
    void free_symplectic_basis(int **eqns, int num_rows)

def symplectic_basis(manifold, verify=False, debug=False):
    """
    Extend the Neumann-Zagier Matrix to one which is symplectic (up to factors of 2)
    using oscillating curves. Verify parameter explicitly tests if the resulting matrix is symplectic.

    >>> M = Manifold("4_1")
    >>> symplectic_basis(M)
    [-1  0 -1 -1]
    [ 2  0 -2  0]
    [-2 -1 -2 -1]
    [ 2  0  0  0]

    <https://arxiv.org/abs/2208.06969>
    """
    def to_byte_str(s):
        return s.encode('utf-8') if type(s) != bytes else s

    def is_symplectic(M):
        """
        Test if the matrix M is symplectic
        :param M: square matrix
        :return: true or false
        """
        n = len(M)

        for i in range(n):
            for j in range(i, n):
                omega = abs(symplectic_form(M[i], M[j]))

                if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                    if omega != 2:
                        return False
                elif omega:
                    return False

        return True

    def symplectic_form(u, v):
        return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i] for i in range(len(u) // 2)])

    cdef int **c_eqns;
    cdef int **g_eqns;
    cdef int num_rows, num_cols, dual_rows;
    cdef int* eqn;
    cdef c_Triangulation *c_triangulation = NULL;

    string = manifold._to_string()
    b_string = to_byte_str(string)
    c_triangulation = read_triangulation_from_string(b_string)

    if c_triangulation is NULL:
        raise ValueError('The Triangulation is empty.')

    peripheral_curves(c_triangulation)

    # Cusp Equations
    gluing_eqns = manifold.gluing_equations()
    h, w = gluing_eqns.dimensions()

    eqns = []
    for i in range(h - 2 * manifold.num_cusps(), h):
        eqns.append([gluing_eqns[i][j] for j in range(w)])

    # Dual Curve Equations
    g_eqns = get_symplectic_basis(c_triangulation, &dual_rows, &num_cols, int(debug))

    for i in range(dual_rows):
        eqns.append([g_eqns[i][j] for j in range(num_cols)])

    free_symplectic_basis(g_eqns, dual_rows)

    # Convert to Neumann Zagier Matrix
    rows = len(eqns)
    retval = [[eqns[i][3 * (j // 2) + j % 2] - eqns[i][3 * (j // 2) + 2] for j in range(rows)] for i in range(rows)]

    if verify:
        if is_symplectic(retval):
            print("Result is symplectic (up to factors of 2)")
        else:
            print("Warning: Result is not symplectic")

    return retval