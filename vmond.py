import numpy as np
import vmond as v
import math

superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

class MatrixSolveError(Exception):
    """ Any error intentionally raised by Gauss or Vandermonde. """
    pass


def Gauss(A: np.array, b: np.array) -> np.array:
    """
    Gaussian elimination procedure (aka Row Reduction) for
    solving linear systems of the form Ax = b.
    Scaled Partial Pivoting (SPP) is used in order to avoid
    numerical errors associated with the Naive Algorithm
    (read C&K sections 7.1-7.2 for more details).
    Assumes array A is square (N x N) and
    solution vector b is N-dimensional.
    Returns values of vector x.
    WARNING: To avoid numerical divide by zero / overflow errors,
    make sure you specify the datatype (dtype) to be np.double!
    References:
        + Phase 1 Row Reduction Alg: C&K 6th ed. p. 267
        + Phase 2 Back Substitution Alg: C&K 6th ed. p. 269
    # - - - - - #
    #  TESTING  #
    # - - - - - #
    Basic Case: Let's get a little more adventurous.
    >> A = np.array([[3, -13, 9, 3], [-6, 4, 1, -18], [6, -2, 2, 4], [12, -8, 6, 10]], dtype=np.double)
    >> b = np.array([-19, -34, 16, 26], dtype=np.double)
    >> Gauss(A, b)
    array([ 3.,  1., -2.,  1.])
    """
    N = len(A)
    b = np.copy(b)

    # Build coefficient vector
    x = np.array([0 for i in range(N)], dtype=np.double)

    # Index list prioritizes which row
    # should be scaled first.
    #
    # This is necessary because Naive
    # Gaussian elimination is prone to
    # yield errors if the natural order
    # is followed (C&K 6th ed., p. 261-262)
    Index = [i for i in range(N)]

    # Scale list keeps track of what scale
    # factors are needed to normalize the
    # coefficients for row i in array A.
    Scale = [0 for i in range(N)]

    # - - - - - - - - - - - - - - - - - - - - - - - #
    #   PHASE 1 - Row Reduction w/ Scaled Pivots    #
    # - - - - - - - - - - - - - - - - - - - - - - - #

    # init step - compute values of scale factors for each row
    for i in range(N):
        SMAX = 0
        for j in range(N):
            SMAX = max(SMAX, abs(A[i][j]))
        Scale[i] = SMAX

    # Row Reduction - go through rows, prioritize which one
    # to use as a "pivot", and eliminate entries in corresponding
    # column in all other rows.
    for k in range(N - 1):
        RMAX = 0
        for i in range(k, N):
            l_i = Index[i]
            s_i = Scale[l_i]
            r = abs(A[l_i][k] / s_i)

            # If scaled ratio is higher than previous encounter,
            # then this row should take higher priority.
            if r > RMAX:
                RMAX = r
                j = i

        # Swap terms to reorganize which rows take higher priority.
        SWP = Index[j]
        Index[j] = Index[k]
        Index[k] = SWP

        l_k = Index[k]
        # Eliminate entries in all other rows != row k.
        for i in range(k + 1, N):
            l_i = Index[i]

            XMULT = (A[l_i][ k])/ (A[l_k][ k])
            A[l_i][ k] = XMULT

            for j in range(k + 1, N):
                A[l_i][ j] = A[l_i][ j] - XMULT * A[l_k][ j]

    # - - - - - - - - - - - - - - - - #
    #   PHASE 2 - Back Substitution   #
    # - - - - - - - - - - - - - - - - #

    # Using the current values of Index and the stored coefficients in A,
    # we will now alter the values of solution vector b
    # to match the manipulations we have already done
    # to the coefficient array A.
    for k in range(N - 1):
        l_k = Index[k]
        for i in range(k + 1, N):
            l_i = Index[i]
            b[l_i] = b[l_i] - A[l_i][ k] * b[l_k]

    l_n = Index[N - 1]
    x[N - 1] = b[l_n] / A[l_n][ N - 1]

    # We are now well equipped to solve for the values
    # of our x-vector. This is the final step of the algorithm.
    for i in range(N - 2, -1, -1):
        l_i = Index[i]
        SUM = b[l_i]
        for j in range(i + 1, N):
            SUM = SUM - A[l_i][ j] * x[j]

        x[i] = SUM / A[l_i][ i]

    return x


def powers(i, j):
    # Helper function to match each x & y pair with
    # their corresponding powers.
    y_pow = ("y", j)
    x_pow = ("x", i)
    if i == 0 and j == 0:
        return (1, 0)
    elif i == 0:
        return y_pow
    elif j == 0:
        return x_pow
    return x_pow + y_pow


def Vandermonde(X, Y, Z):
    C = np.array([])
    N = len(X)
    Ns = int(math.sqrt(N))
    V = np.zeros(shape=(N,N), dtype=np.double)
    #V = [([1] * N) for i in range(N)]
    P = [([1] * Ns) for i in range(Ns)]
    powArray = []

    # Check for valid input. Exit if invalid
    try:
        assert len(X) == len(Y) == len(Z)
    except AssertionError:
        raise NotImplementedError("Tried to call Vandermonde with unbalanced X, Y, and Z arrays!")

    # ******************** Build P: (call by need) ********************
    for i in range(Ns):
        for j in range(Ns):
            P[i][j] = powers  # change back to P_vals if experiment failed?

    # Verify the output of P matrix reflects scratch work.
    # Potential Issue: Is that pattern ok?
    # Run by J & Nate
    for i in range(Ns):
        for j in range(Ns):
            if P[i][j] == 1:
                continue
            else:
                P[i][j] = P[i][j](i,j)

    # Takes care of pattern for building V matrix
    for i in range(Ns):
        for j in range(Ns):
            powArray.append(P[i][j])

    # ******************** Build V ********************
    #This doesn't work with P_vals functions. Change it to powers!
    for i in range(N):
        for j in range(N):
            if powArray[j][0] == 1:
                V[i][j] = (Y[j]**powArray[j][1]) * (X[j]**powArray[j][1])
            elif len(powArray[j]) != 2:
                V[i][j] = (Y[i]**powArray[j][3]) * (X[i]**powArray[j][1])
            elif powArray[j][0] == "y":
                V[i][j] = (Y[i]**powArray[j][1])
            elif powArray[j][0] == "x":
                V[i][j] = (X[i] ** powArray[j][1])

    C = Gauss(V, Z)
    return C
