import numpy as np
import vmond as v
import math

superscript = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

def P_vals(i, j):
    # For Visual purposes.
    # Output is string version of powers.
    y_pow = "y" + f'{j}'.translate(superscript)
    x_pow = "x" + f'{i}'.translate((superscript))
    if i == 0 and j == 0:
        return 1
    elif i == 0:
        return y_pow
    elif j == 0:
        return x_pow
    return x_pow + y_pow

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
    V = [([1] * N) for i in range(N)]
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
            #if (j == i) and (i == 0 and j == 0):
                # Keep the first item 1 (C).
             #   continue
            #else:
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

    #print('\n'.join([" ".join(['{:<5}'.format(item) for item in row]) for row in P]))

    # Takes care of pattern for building V matrix
    for i in range(Ns):
        for j in range(Ns):
            powArray.append(P[i][j])
    #print(powArray)

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

    print('\n'.join([''.join(['{:7}'.format(item) for item in row]) for row in V]))

    C = v.Gauss(V, Z)
    return C


a = np.array([1, 2, 3, 4], dtype=np.double)
b = np.array([1, 2, 3, 4], dtype=np.double)
c = np.array([1, 2, 3, 4], dtype=np.double)
Vandermonde(a, b, c)
