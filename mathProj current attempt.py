import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
%matplotlib notebook

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

def Bivariate(v_matrix: np.array, x: np.double, y: np.double):
    """
    A function that iterates through a matrix of polynomial coefficients and multiplies them
    with the corresponding (x^n)(y^n) term. In this case, the index of the matrix (ith row, jth col)
    corresponds to the powers on x and y.
    param: v_matrix: a matrix of type np.matrix which contains the coefficients for the interpolating polynomial
    param: x: an x value of type np.double to be used in determining the z value at a point
    param: y: a y value of type np.double to be used in determining the z value at a point
    return: None if v_matrix is empty, otherwise return the polynomial evalueated at (x) and (y), return type should
            be np.double. This is evaluated as:
            v_matrix[0][0]*(x^0)*(y^0) + v_matrix[0][1]*(x^0)*(y^1) + ... + v_matrix[n-1][n-1]*(x^(n-1))*(y^(n-1))
            (where n is the number of rows [and columns since it is an NxN matrix] in v_matrix)
    # [[1, 2],        # 1 is element (0,0), 2 is element (0,1), 4 is element (1,1), etc.
    # [3, 4]]
    Ex/
    >>> v_matrix = [[1,2], [3,4]]
    >>> Bivariate(v_matrix, 1, 1)
    10
    >>> Bivariate(v_matrix, 3,4)
    66
    """
    length = len(v_matrix)  
    if length == 0:
        return None

    z = np.double(0)

    c = 0
    r = 0
    for i in range(length):
        if ((i + 1) % math.sqrt(length)) == 0:
            z += v_matrix[i]*(x**r)*(y**c)
            r += 1
            c = 0
        else:
            z += v_matrix[i]*(x**r)*(y**c)
            c += 1
    return z

def Graph_Bivariate(vmatrix: np.matrix, orig_x: np.array, orig_y:np.array, orig_z: np.array):
    """
    A function which plots an original set of (x,y,z) points in a 3D space and overlays that graph with a graph
    of the function output by Bivariate() (this is the interpolating function for the original (x,y,z) points).
    param:  vmatrix: an np.matrix of coefficients of the terms of the interpolating polynomial. These coefficients
            are of the form:
            [[c,y,y^2,y^3,...........................y^(n//4)],         // powers of y increase horizontally, powers of x vertically.
             [x, xy, xy^2,.........................,xy^(n//4)],
             [x^2,(x^2)y,(x^2)y^2,...........................]
             [x^3,(x^3)y,....................................]
             ................................................]
             ................................................]
             [(x^(n//4)),(x^(n//4))y,..., (x^(n//4))(y^(n//4)]]
             These coefficients are used by Bivariate() to compute the corresponding z-value for a given (x,y) pair as
             specified by the interpolating polynomial.
    param: x: an np.array of the original x-coordinates that were interpolated.
    param: y: an np.array of the original y-coordinates that were interpolated.
    param: z: an np.array of the original z-cooordinates that were interpolated.
    return: outputs a graph plotted using Axes3D from mpl_toolkits.mplot3d
    libraries/files used:
    - bivariate.py
    - numpy
    - matplotlib.pyplot
    - mpl_toolkits.mplot3d
    """
    # basic error catching
    # input arrays must be of equal size so that they represent viable points
    if (len(orig_x) != len(orig_y)) or (len(orig_x) != len(orig_z)):
        raise Exception("Input arrays are of unequal length")
    
    # matrix of coeffs should be of size sqrt(n) where n is the number of points/coordinates being used
   # if len(vmatrix) != math.sqrt(len(orig_x)):
    #    raise Exception("Coefficient matrix (vmatrix) is of incorrect size")
    
    # set up the initial graphing space and initialize as a 3d space
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
   
    calculated_z = []
    
    # generate some sample x,y vals to plot, this helps to fill out the curve
    # make sure that the range of value for the x and y's are the same length
    # take the larger of the two ranges for this purpose
    if (orig_x.max() - orig_x.min()) > (orig_y.max() - orig_y.min()):
        # if x is the larger range
        sample_x = np.arange(orig_x.min(),orig_x.max()+1,.1)
        sample_y = np.arange(orig_x.min(), orig_x.max() + 1, .1)
    else:
        sample_x = np.arange(orig_y.min(),orig_y.max()+1,.1)
        sample_y = np.arange(orig_y.min(), orig_y.max() + 1, .1)
        
    # ensure that we include the original points
    sample_x = np.append(sample_x, orig_x)
    sample_y = np.append(sample_y, orig_y)
    
    # Create the array of z values based off of constants in vmatrix and original x and y vals
    for index in range(len(sample_x)):
        calculated_z.append(np.double(Bivariate(vmatrix, sample_x[index],sample_y[index])))
        
    # convert to correct data type
    calculated_z = np.array(calculated_z, dtype=np.double)
    
    # if an odd num of data points, add one to make even
    if len(sample_x) % 2 != 0:
        newX = random.randint(orig_x.min(), orig_x.max())
        sample_x = np.append(sample_x, newX)
        
        newY = random.randint(orig_y.min(), orig_y.max())
        sample_y = np.append(sample_y, newY)
        
        newZ = Bivariate(vmatrix, newX, newY)
        calculated_z = np.append(calculated_z, newZ)
    
    # reshape the data for the purpose of surface plotting
    reshapedx = sample_x.reshape(2, len(sample_x)//2)
    reshapedy = sample_y.reshape(2, len(sample_y)//2)
    length = len(calculated_z) // 2
    reshapedz = calculated_z.reshape(2, length)
    
    # surface plot
    #ax.plot_surface(reshapedx, reshapedy, reshapedz)
    ax.plot_trisurf(sample_x, sample_y, calculated_z)
    
    #ax.plot(sample_x, sample_y, calculated_z)
    
    # plot of original data points
    ax.scatter(orig_x, orig_y, orig_z, c='red')
    
    # rotate the axes and update for an interactive graph
    for angle in range(0, 360):
        ax.view_init(30, angle)
        plt.draw()
        plt.pause(.001)


# ------------------------ #
#   A Few Neat Examples    #
# ------------------------ #

# coeffs calculated using gaussian elim
# here our coeffs are in the order: [[c, y], [x,xy]] (see docstring for example)
# our original x's, y's, and z's
xs = np.array([1,3,4,5], dtype=np.double) 
ys = np.array([1,2,6,4],dtype=np.double)
zs = np.array([1,2,3,4],dtype=np.double)
our_v = Vandermonde(xs, ys, zs)
Graph_Bivariate(our_v, xs, ys,zs)
   
#Ex/2
xs2 = np.array([2,3,4,5], dtype=np.double)
ys2 = np.array([3,4,7,8],dtype=np.double)
zs2 = np.array([4,7,8,12],dtype=np.double)
our_v2 = Vandermonde(xs2,ys2, zs2)
Graph_Bivariate(our_v2, xs2, ys2,zs2)

#Ex/3 9-points (getting spicy)
xs3 = np.array([1,1,1,2,2,2,3,3,3], dtype=np.double)
ys3 = np.array([1,2,3,1,2,3,1,2,3], dtype=np.double)
zs3 = np.array([3.2,4.4,6.5,2.5,4.7,5.8,5.1,3.6,2.9], dtype=np.double)
our_v3 = Vandermonde(xs3,ys3,zs3)
Graph_Bivariate(our_v3, xs3, ys3, zs3)
