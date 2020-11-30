"""
Project: PolySolvy Version X
Math 351 - Numerical Analysis
Under the instruction of Professor Eric Merchant

Team: Nathan, Jaeger, Ronnie

File Description:

bivariate.py uses the functionality of vmond.py to generate a matrix of coefficients for
a bivariate interpolating polynomial. It implements a naive function Bivariate() which combines the
coefficients contained in the matrix with the appropriate (x^n)(y^n) terms to output a z-value
for a given x and y value. Specific procedure is included in function docstring.
"""

#import vmond
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------- #
# Get a list of points to interpolate #
# ----------------------------------- #
# FIXME: an input loop to accept points? or just define as global var?
# X = []
# Y = []
# Z = []
# VMATRIX = vmond.Vandermonde(X, Y, Z)      # this is a numpy matrix

def Bivariate(v_matrix: np.matrix, x: np.double, y: np.double):
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
    if len(v_matrix) == 0:
        return None

    z = 0
    for r in range(len(v_matrix)):
        for c in range(len(v_matrix)):
            z += v_matrix[r,c]*(x**r)*(y**c)
    return z


def Graph_Bivariate(vmatrix: np.matrix, x: np.array, y:np.array):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    z = []

    # Create the array of z values based off of constants in vmatrix and x and y vals
    for index in range(len(x)):
        z.append(Bivariate(vmatrix, np.double(x[index]),np.double(y[index])))
    z = np.array(z)

    #graph surface via trisurf in blue and overlay points in red
    ax.plot_trisurf(x,y,z,linewidth=0, triangles='triangles')
    ax.scatter(x,y,z, c='red')


# ------------------- #
#   A Neat Example    #
# ------------------- #
vm = np.matrix([[0.97500,-5.27500,5.95000],[ -3.92500, 19.82500, -21.55000], [3.40000, -14.70000, 18.50000]])    
xs = np.array([1,1,1,2,2,2,3,3,3])
ys = np.array([1,2,3,1,2,3,1,2,3])

Graph_Bivariate(vm, xs, ys)

