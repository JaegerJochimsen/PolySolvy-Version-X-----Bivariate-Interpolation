import bivariate 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
    # set up the initial graphing space and initialize as a 3d space
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    calculated_z = []
    
    # generate some sample x,y vals to plot, this helps to fill out the curve
    # here the values range from the min to the max + 1 elements in each original array 
    sample_x = np.arange(orig_x.min(),orig_x.max()+1,.1)
    sample_y = np.arange(orig_y.min(), orig_y.max() + 1, .1)

    # ensure that we include the original points
    sample_x = np.append(sample_x, orig_x)
    sample_y = np.append(sample_y, orig_y)

    # Create the array of z values based off of constants in vmatrix and original x and y vals
    for index in range(len(sample_x)):
        calculated_z.append(np.double(bivariate.Bivariate(vmatrix, sample_x[index],sample_y[index])))

    # convert to correct data type
    calculated_z = np.array(calculated_z)

    # a plot of our parametric curve overlayed over our original data points
    ax.plot(sample_x,sample_y,calculated_z)
    ax.scatter(orig_x, orig_y, orig_z, c='red')

# ------------------- #
#   A Neat Example    #
# ------------------- #
# coeffs calculated using gaussian elim
v_matrix = np.matrix([[8/3,-7/3],[0,2/3]])   # here our coeffs are in the order: [[c, y], [x,xy]] (see docstring for example)

# our original x's, y's, and z's
xs = np.array([1,3,4,4]) 
ys = np.array([1,2,1,4])
zs = np.array([1,2,3,4]) 
Graph_Bivariate(v_matrix, xs, ys, zs)
