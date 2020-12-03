import math
import bivariate 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib notebook

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
    if len(vmatrix) != math.sqrt(len(orig_x)):
        raise Exception("Coefficient matrix (vmatrix) is of incorrect size")
    
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
        calculated_z.append(np.double(bivariate.Bivariate(vmatrix, sample_x[index],sample_y[index])))

    # convert to correct data type
    calculated_z = np.array(calculated_z)

    # reshape the data for the purpose of surface plotting
    reshapedx = sample_x.reshape(2, len(sample_x)//2)
    reshapedy = sample_y.reshape(2, len(sample_y)//2)
    reshapedz = calculated_z.reshape(2, len(calculated_z)//2)
    
    # surface plot
    ax.plot_surface(reshapedx, reshapedy, reshapedz)
    
    # a plot of original data points
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
v_matrix = np.matrix([[(9/8), (-11/16)],[(3/8),(3/16)]])   # here our coeffs are in the order: [[c, y], [x,xy]] (see docstring for example)

# our original x's, y's, and z's
xs = np.array([1,3,4,5]) 
ys = np.array([1,2,6,4])
zs = np.array([1,2,3,4]) 
Graph_Bivariate(v_matrix, xs, ys, zs)


# Here is a good example, uses evenly spaced points without repeats
xs2 = np.array([1,2,3,4])
ys2 = np.array([5,6,7])
zs2 = np.array([9,10,11,12])
v_matrix = np.matrix([[-4,3],[-2,0]])
Graph_Bivariate(v_matrix,xs2,ys2,zs2)

