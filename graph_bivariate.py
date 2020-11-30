"""
Project: PolySolvy Version X
Math 351 - Numerical Analysis
Under the instruction of Professor Eric Merchant

Team: Nathan, Jaeger, Ronnie

File Description:
graph_bivariate.py uses the functionality of bivariate.py to generate a 3D surface plot of an interpolating function
and the points it interpolates (of the (x,y,z)). Graph_Bivariate() generates a graph which displays the surface and
runCounter.py
"""
import bivariate 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def Graph_Bivariate(vmatrix: np.matrix, x: np.array, y:np.array, original_z: np.array):
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
v_matrix = np.matrix([[8/3,-7/3],[0,2/3]])    
x = np.array([1,3,4,4]) 
y = np.array([1,2,1,4])   
zs = np.array([1,2,3,4]) 
Graph_Bivariate(v_matrix, x, y, zs)

