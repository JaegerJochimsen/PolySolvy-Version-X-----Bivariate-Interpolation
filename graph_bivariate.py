def Graph_Bivariate(vmatrix: np.matrix, x: np.array, y:np.array, original_z: np.array):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    calculated_z = []
    
    # generate some sample x,y vals to plot, this helps to fill out the curve
    sample_x = np.arange(x.min() - 1,x.max() + 1,.1)
    sample_y = np.arange(y.min() - 1, y.max() + 1, .1)

    # make sure to include our original x's and y's as well
    for i in range(len(x)):
        sample_x = np.append(sample_x, x[i])
        sample_y = np.append(sample_y, y[i])

    # Create the array of z values based off of constants in vmatrix and original x and y vals
    for index in range(len(sample_x)):
        calculated_z.append(np.double(Bivariate(vmatrix, sample_x[index],sample_y[index])))
        
    # convert to correct data type
    calculated_z = np.array(calculated_z)

    # a plot of our parametric curve overlayed over our original data points
    ax.plot(sample_x,sample_y,calculated_z)
    ax.scatter(x,y,original_z, c='red')

# ------------------- #
#   A Neat Example    #
# ------------------- #

# coeffs calculated using gaussian elim
v_matrix = np.matrix([[8/3,-7/3],[0,2/3]])   

# our original x's, y's, and z's
xs = np.array([1,3,4,4]) 
ys = np.array([1,2,1,4])
zs = np.array([1,2,3,4]) 
Graph_Bivariate(v_matrix, xs, ys, zs)
