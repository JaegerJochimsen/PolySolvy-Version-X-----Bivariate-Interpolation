
## Underlying Background and Theory: ##

The problem of finding an interpolating polynomial for a dataset in multiple variables is an area of extensive numerical research. In this Python project, we explore how this may be done practically for bivariate datasets representing a mystery function f(x,y) -> z.

Recall that interpolation is the problem of constructing a polynomial function p belonging to a finite dimensional linear space from a given set of data (Gasca et al., 2000). p is said to interpolate f if for every sample point x_i within the domain of f, p(x_i) = f(x_i) with 0 <= i <= n.

As we have seen in our numerical analysis class this term, this univariate case is fairly straightforward to solve. From a theoretical standpoint, we saw that a unique interpolation polynomial can be determined from the Lagrange form ...

<img src="https://latex.codecogs.com/gif.latex?p(x)=\sum_{i=0}^{n}{y_il_i(x)}"/>
<img src="https://latex.codecogs.com/gif.latex?l_i(x)=\prod_{j=0,j\neq{i}}^{n}{\frac{x-x_j}{x_i-x_j}}"/>

Where i = 0, .., n for n + 1 interpolation nodes, and l_i represents the cardinal polynomial for the ith node. In this respect, the y_i terms act as scaling factors and the l_i acts as a continuous analog to the binary Kronecker delta (Cheney & Kincaid, 2008, p. 126-127).

Other approaches to the interpolation problem include the Hermite and Birkhoff methods - both of which utilitize not only the values of the interpolation nodes at given points but also the specified values of their derivatives (Wikipedia, 2020). The Newton approach we learned in class this term approaches the problem from a recursive standpoint using divided differences (Cheney & Kincaid, 2008).

Specifying the form of the interpolation polynomial in the multivariate case however, is extensively more difficult. It requires even more of a basis in core topics in linear algebra - namely the concepts of abstract vector spaces and subspaces.

In their survey of multivariate interpolation methods, Gasca et al. provide a general definition of the Lagrange form for the multivariate interpolating polynomial:

<img src="https://latex.codecogs.com/gif.latex?p(x)=\sum_{i=0}^{n}{y_i\prod_{j=0,j\neq{i}}^{n}\frac{H_j(x)}{H_j(x_i)}}"/>

Where for each point x_i, a hyperplane H_i containing x_i but no other interpolation nodes x_j is chosen. This geometric interpretation of the Kronecker delta allows us to expand the Lagrange polynomial to a higher dimensional space.

The major caveat with the interpolating problem in multiple variables is there doesn't appear to be a universal guarantee that there is a unique solution for a given dataset. In their 1990 article, Carl de Boor and Amos Ron emphasize that generalizing the interpolation problem to multiple dimensions has a tendency to not preserve all favorable properties. That said, in this project, we learned you could guarantee a unique solution for datasets of certain sizes by reframing the interpolation problem in the context of solving a linear system.

## Computational Approach: ##

In order to find an interpolating polynomial for a given dataset, we decided to utilize Vandermonde matrices. These matrices are specially designed to find the coefficients of a general multivariate polynomial by solving a linear system of equations. Below we specify how the Vandermonde matrix V is constructed from a given dataset.

Given N^2 - dimensional data vectors X, Y, and Z, we are able to define a matrix-valued function V(X, Y) such that V(X, Y) * C = Z for some unknown vector of coefficients C.

Due to it being unclear whether a multivariate interpolation problem will yield a unique solution, we decided that it would be best to ensure that our Vandermonde matrix is square (N^2 x N^2) with no free variables. Due to the nature of bivariate polynomials, this is why it is of the utmost importance that the number of datapoints given to this function is a perfect square. This is the only way to guarantee that a unique solution will be found with our approach.

For the bivariate case, we first construct an N x N matrix of powers P in the form P[i, j] = x^i * y^j.

This matrix P contains all the terms of our multivariate interpolation polynomial. After all terms are calculated, P is then flattened from an N x N matrix to an 1 x N^2 row vector containing all of the same terms.

Our Vandermonde matrix V[X, Y] then uses the flattened form of power matrix P to build each row of the Vandermonde matrix V.

For example, given dataset X = [1, 2, 3, 4] and Y = [1, 2, 3, 4], our power matrix would be a 2 x 2 matrix of the form:

```
[ [ x^0 * y^0, x^0 * y^1],
  [ x^1 * y^1, x^1 * y^1] ]
```

Its flattened form would be:

```
[x^0 * y^0, x^0 * y^1, x^1 * y^0, x^1 * y^1]
```

And V(X, Y) would be a 4 x 4 matrix with the following values:

```
[ [ 1, 1, 1, 1 ],
  [ 1, 2, 2, 4 ],
  [ 1, 3, 3, 9 ],
  [ 1, 4, 4, 16] ]
 ```
 
Once V has been constructed in this fashion, the unknown coefficient vector C can be solved for using a standard Gaussian Elimination / Row Reduction procedure. We use the version with partial pivoting specified in Cheney and Kincaid, 6th ed. on pages 267-269.

## Contents: ##
  - vmond.py contains the Vandermonde() function along with the Gauss() function it calls; these are used to calculate the coefficient vector used by bivariate.py
  - bivariate.py combines the terms from a coefficient vector produced by Vandermonde() with their corresponding (x^n)(y^n) terms (this is dependent on the size of the data set)
  - graph_bivariate.py contains the graphing functionality of the project. It uses bivariate.py to generate appropriate z-values for a generated range of data points (including those which are being interpolated) and then graphs them using plot_trisurf() from Axes3D. It overlays the interpolating function with the original points that are being interpolated.
  - multivar_interp_generic.py is a file that was used to visualize the power matrix and terms of the bivariate function. It has no functional use in EXAMPLE.py
  - EXAMPLE.py the working example which contains all of the code files compiled along with a series of examples. The intention is that these be run in a jupyter notebook using the IPython kernel (this allows for the interactive graphing).
 
PolySolvy-Version-x---Bivariate.py contains all the individual files used by Working_Example.py to interpolate and graph several examples of bivariate functions.
  
## How To Use EXAMPLE.py: ## 
  - Copy all the code in EXAMPLE.py
  - Open web-browser tab and navigate to https://jupyter.org/try
  - Once there, click on "Try Classic Notebook" and then wait for the page to load
  - Once the page has loaded click on "insert" tab at the top menue, from there click on "Insert Cell Below".
  - Paste the entirety of the copied code from EXAMPLE.py into the "In [ ]" cell that appears
  - Hit "Run" on the top menu and enjoy.
    
    *NOTE: to try other examples enter them at the bottom of the cell (near the other examples) and then press "Run" again. Input must follow the same patterning as the other examples and contain a "square" amount of data points (i.e. 1, 4, 9, etc.). The limitation of using jupyter notebook is that it can't handle running more than 9 points currently. Further development for sets of 16+ is needed.
  
 
## References: ##

- de Boor, C., & Ron, A. (1990). On multivariate polynomial interpolation. Constructive Approximation, 6(3), 287–302. https://doi.org/10.1007/BF01890412
- E. Ward Cheney and David R. Kincaid. 2007. Numerical Mathematics and Computing (6th. ed.). Brooks/Cole Publishing Co., USA.
- Gasca, M., & Sauer, T. (2000). Polynomial interpolation in several variables. Advances in Computational Mathematics, 12(4), 377. https://doi.org/10.1023/A:1018981505752
- University of Waterloo Electrical Engineering: https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/05Interpolation/multi/
- Wikipedia, the free encyclopedia, https://en.wikipedia.org/wiki/Main_Page
