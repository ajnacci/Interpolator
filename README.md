# Interpolator
A Python library for function interpolation and approximation. It focuses on flexibility and intuitivity, and can be used to graph, generate, and evaluate various interpolation forms from data or from functions.

# Usage
Currently, there are only 3 methods that the user will typically call.
* **determine_coefficients**(f, N, process, dist="even", extra_args=dict())

    This method returns a list **L** that represents the coefficients of the interpolated polynomial. For example, if **L = [1,2,3]**, then our interpolated polynomial is **1 + 2x + 3x^2**.
    
    The argument **f** represents a python function that will take as input any of the "x" values for interpolation and give the corresponding value f(x). This is the function that we are trying to approximate, although we only care about its values at the "x" interpolation points, so you can fudge the function definition when needed. The argument **N** represents different things for different interpolation methods, so see the details on the methods in the next section to determine what value to provide. The optional argument **dist** refers to the distribution of the "x" interpolation points. You can read more about that in the next section. The optional argument **extra_args** is a dictionary that maps the name of any extra arguments (as a string) to the value of the extra argument. Extra arguments are dependent on the method of interpolation and point distribution, so check the next two sections for details.
    
    This method doesn't raise any errors at the moment unless they arise from the numpy or math libraries.
    
 * **estimate_error**(f, pn, N, ir=2, il=-1, sensitivity=4, iterations=10)
 
     This method returns a tuple whose second value is the approximation of the max absolute difference between the functions f and pn, and whose first value is the "x" at which this max absolute difference is achieved.
     
     This method works by sampling **sensitivity** \* N evenly spaced points, then iteratively subdividing the interval formed by the maximum found error point's neighbors. Technically, it can be run on any interpolation-function pair, although I suggest you only use it on polynomial interpolants with evenly spaced sample points.
     
     The argument **f** represents a python function that maps x to f(x) over the entire interval **[il, il+ir]**. Note that this requires us to know **f**, so we have stricter rules about when we can run this error estimation. The argument **pn** represents a python function that maps x to the interpolated pn(x) over the entire interval **[il, il+ir]**. The optional argument **ir** represents the range of the closed interval we restrict x to, and the optional argument **il** represents the lower boundary of this interval. The optional argument **sensitivity** represents how many subintervals we construct between each point before selecting one to descend into. Increasing this value should generally increase the accuracy and the rate of convergence of the method. The optional argument **iterations** represents how many times the method will descend before returning its best estimate.
     
     This method doesn't raise any errors at the moment unless they arise from the numpy or math libraries.
     
* **eval_funct_at**(f, x)

    This method returns a value that a polynomial (represented as a list of coefficients) takes at a given point x.
    
    The argument **f** represents the polynomial to evaluate, and the argument **x** is the point at which to evaluate it.
    
    This method doesn't raise any errors at the moment.
 
 # Interpolants
 When calling **determine_coefficients**, pass in the desired interpolation method as a string from the list below.
 * **Newton**- Finds the Newton form polynomial interpolant of *degree N*.
 * **Lagrange**- Finds the Lagrange form polynomial interpolant of *degree N*.
 * **Hermite-Fejer 0 slope**- Finds the Hermite-Fejer polynomial interpolant of *degree 2N + 1* given only values of a function, and assuming that the derivative should be 0 at each provided point.
 * **Natural cubic spline**- Finds the natural cubic spline (ie. the cubic spline with second derivative 0 at its two extreme given points x_0 and x_N) over *N intervals* bounded by a total of *N + 1 points*.
 
 # Point distributions
 When calling **determine_coefficients**, pass in the desired point distribution as a string from the list below. Note that passing in a value **N** generally creates **N + 1** points. In version 0.1, there is no user-provided point distribution. However, this will be added soon.
 * **Even**- The argument **extra_args** must contain two extra arguments in order to run this. They are 'interval_lower' (the lower bound for the interval and the value of x_0) and 'interval_width' (the size of the interval). The points are then evenly selected throughout the interval, with x_k = 'interval_lower' + k \* interval_width/N.
 * **Cotangent**- Currently, this will always distribute points on [-1,1], and requires N to be an odd number. It calculates x_k = (1/2)\*cotangent(k\*pi/(N + 1)) for k = [+/-] 1, 2, ..., (N+1)/2.
     
 # Versions
 ## Version 0.1 - Current
 In this version, the package has limited functionality and still needs to be organized better. The goal for Version 0.2 is to add data-point direct entry, organize the methods into sub-modules, implement graphing, and perhaps improve the usefulness of the function representations beyond just a list of coefficients. This will help improve error reporting, as well.
