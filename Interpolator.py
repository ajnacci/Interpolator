import sys

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    sys.exit(1)

import math

import pkg_resources
pkg_resources.require("numpy==1.16.5", "mpmath==1.1.0")

import numpy as np
import mpmath

# Generates points on a region given a distribution scheme.
def point_generator(N, distribution="even", extra_args=dict()):
  xi = None
  if(distribution == "even"):
    interval_lower = extra_args['interval_lower']
    interval_width = extra_args['interval_width']
    xi = [interval_lower + (interval_width*i)/N for i in range(N+1)]
  if(distribution == "cotangent"):
    xi = [float(mpmath.cot((i+(math.copysign(1, i-(N+1)/2)-N)/2)*math.pi/(N+1))/2) for i in range(N+1)]
  return sorted(xi)

# Efficiently calculate each f[x_0 ... x_v]
def memo_c(xi, f, N):
  c = np.zeros((N+1,N+1), dtype="float64")
  for i in range(N+1):
    c[i][0] = float(f(xi[i]))
  
  #m, n-m
  for i2 in range(N): #n-m-1
    i = i2+1 #n-m
    for j in range(N-i+1): #m, capped where m+n-m > N, n > N since not needed
      c[j][i] = float((c[j+1][i-1] - c[j][i-1])/(xi[j+i]-xi[j]))
      
  return c

# Define a whole bunch of polynomial functions. Polynomials will be
# represented as lists of coefficients.

# Multiply a function by a scalar
def const_mult_funct(f, c):
  fn = len(f)
  h = list()
  for i in range(fn):
    h.append(f[i]*c)
    
  return h

# Multiply two functions together
def multiply_funct(f, g):
  fn = len(f)
  gn = len(g)
# x^2 + x + 1 all squared, [1,1,1]*[1,1,1]
# fn = gn = 3
  h = list()
  for i in range(fn+gn-1): # 3 + 3 - 1 = 5
    h.append(0) # h = [0]
    j = min(fn-1,i) # j = min(2, 0) = 0
    while i - j < gn and j >= 0: #while 0 - 0 < 3
      h[i] = h[i] + f[j]*g[i-j] # h[0] = h[0] + f[0]g[0-0]
      j = j-1 #
          
  return trim_function(h)

# Generate f^0, f^1, ..., f^n and return these polynomials
# in a list such that result[k] = f^k
def exponentiate_list_funct(f, n):
  acc = list()
  acc.append(list())
  acc[0].append(1)
  for i2 in range(n):
    i = i2+1
    acc.append(multiply_funct(acc[i2], f))
    
  return acc

# Just return f^n
def exponentiate_funct(f, n):
  return exponentiate_list_funct(f, n)[n]

# Find f+g and return it
def add_funct(f, g):
  fn = len(f)
  gn = len(g)
  
  if fn >= gn:
    h = f
    for i in range(gn):
      h[i] = h[i] + g[i]
  else:
    h = g
    for i in range(fn):
      h[i] = h[i] + f[i]
  
  return trim_function(h)

# Find f(g(x)) and return it
def compose_funct(f, g):
  fn = len(f)
  ex = exponentiate_list_funct(g, fn)
  h = list()
  h.append(0)
  for i in range(fn):
    h = add_funct(h, multiply_funct([f[i]], ex[i]))
    
  return trim_function(h)

# Find f(x - r) and return it
def shift_funct(f, r):
  return compose_funct(f, [-r, 1])

# Calculate f(x) for a given value x
def eval_funct_at(f, x):
  fn = len(f)
  acc = 0
  for i in range(fn):
    acc = acc + f[i]*(x**i)
    
  return acc

# Find the derivative of a polynomial f. This was only used
# for testing purposes.
def derivative(f):
  h = list()
  fn = len(f)
  for i2 in range(fn-1):
    i = i2+1
    h.append(i*f[i])
    
  return h

# Trim leading 0 coefficients
def trim_function(f):
  fn = len(f)
  h = list()
  hn = fn
  for i in range(fn):
    k = fn-1-i
    if f[k] == 0.0:
      hn = hn - 1
    else:
      break
  
  for i in range(hn):
    h.append(f[i])
  
  return h

# Print a readable version of the polynomial f
def display_funct(f2):
  f = trim_function(f2)
  fn = len(f)
  acc = ""
  for i in range(fn):
    k = fn-i-1
    V = f[k]
    if i != 0:
      V = abs(V)
    acc = acc + str(V)
    if i != fn-1:
      acc = acc + "x"
    if k > 1:
      acc = acc + "^" +str(k)
    if i != fn-1:
      if f[k-1] >= 0:
        acc = acc + " + "
      else:
        acc = acc + " - "
  print(acc)
  
# Generate the nth term of the newton form
# for polynomial interpolation
def newton_part(c, xi, n):
  h = list()
  h.append(1)
  for i in range(n):
    h = multiply_funct(h,[-1*xi[i], 1])
    
  h = const_mult_funct(h, c[0][n])
  return h

# Generate the Newton form polynomial interpolant given
# f[x_i...x_j], a list of points, and the degree
# of the interpolant.
def full_newton(c, xi, N):
  h = list()
  h.append(0)
  
  for i in range(N+1):
    h = add_funct(h, newton_part(c, xi, i))
    
  return h

# Generate the nth term of the lagrange form
# for polynomial interpolation
def lagrange_part(xi, n, N):
  h = list()
  h.append(1)
  for i in range(N+1):
    if i != n:
      h = multiply_funct(h, [-1*xi[i], 1])
      h = const_mult_funct(h, 1/(xi[n]-xi[i]))
  
  return h

# Generate the Lagrange form polynomial interpolant given
# f, a list of points, and the degree
# of the interpolant.
def full_lagrange(f, xi, N):
  h = list()
  h.append(0)
  
  for i in range(N+1):
    h = add_funct(h, const_mult_funct(lagrange_part(xi, i, N), f(xi[i])))
    
  return h

# Generate the H-F form polynomial interpolant's nth
# term using f, a list of points, and N for degree (2N+1).
# This function only computes the terms which are concerned
# with the value of f
def hermite_value_part(f, xi, n, N):
  lag = lagrange_part(xi, n, N)
  lagrange_der_val = eval_funct_at(derivative(lag), xi[n])
  h = const_mult_funct([-1*xi[n], 1],-2*lagrange_der_val)
  h = add_funct(h, [1])
  h = multiply_funct(h, exponentiate_funct(lag, 2))
  return h

# Generate the H-F form polynomial interpolant's nth
# term using fd (derivative of f), a list of points,
# and N for degree (2N+1). This function only computes 
# the terms which are concerned with the value of f
def hermite_der_part(fd, xi, n, N):
  lag = lagrange_part(xi, n, N)
  h = multiply_funct([-1*xi[n], 1], exponentiate_funct(lag, 2))
  return const_mult_funct(h, fd(xi[n]))

# Generate the H-F polynomial given that HF'=0 at each
# given point. Requires the function to approximate f, 
# the list of points xi, and N which corresponds to a
# degree 2N+1 polynomial
def hermite_fejer_noslope(f, xi, N):
  h = list()
  h.append(0)
  for i in range(N+1):
    h = add_funct(h, const_mult_funct(hermite_value_part(f, xi, i, N), f(xi[i])))
    
  return h

# Return a function which can be used to generate an array representing
# theta minus in the mass matrix
def theta_minus(xi,N):
  xi = np.asarray(xi, dtype="float64")
  return lambda k : np.where(k==0, 0, np.where(k<N, (xi[k+1]-xi[k])/(xi[k+1]-xi[k-1]), -100))

# Return a function which can be used to generate an array representing
# theta plus in the mass matrix
def theta_plus(xi,N):
  xi = np.asarray(xi, dtype="float64")
  return lambda k : np.where(k==0, 0, np.where(k<N, (xi[k]-xi[k-1])/(xi[k+1]-xi[k-1]), -100))

# Return a function which can be used to generate the mass matrix
def linear_system_gen(xi, N):
  tmf = np.vectorize(theta_minus(np.array(xi),N))
  tpf = np.vectorize(theta_plus(np.array(xi),N))
  tm = np.fromfunction(tmf, (N+1,), dtype="int32")
  tp = np.fromfunction(tpf, (N+1,), dtype="int32")
  return lambda r, c : np.where(r==c, 2/3, np.where((r==0 or r==N) and (abs(r-c)==1), 1/3, np.where(r==c-1, tp[r]/3, np.where(r-1==c, tm[r]/3, 0))))
  
# Return a function which can be used to generate an array representing
# the right hand side of the linear equation we solve for the slopes.
def r_gen(xi, c, N):
  tmf = np.vectorize(theta_minus(np.array(xi),N))
  tpf = np.vectorize(theta_plus(np.array(xi),N))
  tm = np.fromfunction(tmf, (N+1,), dtype="int32")
  tp = np.fromfunction(tpf, (N+1,), dtype="int32")
  return lambda k : np.where(k==0, c[0][1], np.where(k<N, tm[k]*c[k-1][1]+
                                      tp[k]*c[k][1], c[N-1][1]))

# Generate the cubic spline in the region between x[i] and x[i+1]
def gen_spline_part(f, s, xi, n):
  g = [-1*xi[n], 1] #x - x_nu
  h = [-1*xi[n+1], 1] #x - x_(nu+1)
  
  dx = xi[n+1] - xi[n]
  
  # acc1 = 2(x-x_nu)/(x_(nu+1) - x_nu)
  acc1 = const_mult_funct(g, 2/dx)
  acc1 = add_funct(acc1, [1])
  acc1 = const_mult_funct(acc1, f(xi[n]))
  acc1 = add_funct(acc1, const_mult_funct(g, s[n]))
  h_temp = const_mult_funct(h, 1/dx)
  other = exponentiate_funct(h_temp, 2)
  acc1 = multiply_funct(acc1, other)
  
  # acc2 = -2(x-x_(nu+1))/(x_(nu+1) - x_nu)
  acc2 = const_mult_funct(h, -2/dx) 
  acc2 = add_funct(acc2, [1])
  acc2 = const_mult_funct(acc2, f(xi[n+1]))
  acc2 = add_funct(acc2, const_mult_funct(h, s[n+1]))
  g_temp = const_mult_funct(g, 1/dx)
  other = exponentiate_funct(g_temp, 2)
  acc2 = multiply_funct(acc2, other)
  
  return add_funct(acc1, acc2)
  
# Generate a list of the cubic form of the spline over each
# region. spline_list[i] represents the spline on [x[i], x[i+1]]
def generate_splines(f, s, xi, N):
  spline = lambda n : gen_spline_part(f,s,xi,n)
  
  spline_list = list()
  for i in range(N):
    spline_list.append(spline(i))
    
  return spline_list
  
def determine_coefficients(f, N, process, dist, extra_args=dict()):
  xi = point_generator(N, distribution=dist, extra_args=extra_args)
  if process == "Newton":
    # c[i][j] will represent f[x_(i)...x_(i+j)]. We calculate it efficiently
    # using memoization all at once.
    c = memo_c(xi, f, N)
    
    # returns a "polynomial", which is just a list of coefficients
    # that represents the newton form interpolant.
    return full_newton(c, xi, N)
  
  if process == "Lagrange":
    return full_lagrange(f, xi, N)
    
  if process == "Hermite Fejer 0 slope":
    return hermite_fejer_noslope(f, xi, N)
  
  if process == "Natural cubic spline":
    # Needed for some numpy functions, but won't be used in calculations
    xi.append(-100)
    
    # c[i][j] will represent f[x_(i)...x_(i+j)]. We calculate it efficiently
    # using memoization all at once.
    c = memo_c(xi, f, N)

    # Generate the mass matrix
    lsg = np.vectorize(linear_system_gen(xi,N))
    linsys = np.fromfunction(lsg, (N+1, N+1), dtype="int32")
    print(linsys)
    
    # Generate the RHS
    rr = np.vectorize(r_gen(xi, c, N))
    b = np.fromfunction(rr, (N+1,), dtype="int32")
    
    # Solve for the s vector of slopes
    s = np.linalg.solve(linsys, b)
  
    # Generates a list of "polynomials", which are just stored as lists of 
    # coefficients.
    # spline[i] gives the interpolant between xi[i] and xi[i+1]
    # spline will thus have N parts
    spline = generate_splines(f, s, xi, N)
    return spline
    
# Estimate the maximum error between f and pn on the specified interval.
# This code generally works best for evenly spaced points.
def estimate_error(f, pn, N, ir=2, il=-1, sensitivity=4, iterations=10):
  vals = list()
  interval_range = ir
  interval_lower = il
  sig_point = 0
  sig_value = 0
  # test N*sensitivity evenly distributed points, then subdivide the interval
  # around the point with the largest magnitude and repeat
  for _ in range(iterations):
    test_points = [(interval_range*i)/(sensitivity*N) + interval_lower 
                   for i in range(sensitivity*N)]

    for ele in test_points:
      val = pn(ele) - f(ele)
      if abs(sig_value) < abs(val):
        sig_point = ele
        sig_value = abs(val)
    
    # I put these in a list for testing purposes
    vals.append((sig_point, sig_value))
    
    interval_range = 2*interval_range/(sensitivity*N)
    interval_lower = sig_point - interval_range/2
    
  return vals[iterations-1]