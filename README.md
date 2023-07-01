# greens_examples
Example 2d integration problems using Green's Theorem with arbitrary domains and integrands

Two example problems using the "C" shaped domain with integrands f(x,y)=x^2+y^2 and f(x,y)=cos(1/2(x^2+y^2)).
The exact solutions are 35pi and 7pi/4(sin(9/2)-sin(1/2)) respectively.
The double integral is converted to a line integral using Green's Theorem.
The new integrand is calculated by numerically solving an ODE with arbitrary initial condition using Euler's method.
The ODE approximation is given on a grid and then interpolated.
The line integral is approximated with equally spaced quadrature.
There is a first order finite difference approximation used.

getpaths(n1,n2,n3,n4) generates the paths that define the domain with ni points for the ith path

plotdomain_ex1() and plotdomain_ex2() plot the domain and the integrand for the two examples

plot_integrand_ex1() and plot_integrand_ex2() plots the converted integrand for the line integral
The integrand is approximated on a grid using Euler's method to solve an ODE with arbitrary initial conditions.
The grid approximation is interpolated. The grid and interpolation are plotted together.

greens_ex1(n1,n2,n3,n4,n) and greens_ex2(n1,n2,n3,n4,n) approximate the integral and report relative error.
ni points are used for the quadrature on the ith path. 
n points are used in the x and y directions to approximate the ODE to convert the integrand.
A first order finite difference approximation is used.
Overall convergence is first order. 
n1=n2=n3=n4=n=500 is a good initial choice to try the code.
