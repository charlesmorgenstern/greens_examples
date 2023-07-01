#####################################################
#####################################################
using Plots
using Interpolations
#####################################################
#####################################################
function getpaths(n1,n2,n3,n4) #generate four paths that define the c domain
#paths defined in counterclockwise orientation for Green's Thm

paths=Dict();

t=collect(range(pi/4,2*pi,n1)); #outer circle of radius 3
x=3*cos.(t)
y=3*sin.(t)
p1=Matrix{Float64}(undef,n1,2)
p1[:,1]=x
p1[:,2]=y
paths["p1"]=p1;


x=collect(range(3,1,n2));  #horizontal line
p2=Matrix{Float64}(undef,n2,2)
p2[:,1]=x
p2[:,2].=0.0
paths["p2"]=p2;

t=t[end:-1:1] #inner circle of radius one
t=collect(range(2*pi,pi/4,n3));
x=cos.(t)
y=sin.(t)
p3=Matrix{Float64}(undef,n3,2)
p3[:,1]=x
p3[:,2]=y
paths["p3"]=p3;


x=collect(range(cos(pi/4),3*cos(pi/4),n4)); #line of slope 1
p4=Matrix{Float64}(undef,n4,2)
p4[:,1]=x
p4[:,2]=x
paths["p4"]=p4;

return paths

end

###############################################################
###############################################################
#EXAMPLE 1: integrand is f(x,y)=x^2+y^2
################################################################
################################################################
function plotdomain_ex1() 
#plots the domain legend shows indexing of paths
#plots the integrand contours for ex1


#clear current plot so they don't overlay
Plots.CURRENT_PLOT.nullableplot = nothing

#get paths and plot them
paths=getpaths(100,100,100,100) 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]
fig=plot(p1[:,1],p1[:,2], linewidth=5)
plot!(p2[:,1],p2[:,2], linewidth=5)
plot!(p3[:,1],p3[:,2], linewidth=5)
plot!(p4[:,1],p4[:,2], linewidth=5)

# Define the function
f(x, y) = x^2 + y^2

# Create a grid of x and y values for plotting
x = range(-3.5, 3.5, length=100)
y = range(-3.5, 3.5, length=100)

# Evaluate the function on the grid
z = [f(i, j) for i in x, j in y]

# Create the contour plot
contour!(x, y, z, levels=100, color=:plasma, xlabel="x", ylabel="y", title="domain and integrand f(x,y)=x^2+y^2")

# show the plot
return fig
end

#####################################################
#####################################################

function plot_integrand_ex1()
#calculate the new integrand for the line integral
#using Euler's method to solve an ODE
#interpolate the integrand
#plot the grid from the ODE method and the interpolated function


#clear current plot so they don't overlay
Plots.CURRENT_PLOT.nullableplot = nothing

#define number of points in x and y for Euler's method
nx = 20
ny = 20

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get vectors of x and y values
x = range(-3.5, 3.5, length=nx)
y = range(-3.5, 3.5, length=ny)

#initialize arrays
z = Array{Float64}(undef, nx, ny)
xx = Array{Float64}(undef, nx, ny)
yy = Array{Float64}(undef, nx, ny)
xx[1,:] .= x[1]
yy[1,:] = y

#arbitrary initial condition to make sol'n unique
z[:,1] .= 0 

#euler's method
for i=2:nx
 for j=1:ny
  xx[i,j]=x[i]
  yy[i,j]=y[j]
  z[j,i]=z[j,i-1]+hx*f(x[i-1],y[j])
 end
end
z=z'

#interpolate grid
itp = LinearInterpolation((x, y), z)

#fine grid for plotting
x2 = range(extrema(x)..., length=300)
y2 = range(extrema(y)..., length=300)
z2 = [itp(x,y) for y in y2, x in x2]

#plot
fig=heatmap!(x2, y2, z2, title="Interpolated heatmap of new integrand")
scatter!(xx,yy,zcolor = z,legend=false)


end

################################################################
################################################################

function greens_ex1(n1,n2,n3,n4,n)  
#ni points for the ith path
#n segments in x and y for Euler's method to get integrand
#numeric approx of integral using Green's Thm 
#integral can be varified analytically by hand 
#to be 35pi. We make no assumptions of what the curves or integrand are.
#We only assume the points are equally spaced in some parameter t.
#We let t range from 0 to 1 for each curve.
#We assume the integrand is a function f(x,y) with continuous derivatives
#Actual integrand for this example is f(x,y)=x^2+y^2

#set parameters for euler's method
nx = n 
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=x.^2+y.^2

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for euler's method
z = Array{Float64}(undef, nx, ny)

#euler's method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  z[j,i]=z[j,i-1]+hx*f(x1[i-1],y1[j])
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line 
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

int1=0.0 #line integral for outer circle
for i=2:n1
x=p1[i,1]
y=p1[i,2]
dy=(p1[i,2]-p1[i-1,2])/h1
int1+=itp(x,y)*dy
end 
int1*=h1

int2=0.0 #line integral for horizontal line
for i=2:n2
x=p2[i,1]
y=p2[i,2]
dy=(p2[i,2]-p2[i-1,2])/h2
int2+=itp(x,y)*dy
end 
int2*=h2

int3=0.0 #line integral for inner circle
for i=2:n3
x=p3[i,1]
y=p3[i,2]
dy=(p3[i,2]-p3[i-1,2])/h3
int3+=itp(x,y)*dy
end 
int3*=h3


int4=0.0 #line integral for line w/ slope 1
for i=2:n4
x=p4[i,1]
y=p4[i,2]
dy=(p4[i,2]-p4[i-1,2])/h4
int4+=itp(x,y)*dy
end 
int4*=h4

int=(int1+int2+int3+int4)


exactarea=35*pi

display("Exact integral is:")
display(exactarea)
display("Numeric approx of Green's Thm integral is:")
display(int)
display("Relative error is:")
err=abs((exactarea-int)/exactarea)
display(err)
end

###############################################################
###############################################################
#EXAMPLE 2: integrand is f(x,y)=cos(1/2(x^2+y^2))
################################################################
################################################################
function plotdomain_ex2() 
#plots the domain legend shows indexing of paths
#plots the integrand contours for ex2


#clear current plot so they don't overlay
Plots.CURRENT_PLOT.nullableplot = nothing

#get paths 
paths=getpaths(100,100,100,100) 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

# Define the function
f(x, y) =cos(1/2*( x^2 + y^2))

# Create a grid of x and y values for plotting
x = range(-3.5, 3.5, length=100)
y = range(-3.5, 3.5, length=100)

# Evaluate the function on the grid
z = [f(i, j) for i in x, j in y]

# Create the contour plot and paths
fig=contour(x, y, z, levels=50, color=:plasma, xlabel="x", ylabel="y", title="domain and integrand f(x,y)=cos(1/2(x^2+y^2))")
plot!(p1[:,1],p1[:,2], linewidth=5,c=:blue)
plot!(p2[:,1],p2[:,2], linewidth=5)
plot!(p3[:,1],p3[:,2], linewidth=5)
plot!(p4[:,1],p4[:,2], linewidth=5)

# show the plot
return fig
end

#####################################################
#####################################################

function plot_integrand_ex2()
#calculate the new integrand for the line integral
#using Euler's method to solve an ODE
#interpolate the integrand
#plot the grid from the ODE method and the interpolated function

#clear current plot so they don't overlay
Plots.CURRENT_PLOT.nullableplot = nothing

#define number of points in x and y for Euler's method
nx = 20
ny = 20

#step size
hx = 7/(nx-1)

#define function
f(x,y)=cos(1/2*(x.^2+y.^2))

#get vectors of x and y values
x = range(-3.5, 3.5, length=nx)
y = range(-3.5, 3.5, length=ny)

#initialize arrays
z = Array{Float64}(undef, nx, ny)
xx = Array{Float64}(undef, nx, ny)
yy = Array{Float64}(undef, nx, ny)
xx[1,:] .= x[1]
yy[1,:] = y

#arbitrary initial condition to make sol'n unique
z[:,1] .= 0 

#euler's method
for i=2:nx
 for j=1:ny
  xx[i,j]=x[i]
  yy[i,j]=y[j]
  z[j,i]=z[j,i-1]+hx*f(x[i-1],y[j])
 end
end
z=z'

#interpolate grid
itp = LinearInterpolation((x, y), z)

#fine grid for plotting
x2 = range(extrema(x)..., length=300)
y2 = range(extrema(y)..., length=300)
z2 = [itp(x,y) for y in y2, x in x2]

#plot
fig=heatmap!(x2, y2, z2, title="Interpolated heatmap of new integrand")
scatter!(xx,yy,zcolor = z,legend=false)


end

################################################################
################################################################

function greens_ex2(n1,n2,n3,n4,n)  
#ni points for the ith path
#n segments in x and y for Euler's method to get integrand
#numeric approx of integral using Green's Thm 
#integral can be varified analytically by hand 
#to be 7pi/4(sin(9/2)-sin(1/2)). 
#We make no assumptions of what the curves or integrand are.
#We only assume the points are equally spaced in some parameter t.
#We let t range from 0 to 1 for each curve.
#We assume the integrand is a function f(x,y) with continuous derivatives
#Actual integrand for this example is f(x,y)=cos(1/2(x^2+y^2))

#set parameters for euler's method
nx = n 
ny = n

#step size
hx = 7/(nx-1)

#define function
f(x,y)=cos(1/2*(x.^2+y.^2))

#get x and y points
x1 = range(-3.5, 3.5, length=nx)
y1 = range(-3.5, 3.5, length=ny)

#initialize array for euler's method
z = Array{Float64}(undef, nx, ny)

#euler's method for solving ODE to get new integrand
z[:,1] .= 0 #arbitrary initial condition to make sol'n unique
for i=2:nx
 for j=1:ny
  z[j,i]=z[j,i-1]+hx*f(x1[i-1],y1[j])
 end
end
z=z'

#interpolate sol'n
itp = LinearInterpolation((x1, y1), z) #interpolate solution

#get paths for line integrals
paths=getpaths(n1,n2,n3,n4) #load paths with ni points 
p1=paths["p1"]
p2=paths["p2"]
p3=paths["p3"]
p4=paths["p4"]

h1 = 1/(n1-1) #spacing for parameter of outer circle
h2 = 1/(n2-1) #spacing for parameter of bottom line 
h3 = 1/(n3-1) #spacing for parameter of inner circle
h4 = 1/(n4-1) #spacing for parameter of top line

int1=0.0 #line integral for outer circle
for i=2:n1
x=p1[i,1]
y=p1[i,2]
dy=(p1[i,2]-p1[i-1,2])/h1
int1+=itp(x,y)*dy
end 
int1*=h1

int2=0.0 #line integral for horizontal line
for i=2:n2
x=p2[i,1]
y=p2[i,2]
dy=(p2[i,2]-p2[i-1,2])/h2
int2+=itp(x,y)*dy
end 
int2*=h2

int3=0.0 #line integral for inner circle
for i=2:n3
x=p3[i,1]
y=p3[i,2]
dy=(p3[i,2]-p3[i-1,2])/h3
int3+=itp(x,y)*dy
end 
int3*=h3


int4=0.0 #line integral for line w/ slope 1
for i=2:n4
x=p4[i,1]
y=p4[i,2]
dy=(p4[i,2]-p4[i-1,2])/h4
int4+=itp(x,y)*dy
end 
int4*=h4

int=(int1+int2+int3+int4)


exactarea=((7*pi)/4)*(sin(9/2)-sin(1/2))

display("Exact integral is:")
display(exactarea)
display("Numeric approx of Green's Thm integral is:")
display(int)
display("Relative error is:")
err=abs((exactarea-int)/exactarea)
display(err)
end

###############################################################
###############################################################

