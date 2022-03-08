from sympy import *
import math

x = symbols('x')
y = symbols('y')

x_range = (x,0,1)
y_range = (y,0,1)

print "f(x,y) = x   : ", integrate(x, x_range, y_range)
print "f(x,y) = y   : ", integrate(y, x_range, y_range)
print "f(x,y) = x*y : ", integrate(x*y, x_range, y_range)





