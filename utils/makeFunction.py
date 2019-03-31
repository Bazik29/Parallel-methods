import math

import mpmath
import scipy
import scipy.integrate as integrate
import scipy.optimize

from jinja2 import Environment, FileSystemLoader, Template

import sympy
from sympy.printing.cxxcode import cxxcode
from sympy.utilities.lambdify import lambdify

x = sympy.Symbol('x')
y = sympy.Symbol('y')

foo_a = 1
foo_b = 2
f_foo = x * (sympy.pi/2 - sympy.atan(x/2))

bar_a = 1
bar_b = 3
bar_c = 0
bar_d = 3

condition = "((y < x) && (y > x / 2))"
f_bar = pow(x + y, 2) / x

foo_res = integrate.quad(lambdify(x, f_foo, 'scipy'), foo_a, foo_b)[0]

def tmp(x_, y_):
    if (y_ < x_) and (y_ > x_ / 2):
        return pow(x_ + y_, 2) / x_
    else:
        return 0

bar_res = integrate.dblquad(tmp, bar_a, bar_b, lambda q: bar_c, lambda q: bar_d)
bar_res = 13.36 #?

foo_x = f_foo.diff(x)
foo_xx = f_foo.diff(x, 2)
foo_xxx = f_foo.diff(x, 3)
foo_xxxx = f_foo.diff(x, 4)

foo_max_1 = -abs(foo_x)
foo_max_1 = lambdify(x, foo_max_1, 'scipy')
foo_max_1 = scipy.optimize.fminbound(foo_max_1, foo_a, foo_b)

foo_max_2 = -abs(foo_xx)
foo_max_2 = lambdify(x, foo_max_2, 'scipy')
foo_max_2 = scipy.optimize.fminbound(foo_max_2, foo_a, foo_b)

foo_max_3 = -abs(foo_xxx)
foo_max_3 = lambdify(x, foo_max_3, 'scipy')
foo_max_3 = scipy.optimize.fminbound(foo_max_3, foo_a, foo_b)

foo_max_4 = -abs(foo_xxxx)
foo_max_4 = lambdify(x, foo_max_4, 'scipy')
foo_max_4 = scipy.optimize.fminbound(foo_max_4, foo_a, foo_b)

bar_x = f_bar.diff(x)
bar_y = f_bar.diff(y)
bar_xx = f_bar.diff(x, x)
bar_yy = f_bar.diff(y, y)
bar_xy = f_bar.diff(x, y)

foo = cxxcode(f_foo, standard='C++11')
foo_x = cxxcode(foo_x, standard='C++11')
foo_xx = cxxcode(foo_xx, standard='C++11')
foo_xxx = cxxcode(foo_xxx, standard='C++11')
foo_xxxx = cxxcode(foo_xxxx, standard='C++11')

bar = cxxcode(f_bar, standard='C++11')
bar_x = cxxcode(bar_x, standard='C++11')
bar_y = cxxcode(bar_y, standard='C++11')
bar_xx = cxxcode(bar_xx, standard='C++11')
bar_yy = cxxcode(bar_yy, standard='C++11')
bar_xy = cxxcode(bar_xy, standard='C++11')

dict_names = {
    "foo_a": foo_a,
    "foo_b": foo_b,
    "foo": foo,
    "foo_x": foo_x,
    "foo_xx": foo_xx,
    "foo_xxx": foo_xxx,
    "foo_xxxx": foo_xxxx,
    "bar_a": foo_a,
    "bar_b": bar_b,
    "bar_c": bar_c,
    "bar_d": bar_d,
    "condition": condition,
    "bar": bar,
    "bar_x": bar_x,
    "bar_y": bar_y,
    "bar_xx": bar_xx,
    "bar_yy": bar_y,
    "bar_xy": bar_xy,
    "foo_max_1": foo_max_1,
    "foo_max_2": foo_max_2,
    "foo_max_3": foo_max_3,
    "foo_max_4": foo_max_4,
    "foo_res": foo_res,
    "bar_res": bar_res
}

env = Environment(loader=FileSystemLoader('./utils'))
template = env.get_template('function.tmpl')

my_file = open("function.h", "w")
my_file.write(template.render(dict_names))
my_file.close()

print("Done")
