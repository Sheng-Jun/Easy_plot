# Easy_plot
A plot script of explicit functions: ℝ → ℝ  

I wrote this script as a shell command for quickly plotting some functions,  
which could plot linear, semilogx, semilogy, loglog scale in Cartesian coordinate, and also could plot in counterclockwise/clockwise Polar coordinate with different offset of the origin point.

Moreover, I wrote it as def Easy_Plot() in script, so people could write a python script and import this script to plot diagrams. In this way, people could assign the title, labels for x and y-axes, even Latex form string and assign the style of lines. Besides drawing normal functions, people could directly input an array of domain data and serial arrays of range data to plot diagrams!

[The manual for using as a shell command]

python easy_plot.py [-d] [-D] [-r] [-R] [-P] [-C] [-o] [-p] [-e] [-i] [-A] [-N] 'equation 1 of x' 'equation 2 of x'...\n
                                                         (several functions [in numpy] marked by single/double qutoes)
 
Options:
-d a,b
  Let the domain to be [a,b]; the default is [-10,10].
-D a,b,B
  Similar to '-d', but the log scale of the base B, the default for B is 10.
-r c,d
  Let the range to be [c,d]; the default is [-10,10].
  If input -r'auto' then it'll be auto scaled.
-R c,d,B
  Similar to '-r', but the log scale of the base B, the default for B is 10.
  If input -R'auto' then it'll be auto scaled.
-P c,d
  Polar coordinate in the counterclockwise direction with the range = [c,d].
  If input -P'auto' then it'll be auto scaled.
-C c,d
  Similar to '-P', but in the clockwise direction.
-o t
  Let the offset for the origin be t radians, the default is 0.
-p N
  Let the size of the partition of the domain to be N; if not, the default is 500.
-e p,q...
  Exclude points p, q... in the domain.
-i n
  Let 'inf' and '-inf' to be n, otherwise printing a warning.
-A
  Print all function values of the last one fuction, the 2nd contains abs(negative), if semilogy/loglog.
-N
  Don't print the figure.
-h
  Print this help information.
This script provided some constants, please look them up.

e.g.
sj@machine ~$ python easy_plot.py -d 0,5**2 -r -5,5 -p 200*2 -e 0,2*10 1/x '-e**(1/x)' '5./2*cos(x/pi)'

The manual for using as a python module

[fig, ax, domain, func1, func2] = Easy_Plot(equ_list, *args, **kwargs)
This function will return a list: [fig, ax, domain, func1, func2].
  A) fig is a matplotlib.figure.Figure object, but if user inputs their
    own AxesSubplot, fig will be a float number 0.
  B) ax is a matplotlib.axes.AxesSubplot object.
  C) domain is a array of points of the domain.
  D) func1 is a list of arrays of function values of equ_list.
    But the values beyond the range (when [a]semilogy/loglog, and [b]polar)
    will be modified, and the original one will be stored in
    the corresponding array in list func2.
    In case [a], func2 contains abs(negative function value).
    In case [b], func2 contains negative function values.
    Otherwise the corresponding items in func2 will be a string, 'useless'.

Some other modules also be imported:
  from numpy import *
  from scipy import special
  import matplotlib.pyplot as plt
  from matplotlib import rc
  from getopt import getopt
  from sys import argv

Here add some variables, and Latex strings are allowed:
Necessary
equ_list
 A list, which can contain strings are functions of x, or/and arrays(numpy) that contains function values. And you can only give a string or an array, if just one function. If equ_list just are some arrays, you can directly see the shape without giving domain and range, but the dimension of them must be the same.

args
equ_label
 A list of label of functions, e.g.['$x$']; otherwise they'll be ordinal numbers.
style
 Cartesian coordinate: 0: linear plot(default), 1: semilogx, 2: semilogy, 3: loglog.
 Polar coordinate: 4: counterclockwise, 6: clockwise.
  If some part is beyond the range, it'll be marked by another label automatically.
'debug'
 Print the debug info.

kwargs
[The default domain = [-10,10], range = [-10,10], and base of log-scale = 10]
x/domain = [a, b]
 Let the domain be [a,b]. If log-scale, the base will be 10.
x/domain = [a, b, c]
 Let the domain be [a,b], and the base of log-scale be c.
x/domain = an array(numpy) of data points. In this case, the x-scaling is auto.
x/domain = [(an array), a, b]
 Similar to the above, but limit the domain to be [a, b]. If log-scale, the base will be 10.
x/domain = [(an array), c]
 Similarly, but just assign the base of log-scale.
x/domain = [(an array), a, b, c]
 Similarly, but assign the limitation and the base.
y/range = [a, b], y/range = [a, b, c] All similar to the above.
y/range = 'auto'
 The y-scaling is auto. But if equ_list just are some arrays, the y-scaling is already auto.
y/range = ['auto', c]
 The y-scaling is auto, and assign the base of log-scale.
offset = t
 Let the offset for the origin be t radians, the default is 0. Operators are allowed.
ax = a matplotlib.axes.AxesSubplot object
 By this, user create a figure with a set of subplots, then use Easy_Plot() to plot one of them.
 For example,
     fig, (ax1, ax2) = subplots(2, 1, sharex=True)
     a = Easy_Plot('x', ax=ax1,...)
     b = Easy_Plot('1/x', ax=ax2,...)
ls = a string if one plot, or a list of strings if more than one plot.
 Assign the style, or/and color of lines in python form, like 'ro' for one plot or ['ro', 'k-'] for two plots. If doesn't assign or assign '', means set by python. But the automatic colors for different parts won't work.
color = a string if one plot, or a list of strings if more than one plot.
 Besides 'red' and 'blue' etc., HTML color codes are also fine.
 But the user is suggested not to use ls and color to assign the color of lines at the same time.
alpha = a float number between 0. and 1. if one plot, or a list of float numbers if more than one plot.
loc = an integer.
 Set the location of legend.
p/partition = N
 Let the size of the parition of the domain to be N; If not, the default is 500. Operators are allowed.
ex/exclusion = a number, or a list of numbers.
 A list of points to be excluded. Operators are allowed.
i/inf = n
 If don't want to get a warning when the value is 'inf' or '-inf', then just let them to be n. Operators are allowed.
title = '...', xlabel = '...', ylabel = '...'
 Assign the title, xlabel, and ylabel. But the labels is only for the Cartesian.

Constant values list: print(const)

e.g.
result = Easy_Plot(['100*(10**4/(10**8+x**2)**0.5)**5'],  ['$Plummer\'s\ model\ with\ 5$'],  3,
                              x = [10**-2, 10**10],  y = [10**-4, 10**3],
                              title = '$An\ empirical\ model$',  xlabel = '$Radius$', ylabel = '$Density$')
result2 = Easy_Plot('x/2/pi',  '$spiral$',  6,  x = [0, 4*pi],  y = [0, 2],  ls = 'r',  offset = 1.5*pi,  title = '$Sample$')
result3 = Easy_Plot(['abs(x)',  array([2.1, 2.5, 1.2, 3.5])],  ['test', 'data'],  4,
                              x = array([-0.5*pi, pi, 0.5*pi, 0.75*pi]),  y = [0, 4],  ls = ['k-', 'ro'])
result4 = Easy_Plot(['1/x',  'sin(1/x)'],  ['$\\frac{1}{x}$',  'oscillator'],  2,  x = [-10, 10],  y = [0.1, 1.3],  ls = ['k-', ''])
