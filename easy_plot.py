#! /usr/bin/env python
#-*- coding:utf-8 -*-
#Written on 1/12/2013 by Sheng-Jun Lin
#Version 2 on 19/12/2013, re-write some parts of the domain and range, and can input a matplotlib.axes.AxesSubplot object
#Version 3 on 15/03/2014, add functions for plotting 1-dimension arrays without giving domain, and set location of legend, and write more comments in code

#The basic package for this script are numpy and matplotlib.
#If you don't have the others, you can just comment them.
from __future__ import print_function
from scipy import special
import matplotlib.pyplot as plt
from matplotlib import rc
from getopt import getopt
from sys import argv
import copy
from const import *

#If you do/don't have command: "latex" in shell, you could uncomment/comment this line.
rc('text', usetex=True)

func = 'integral(equ, x0, x1, p)\n\
Bv(v, T)\n\
Bvc(v, T)\n\
Bl(lamda, T)\n\
Blc(lamda, T)\n\
Normal(x, mean, sd)\n\
Maxwell(x, a)\n\
Gammad(x, k, theta)\n\
Chisq(x, k)'

def integral(equ, x0, x1, p):
	"""
Integation of 'equ'(a string) from x0 to x1 with a partition of p points
	"""
	x0 = float(x0)
	x1 = float(x1)
	p = float(p)
	return sum(Easy_Plot(equ, x=[x0+1./(2*p), x1-1./(2*p)], partition=p)[3][0])*(x1-x0)/p

def plummer(r, d0, r0):
	return d0/(1.+(r/r0)**2)**2

def r0(M, d0, db):
	a = sqrt(sqrt(d0/db)-1)
	return (M/(2*pi*d0))**(1./3.)*(arctan(a)-a/sqrt(d0/db))**(-1./3.)

def rb(M, d0, db):
	a = sqrt(sqrt(d0/db)-1)
	return a*r0(M, d0, db)

def Bv(v, T):
	"""
Planck function of frequency and temperature in SI unit
	"""
	I = 2*h*v**3/(c**2*(exp(h*v/(k*T))-1))
	return I

def Bvc(v, T):
	"""
Planck function of frequency and temperature in cgs unit
	"""
	I = 2*hh*v**3/(cc**2*(exp(hh*v/(kk*T))-1))
	return I

def Bl(lamda, T):
	"""
Planck function of wave length and temperature in SI unit
	"""
	I = 2*h*c**2/(lamda**5*(exp(h*c/(lamda*k*T))-1))
	return I

def Blc(lamda, T):
	"""
Planck function of wave length and temperature in cgs unit
	"""
	I = 2*hh*cc**2/(lamda**5*(exp(hh*cc/(lamda*kk*T))-1))
	return I

def Normal(x, h, mean, sd):
	"""
Normal distribution
	"""
	mean = float(mean);sd = float(sd)
        sd = sd/(2.*sqrt(2.*log(2.)))
	return h*exp(-(x-mean)**2/(2*sd**2))
	#return 1./(sd*sqrt(2*pi))*exp(-(x-mean)**2/(2*sd**2))

def Maxwell(x, a):
	"""
Maxwell-Boltzmann distribution\n\
Idea gas: a=sqrt(kT/m)
	"""
	return sqrt(2/pi)*x**2*exp(-x**2/(2.*a**2))/a**3

def Gammad(x, k, theta):
	"""
Gamma distribution
	"""
	k = float(k);theta = float(theta)
	return x**(k-1)*exp(-x/theta)/(special.gamma(k)*theta**k)

def Chisq(x, k):
	"""
Chi-squared distribution
	"""
	k = float(k)
	return 1./(2**(k/2)*special.gamma(k/2))*x**(k/2-1)*exp(-x/2)

def Poisson(x, lamda):
	"""
Poisson
	"""
	if x%1 > 0.5:
		x = int(x+1.)
	else:
		x = int(x)
	return lamda**x*exp(-lamda)/math.factorial(x)

#------------------------------------

def Easy_Plot(equ_list, *args, **kwargs):
	"""
This function will return a list: [fig, ax, domain, func1, func2],\n\
A) fig is a matplotlib.figure.Figure object, but if user inputs their\n\
   own AxesSubplot, fig will be a float number 0.\n\
B) ax is a matplotlib.axes.AxesSubplot object.\n\
C) domain is a array of points of the domain.\n\
D) func1 is a list of arrays of function values of equ_list.\n\
   But the values beyond the range (when [a]semilogy/loglog, and [b]polar)\n\
   will be modified, and the original one will be stored in \n\
   the corresponding array in list func2.\n\
   In case [a], func2 contains abs(negative function value).\n\
   In case [b], func2 contains negative function values.\n\
   Otherwise the corresponding items in func2 will be a string, 'useless'.\n\n\
Some other modules also be imported:\n\
   from numpy import *\n\
   from scipy import special\n\
   import matplotlib.pyplot as plt\n\
   from matplotlib import rc\n\
   from getopt import getopt\n\
   from sys import argv\n\n\
This script also provides a constant values list:  print(const)\n\
And my function list:  print(func)\n\n\
Here add some variables, and Latex strings are allowed:\n\n\
[Necessary]\n\
   equ_list           A list, which can contain strings are functions of x,\n\
                         or/and arrays(numpy) that contains function values.\n\
                         And you can only give a string or an array,\n\
                         if just one function.\n\n\
                         If equ_list just are some arrays, you can directly\n\
                         see the shape without giving domain and range, but\n\
                         the dimension of them must be the same.\n\n\
[args]\n\
   equ_label          A list of label of functions, e.g.['$x$'];\n\
                      otherwise they'll be ordinal numbers.\n\n\
   style              Cartesian coordinate:\n\
                         0: linear plot(default), 1: semilogx,\n\
                         2: semilogy, 3: loglog.\n\
                      Polar coordinate:\n\
                         4: counterclockwise,\n\
                         6: clockwise.\n\
                      If some part is beyond the range, it'll be marked\n\
                      by another label automatically.\n\n\
   'debug'            Print the debug info.\n\n\
[kwargs]\n\
   x/domain = [a, b]     Let the domain be [a,b]. If log-scale, the base will be 10.\n\
   x/domain = [a, b, c]  Let the domain be [a,b], and the base of log-scale be c.\n\
   x/domain = an array(numpy) of data points. In this case, the x-scaling is auto.\n\
   x/domain = [(an array), a, b]\n\
                         Similar to the above, but limit the domain to be [a, b].\n\
                         If log-scale, the base will be 10.\n\
   x/domain = [(an array), c]\n\
                         Similarly, but just assign the base of log-scale.\n\
   x/domain = [(an array), a, b, c]\n\
                         Similarly, but assign the limitation and the base.\n\n\
   y/range = [a, b], y/range = [a, b, c]    All similar to the above.\n\
   y/range = 'auto'      The y-scaling is auto. But if equ_list just are some arrays,\n\
                         the y-scaling is already auto.\n\
   y/range = ['auto', c] The y-scaling is auto, and assign the base of log-scale.\n\n\
   offset = t\n\
        Let the offset for the origin be t radians, the default is 0. Operators are allowed.\n\n\
   p/partition = N\n\
        Let the size of the parition of the domain to be N.\n\
        If not, the default is 500.\n\
        Operators are allowed.\n\n\
   ex/exclusion = a number, or a list of numbers.\n\
        A list of points to be excluded.\n\
	Operators are allowed.\n\n\
   i/inf = n\n\
        If don't want to get a warning when the value is 'inf' or '-inf',\n\
        then just let them to be n.\n\
	Operators are allowed.\n\n\
   ax = a matplotlib.axes.AxesSubplot object\n\
        By this, user can create a figure with a set of subplots, then use\n\
        Easy_Plot() to plot one of them.\n\n\
   ls = a string if one plot, or a list of strings if more than one plot.\n\
        Assign the style, or/and color of lines in python form,\n\
        like 'ro' for one plot or ['ro', 'k-'] for two plots.\n\
        If doesn't assign or assign '', means set by python.\n\
        But the automatic colors for different parts won't work.\n\n\
   color = a string if one plot, or a list of strings if more than one plot.\n\
        Besides 'red' and 'blue' etc., HTML color codes are also fine.\n\
        But the user is suggested not to use ls and color to assign the color\n\
        of lines at the same time.\n\n\
   alpha = a float number between 0. and 1. if one plot,\n\
        or a list of float numbers if more than one plot.\n\
   loc = an integer.\n\
	Set the location of legend.\n\n\
   title = '...', xlabel = '...', ylabel = '...'\n\
        Assign the title, xlabel, and ylabel. But the labels is only\n\
        for the Cartesian.\n\n\
The type of variables and their default values:\n\
	====================================================\n\
	Option     Accepted assignment     Default  \n\
	====================================================\n\
	equ_list   list of string(s)/array(s)  [NECESSARY!!]\n\
		   string/array\n\
	----------------------------------------------------\n\
	style      int                     0\n\
	----------------------------------------------------\n\
	equ_label  list of string(s)       Ordinal number\n\
		   string\n\
	----------------------------------------------------\n\
	ls         list of string(s)       python default\n\
		   string\n\
	----------------------------------------------------\n\
	color      list of string(s)       python default\n\
		   string\n\
	----------------------------------------------------\n\
	alpha      list of number(s)       python default\n\
		   number\n\
	----------------------------------------------------\n\
	loc	   int                     0\n\
	----------------------------------------------------\n\
	domain     bracket of numbers -*   [-10,10]\n\
		   array\n\
	----------------------------------------------------\n\
	range      bracket of numbers -*   [-10,10],\n\
		   'auto'(a string)           base of log-scale = 10\n\
	----------------------------------------------------\n\
	offset     number -*               0\n\
	----------------------------------------------------\n\
	partition  int -*                  500\n\
	----------------------------------------------------\n\
	exclusion  list of number(s) -*    (empty)\n\
	----------------------------------------------------\n\
	inf        number -*               (empty)\n\
	----------------------------------------------------\n\
        title\n\
        xlabel     string                  (empty)\n\
        ylabel\n\
	====================================================\n\
	Suffix -*: Given with mathematical operators are allowed\n\n\
e.g. result  = Easy_Plot(['100*(10**4/(10**8+x**2)**0.5)**5'],\n\
	               ['$Plummer\\\'s\\ model\\ with\\ 5$'], 3,\n\
        	       x = [10**-2, 10**10], y = [10**-4, 10**3],\n\
	               title = '$An\\ empirical\\ model$',\n\
        	       xlabel = '$Radius$', ylabel='$Density$')\n\n\
     result2 = Easy_Plot('x/2/pi', '$spiral$', 6, x = [0, 4*pi], y = [0, 2],\n\
	               ls = 'r', offset = 1.5*pi, title = '$Sample$')\n\n\
     result3 = Easy_Plot(['abs(x)', array([2.1, 2.5, 1.2, 3.5])], ['test', 'data'], 4,\n\
	               domain = array([-0.5*pi,pi,0.5*pi,0.75*pi]), y2 = 'auto',\n\
  	               ls = ['k-', 'ro'])\n\n\
     result4 = Easy_Plot(['1/x', 'sin(1/x)'], ['$\\\\frac{1}{x}$', 'oscillator'], 2,\n\
	               x = [-10, 10], y = [0.1, 1.3], ls=['k-', ''])\n
	"""

	#[1] Default values
	style = 0; label_list = []
	domain = []; domain_i = -10.0; domain_j = 10.0; range_i = -10.0; range_j = 10.0
        b_d = 10.0; b_r = 10.0; os = 0.; p_domain = 500; exclusion = []; INF = nan
	title = None; xl = None; yl = None; pattern = []; alpha = []; col = []; location = 0
	fig = 0.; ax = None; debug = False; PLOT = True; single_ls = dict()
	#====================================================================================#

	#[2-1] Get input from user.
	#Deal with args at first, because of detecting if 'debug' mode or not.
	for item in args:
		if isinstance(item, int):
			style = item
		if isinstance(item, list):
			if len(item) != len(equ_list):
				print('Easy_Plot: [Error] The length of label list is wrong!')
				return 0
			label_list = item
		if isinstance(item, basestring):
			if item == 'debug':
				debug = True
			elif item == 'nonplot':
				PLOT = False
			else:
				label_list = [item]
		if isinstance(item, dict):
			single_ls = item

	if debug: print('Easy_Plot: ---Starting the plotting---')
	#====================================================================================#

	#[2-2] Get the main part of input
	if isinstance(equ_list, list):
		equ_string = 0
		old_range = nan
		for item in equ_list:
			if isinstance(item, ndarray):
				p_range = len(item)
				if (not isnan(old_range)) and (old_range-p_range != 0):
					print('Easy_Plot: [Error] Please input right thing!')
					return 0
				old_range = len(item)
			elif isinstance(item, float):
				print('Easy_Plot: [Error] Please input right thing! \
				Functional values should be stored in an array.')
				return 0
			else:
				equ_string += 1
		if equ_string == 0:
			domain = arange(p_range); domain_i = nan; range_i = nan; p_domain = p_range
			if debug: print('Easy_Plot: [INIT -1] auto-generating domain if do not assign.')
	elif isinstance(equ_list, basestring) or isinstance(equ_list, ndarray):
		if isinstance(equ_list, ndarray):
			p_range = len(equ_list)
			domain = arange(1,p_range+1); domain_i = nan; range_i = nan; p_domain = p_range
			if debug: print('Easy_Plot: [INIT -1] auto-generating domain if do not assign.')
		equ_list = [equ_list]
	else:
		print('Easy_Plot: [Error] Please input a list of functions (string) or functional values (array)!')
		return 0

	if debug: print('Easy_Plot: [INIT 0-1] equ_list={0}\nEasy_Plot: [INIT 0-2] args={1}\n\
Easy_Plot: [INIT 0-3] kwargs={2}'.format(equ_list, args, kwargs))
	#====================================================================================#

	#[2-3] Get the kwargs part of input
	for item in kwargs:
		if item == 'domain' or item == 'x':
			if isinstance(kwargs[item], list):
				if not isinstance(kwargs[item][0], ndarray):
					domain_i = float(kwargs[item][0])
					domain_j = float(kwargs[item][1])
					if debug: print('Easy_Plot: [INIT 1] domain=[{0},{1}]'\
							.format(domain_i, domain_j))
					if len(kwargs[item]) == 2:
						pass
					elif len(kwargs[item]) == 3:
						b_d = float(kwargs[item][2])
						if debug: print('Easy_Plot: [INIT 1-2] base in domain={0}'\
								.format(b_d))
					else:
						print('Easy_Plot: [Error] The domain is something wrong!')
						return 0
				else:
					domain = kwargs[item][0]
					p_domain = len(domain)
					if debug: print('Easy_Plot: [INIT 1] user assigns the domain data.')
					if len(kwargs[item]) == 1:
						domain_i = nan
						if debug: print('Easy_Plot: [INIT 1-2] auto scaling domain.')
					elif len(kwargs[item]) == 2:
						domain_i = nan
						b_d = float(kwargs[item][1])
						if debug: print('Easy_Plot: [INIT 1-2] auto scaling domain, base={0}'\
								.format(b_d))
					elif len(kwargs[item]) == 3:
						domain_i = float(kwargs[item][1])
						domain_j = float(kwargs[item][2])
						if debug: print('Easy_Plot: [INIT 1-2] domain=[{0},{1}]'\
								.format(domain_i, domain_j))
					elif len(kwargs[item]) == 4:
						domain_i = float(kwargs[item][1])
						domain_j = float(kwargs[item][2])
						b_d = float(kwargs[item][3])
						if debug: print('Easy_Plot: [INIT 1-2] domain=[{0},{1}], base={2}'\
								.format(domain_i, domain_j, b_d))
					else:
						print('Easy_Plot: [Error] The domain is something wrong!')
						return 0
			elif isinstance(kwargs[item], ndarray):
				domain = kwargs[item]
				p_domain = len(domain)
				domain_i = nan
				if debug: print('Easy_Plot: [INIT 1] user assigns the domain data, and auto scaling.')
			else:
				print('Easy_Plot: [Error] The domain is something wrong!')
				return 0
		elif item == 'range' or item == 'y':
			if isinstance(kwargs[item], list):
				if kwargs[item][0] != 'auto':
					range_i = float(kwargs[item][0])
					range_j = float(kwargs[item][1])
					if debug: print('Easy_Plot: [INIT 2] range=[{0},{1}]'\
							.format(range_i, range_j))
					if len(kwargs[item]) == 2:
						pass
					elif len(kwargs[item]) == 3:
						b_r = float(kwargs[item][2])
						if debug: print('Easy_Plot: [INIT 2-2] base in range={0}'\
								.format(b_r))
					else:
						print('Easy_Plot: [Error] The range is something wrong!')
						return 0
				else:
					range_i = nan
					b_r = float(kwargs[item][1])
					if debug: print('Easy_Plot: [INIT 2] auto scaling range, base={0}'\
							.format(b_r))
			elif kwargs[item] == 'auto':
				range_i = nan
				if debug: print('Easy_Plot: [INIT 2] auto scaling range.')
			else:
				print('Easy_Plot: [Error] The range is something wrong!')
				return 0
		elif item == 'p' or item == 'partition':
			p_domain = int(kwargs[item])
		elif item == 'ex' or item == 'exclusion':
			exclusion = kwargs[item]
			if not isinstance(exclusion, list): exclusion = [exclusion]
		elif item == 'offset':
			os = float(kwargs[item])
		elif item == 'title':
			title = str(kwargs[item])
		elif item == 'xlabel':
			xl = str(kwargs[item])
		elif item == 'ylabel':
			yl = str(kwargs[item])
		elif item == 'inf' or item == 'i':
			INF = float(kwargs[item])
		elif item == 'ls':
			if isinstance(kwargs[item], basestring):
				pattern = [kwargs[item]]
			elif isinstance(kwargs[item], list):
				pattern = kwargs[item]
			else:
				print('Easy_Plot: [Error] Wrong assignment for ls!')
				return 0
		elif item == 'color':
			if isinstance(kwargs[item], basestring):
				col = [kwargs[item]]
			elif isinstance(kwargs[item], list):
				col = kwargs[item]
			else:
				print('Easy_Plot: [Error] Wrong assignment for color!')
				return 0
		elif item == 'alpha':
			if isinstance(kwargs[item], basestring):
				alpha = [kwargs[item]]
			elif isinstance(kwargs[item], list):
				alpha = kwargs[item]
			else:
				print('Easy_Plot: [Error] Wrong assignment for alpha!')
				return 0
		elif item == 'loc':
			location = int(kwargs[item])
		elif item == 'ax':
			ax = kwargs[item]
		else:
			pass
	#====================================================================================#

	#[3] Check for different coordinate.
	# 	 0      1     2    3     4        6
	#	 linear logx  logy log   counter- clockwise
	# range               V    V
	# domain lin    V,log lin  V,log lin      lin
	# ax     lin    lin   lin  lin   polar    polar
	if style == 2 or style == 3:
		if b_r <= 0 or b_r == 1:
			print('Easy_Plot: [Error] Wrong base! Base of range = {0}'.format(b_r))
			return 0
		if range_i <= 0 or range_j <= 0:
			print('Easy_Plot: [Error] Wrong range!')
			return 0
	if style%2 == 1:
		if b_d <= 0 or b_d == 1:
			print('Easy_Plot: [Error] Wrong base! Base of domain = {0}'.format(b_d))
			return 0
		if domain_i <= 0 or domain_j <= 0:
			print('Easy_Plot: [Error] Wrong domain!')
			return 0
		if domain==[]:
			domain = logspace(log(domain_i)/log(b_d), log(domain_j)/log(b_d), p_domain, base = b_d)
	else:
		if domain==[]:
			domain = linspace(domain_i, domain_j, p_domain)

	if ax == None and PLOT:
		if style <= 3:
			fig, ax = plt.subplots()
		else:
			fig, ax = plt.subplots(subplot_kw = dict(polar = True))
	#====================================================================================#

	#[4] Delete exclusive points from domain.
	#If the points aren't exactly in domain but in between, then insert 'nan'
	if exclusion:
		domain = list(domain)
		k = 0
		for ex in exclusion:
			if isinstance(ex, basestring): ex = float(eval(ex))
			if ex in domain:
				domain[domain.index(ex)] = nan
			elif domain_i < ex < domain_j:
				for j in xrange(len(domain)):
					if domain[j] > ex: domain.insert(j, nan); k+=1; break
			elif domain_j < ex < domain_i:
				for j in xrange(len(domain)):
					if domain[j] < ex: domain.insert(j, nan); k+=1; break
			else:
				pass
		domain = array(domain)
		p_domain += k
	#====================================================================================#

	#[5] Generate the range in looping all functions, or just specified by user's arrays.
	#In 2nd case, we will skip then plot directly.
	#	   linear logx  logy log   counter- clockwise
	# 	   0      1     2    3     4        6
	#check 0                V    V
	#check neg              V    V     V        V
	#check INF user   user  user user  user     user
	#
	#If the range are specified by user's arrays, style will change
	#          -0.5   -1.5  -2.5 -3.5  -4.5     -6.5
	func1 = [] #A list for saving arrays which are original range points
	func2 = [] #A list for saving arrays which are auto-modified range points, but useless when lin, logx
	for i in xrange(len(equ_list)):
		func_expr = 'function[k] = {0}'.format(equ_list[i])
		func_exec = compile(func_expr, '', 'exec')
		function = empty(p_domain, float)
		function2 = empty(p_domain, float)
		if isinstance(equ_list[i], ndarray):
			if p_domain != p_range:
				print('Easy_Plot: [Error] The number of elements of one of arrays is wrong!')
				return 0
			function = equ_list[i]
			style = -style-0.5
			if debug: print('Easy_Plot: [INIT 3-No.{0}] user assigns the range data.'.format(i))
		if debug: print('Easy_Plot: [style: No.{0}-1]={1}'.format(i, style))
		# for quick colour and linestyle
		if not pattern:
			subpattern = ''
		elif len(pattern) == 1:
			subpattern = pattern[0]
		else:
			subpattern = pattern[i]
		# for alpha
		if not alpha:
			subalpha = 1.
		elif len(alpha) == 1:
			subalpha = alpha[0]
		else:
			subalpha = alpha[i]
		# for 16-colour and dict form input
		if not col:
			subcol = dict()
		elif len(col) == 1:
			subcol = dict(color=col[0])
		else:
			subcol = dict(color=col[i])
		if not isnan(INF):
			function3 = empty(p_domain, float)
			seterr(over='ignore')

		#N(partition) = N(normal points)+N(domain:nan)+N(range:neg)+N(range:inf)
		#N(range as function1) = N(plotable)+N(range:inf)
		k = 0
		n_valid_domain = p_domain; n_neg = 0; n_inf = 0
		n_valid_range = p_domain
		if style == 0 or style == 1:
			for x in domain:
				exec(func_expr)
				if isnan(x): n_valid_domain -= 1
				if not isnan(INF):
					if isinf(function[k]):
						function3[k] = INF
						function[k] = nan
						n_inf += 1
						n_valid_domain -= 1
					else:
						function3[k] = nan
				if isnan(function[k]):
					n_valid_range -= 1
				k += 1
		if abs(style) >= 2:
			for x in domain:
				if style > 0: exec(func_exec)
				if not isnan(INF):
					if isinf(function[k]):
						function3[k] = INF
						function[k] = nan
						n_inf += 1
						n_valid_domain -= 1
					else:
						function3[k] = nan
				if isnan(x):
					n_valid_domain -= 1
				elif function[k] < 0:
					function2[k] = -function[k]
					function[k] = nan
					n_neg += 1
					n_valid_domain -= 1
				elif function[k] == 0 and abs(style) < 4:
					function[k] = nan
					function2[k] = nan
					n_neg += 1
					n_valid_domain -= 1
				else:
					function2[k] = nan
				if isnan(function[k]):
					n_valid_range -= 1
				k += 1
		if style < 0: style = -int(style)
		if debug: print('Easy_Plot: [style: No.{0}-2]={1}'.format(i, style))
		func1.append(function)
		if style == 0 or style == 1:
			func2.append('useless')
		else:
			func2.append(function2)

		#If there is '$' in label, then will be in latex.
		if PLOT:
			plot_check = 0
			if not label_list:
				if i == 0:
					l1 = 1;l2 = 'st'
				elif i == 1:
					l1 = 2;l2 = 'nd'
				elif i == 2:
					l1 = 3;l2 = 'rd'
				else:
					l1 = i+1;l2 = 'th'
				if n_valid_domain != 0 and n_valid_range != 0:
					tmp = dict(dict(label=r"{0}{1}".format(l1, l2),\
					alpha=subalpha).items()+single_ls.items()+subcol.items())
					ax.plot(domain, function, subpattern, **tmp)
					plot_check += 1
				if n_neg != 0:
					if style < 4:
						tmp = dict(dict(label= r"abs({0}{1})".format(l1, l2),\
						alpha=subalpha).items()+single_ls.items()+subcol.items())
						ax.plot(domain, function2, subpattern, **tmp)
					else:
						function2 = -function2
						tmp = dict(dict(label= r"neg({0}{1})".format(l1, l2),\
						alpha=subalpha).items()+single_ls.items()+subcol.items())
						ax.plot(domain, function2, subpattern, **tmp)
					plot_check += 1
				if n_inf != 0:
					tmp = dict(dict(label= r"+-inf({0}{1})".format(l1, l2),\
					alpha=subalpha).items()+single_ls.items()+subcol.items())
					ax.plot(domain, function3, subpattern, **tmp)
					plot_check += 1
			else:
				if n_valid_domain != 0 and n_valid_range != 0:
					tmp = dict(dict(label=r"{0}".format(label_list[i]),\
					alpha=subalpha).items()+single_ls.items()+subcol.items())
					ax.plot(domain, function, subpattern, **tmp)
					plot_check += 1
				if n_neg != 0:
					if '$' in label_list[i]:
						if style < 4:
							tmp = dict(dict(label=r"$\vert${0}$\vert$"\
							.format(label_list[i]),\
							alpha=subalpha).items()+single_ls.items()+subcol.items())
							ax.plot(domain, function2, subpattern, **tmp)
						else:
							function2 = -function2
							tmp = dict(dict(label=r"$-${0}"\
							.format(label_list[i]),\
							alpha=subalpha).items()+single_ls.items()+subcol.items())
							ax.plot(domain, function2, subpattern, **tmp)
					else:
						if style < 4:
							tmp = dict(dict(label= r"abs({0})"\
							.format(label_list[i]),\
							alpha=subalpha).items()+single_ls.items()+subcol.items())
							ax.plot(domain, function2, subpattern, **tmp)
						else:
							function2 = -function2
							tmp = dict(dict(label=r"neg({0})"\
							.format(label_list[i]),\
							alpha=subalpha).items()+single_ls.items()+subcol.items())
							ax.plot(domain, function2, subpattern, **tmp)
					plot_check += 1
				if n_inf != 0:
					if '$' in label_list[i]:
						tmp = dict(dict(label= r"inf({0})".format(label_list[i]),\
						alpha=subalpha).items()+single_ls.items()+subcol.items())
						ax.plot(domain, function3, subpattern, **tmp)
					else:
						tmp = dict(dict(label= r"$\pm inf(${0}$)$".format(label_list[i]),\
						alpha=subalpha).items()+single_ls.items()+subcol.items())
						ax.plot(domain, function3, subpattern, **tmp)
					plot_check += 1
			if plot_check == 0: print('Easy_Plot: Nothing to do with func. No.{0} in equ_list!'.format(i+1))
		if debug: print('Easy_Plot: fig = {0}\n\
Easy_Plot: >>Finished plotting the func. No.{1}'.format(fig, i+1))
	if not isnan(INF): seterr(over='warn')
	#====================================================================================#

	#[6] Set x and y-axis, title, legend.
	if PLOT and plot_check != 0:
		if not isnan(domain_i) and style < 4:
			ax.set_xlim(domain_i, domain_j)
		else:
			ax.set_xlim(auto=True)
		if not isnan(range_i):
			ax.set_ylim(range_i, range_j)
		else:
			ax.set_ylim(auto=True)

		if style == 0:
			pass
		elif style == 1:
			ax.set_xscale('log', basex = b_d)
		elif style == 2:
			ax.set_yscale('log', basey = b_r)
		elif style == 3:
			ax.set_xscale('log', basex = b_d)
			ax.set_yscale('log', basey = b_r)
		elif style == 4:
			ax.set_theta_offset(os)
		else:
			ax.set_theta_offset(os)
			ax.set_theta_direction(-1)

		if title != None: ax.set_title(title, fontdict = dict(fontsize = 'x-large', va = 'baseline'))
		if xl != None and style < 4: ax.set_xlabel(xl)
		if yl != None and style < 4: ax.set_ylabel(yl)
		if (n_neg != 0 or n_inf != 0) or len(equ_list) != 1 or label_list: ax.legend(loc = location)
		ax.grid(True)
	#====================================================================================#
	if debug: print('Easy_Plot: ---The end of plotting---\n')
	return [fig, ax, domain, func1, func2]

def PN(input_array):
	"""
Separtely output the positive and negative parts of input array.
	"""
	input = copy.deepcopy(input_array)
	result = Easy_Plot(input, 2, 'nonplot')
	return [result[3][0], result[4][0]]

def Vel(input_array, time, **kwargs):
	"""
Specially for plotting loglog velocity vs radius diagram.
	"""
	result = PN(input_array)
	p_data = result[0]
	n_data = result[1]
	p_outward = Easy_Plot(p_data, 3 ,'$outward, {0}$'.format(time), alpha=0.5, ls='b.' ,**kwargs)
	n_inward = Easy_Plot(n_data, 3 ,'$inward, {0}$'.format(time), ax=p_outward[1], alpha=0.5, ls='g.', **kwargs)
	return n_inward

def main():
	
	#opts is a list of 2-tuple; equ_list is a list
	opts, equ_list = getopt(argv[1:], "hd:e:r:p:D:R:P:C:o:i:AN")

	#default values
	style = 0; domain_i = -10.0; domain_j = 10.0; range_i = -10.0; range_j = 10.0
	p_domain = 500; ex = []; b_d = 10.0; b_r = 10.0; os = 0.; INF = nan; A = False; N = True

	for opt, value in opts:
		if opt == '-h':
			print("easy_plot.py: A Python script to plot explicit real functions\n\
              of one variable, i.e. f:R->R.\n\
-----Using as a script in shell-----\n\
python easy_plot.py [-d] [-D] [-r] [-R] [-P] [-C] [-o] [-p] [-e]\n\
 [-i] [-A] [-N]  'equation 1 of x' 'equation 2 of x'...\n\
 (several functions [in numpy] marked by single/double quotes)\n\n\
   Options:\n\
   -d a,b    Let the domain to be [a,b];\n\
             the default is [-10,10].\n\
   -D a,b,B  Similar to '-d', but the log scale of the base B,\n\
             the default for B is 10.\n\
   -r c,d    Let the range to be [c,d]; the default is [-10,10].\n\
             If input -r'auto' then it'll be auto scaled.\n\
   -R c,d,B  Similar to '-r', but the log scale of the base B,\n\
             the default for B is 10.\n\
             If input -R'auto' then it'll be auto scaled.\n\
   -P c,d    Polar coordinate in the counterclockwise direction\n\
             with the range = [c,d].\n\
             If input -P'auto' then it'll be auto scaled.\n\
   -C c,d    Similar to '-P', but in the clockwise direction.\n\
   -o t      Let the offset for the origin be t radians,\n\
             the default is 0.\n\
   -p N      Let the size of the partition of the domain\n\
             to be N; if not, the default is 500.\n\
   -e p,q... Exclude points p, q... in the domain.\n\
   -i n      Let 'inf' and '-inf' to be n,\n\
             otherwise printing a warning.\n\
   -A        Print all function values of the last one function,\n\
             the 2nd contains abs(negative), if semilogy/loglog,\n\
             or just negative part, if polar.\n\
   -N        Don't print the figure.\n\
   -h        Print this help information.\n\n\
This script provided some constants, please look them up.\n\
e.g. python easy_plot.py -d 0,5**2 -r -5,5 -p 200*2 -e 0,2*10 \\\\\n\
                         1/x '-e**(1/x)' '5./2*cos(x/pi)'\n\n\
-----Using as a module in python-----\n\
In python, after typing:\n\n\
>>> import sys\n\
>>> sys.path.append('the path of this script')\n\
>>> from easy_plot import *\n\n\
Then you could use the function: Easy_Plot(equ_list, equ_label, style)\n\
For further infomation, please look up help(Easy_Plot)\n")
			exit()
		elif opt == '-d':
			tmp = value.split(",")
			domain_i = float(eval(tmp[0]))
			domain_j = float(eval(tmp[1]))
		elif opt == '-r':
			tmp = value.split(",") + [None]
			range_i = float(eval(tmp[0])) if (tmp[0]!='auto') else nan
			if tmp[1] != None: range_j = float(eval(tmp[1]))
		elif opt == '-p':
			p_domain = int(eval(value))
		elif opt == '-e':
			ex = value.split(",")
		elif opt == '-D':
			tmp = value.split(",") + [None]
			if tmp[2] != None: b_d = float(eval(tmp[2]))
			domain_i = float(eval(tmp[0]))
			domain_j = float(eval(tmp[1]))
			style += 1
		elif opt == '-R':
			tmp = value.split(",") + [None, None]
			range_i = float(eval(tmp[0])) if (tmp[0]!='auto') else nan
			if tmp[1] != None: range_j = float(eval(tmp[1]))
			if tmp[2] != None: b_r = float(eval(tmp[2]))
			style += 2
		elif opt == '-i':
			INF = float(eval(value))
		elif opt == '-P':
			style = 4
			tmp = value.split(",") + [None]
			range_i = float(eval(tmp[0])) if (tmp[0]!='auto') else nan
			if tmp[1] != None: range_j = float(eval(tmp[1]))
		elif opt == '-C':
			style = 6
			tmp = value.split(",") + [None]
			range_i = float(eval(tmp[0])) if (tmp[0]!='auto') else nan
			if tmp[1] != None: range_j = float(eval(tmp[1]))
		elif opt == '-o':
			os = float(eval(value))
		elif opt == '-A':
			A = True
		elif opt == '-N':
			N = False
		else:
			pass

	if not equ_list:
		print('[Error] You did NOT assign any functions!')
		exit()

	number = len(equ_list)-1
	M = Easy_Plot(equ_list, style, x = [domain_i, domain_j, b_d], y = [range_i, range_j, b_r], partition = p_domain, exclusion = ex, offset = os, inf = INF)
	if A:
		print('function={}'.format(M[3][number]))
		if style > 1: print('function2={}'.format(M[4][number]))
	if N:
		plt.show()

if __name__ == '__main__':
    main()
