import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy, pylab

max_depth = 0
func_calls = 0
good_interval_length = []

def recursion_function(math_funct, interval, tolerance, depth=1, coarse_sum=None):
	global func_calls
	global max_depth
	global good_interval_length
	midpoint = (interval[0] + interval[1]) / 2.0
	if coarse_sum is None:
		coarse_sum = math_funct(midpoint) * (interval[1] - interval[0]) 
		func_calls += 1
	fine_sum_left = math_funct((interval[0] + midpoint) / 2.0) * (midpoint - interval[0])
	fine_sum_right = math_funct((interval[1] + midpoint) / 2.0) * (interval[1] - midpoint)
	fine_sum = fine_sum_left + fine_sum_right
	func_calls += 2
	newdepth = depth + 1
	if depth >= 30:
		max_depth = 30
		good_interval_length.append([midpoint, -math.log(interval[1] - interval[0])])
		return fine_sum
	if depth < 4:
		return recursion_function(math_funct, [interval[0], midpoint], tolerance / 2.0, depth=newdepth, coarse_sum=fine_sum_left) + recursion_function(math_funct, [midpoint, interval[1]], tolerance / 2.0, depth=newdepth, coarse_sum=fine_sum_right)
	if abs(float(coarse_sum) - float(fine_sum)) <= 3 * tolerance:
		if depth > max_depth:
			max_depth = depth
		good_interval_length.append([midpoint, -math.log(interval[1] - interval[0])])
		return fine_sum
	else:
		return recursion_function(math_funct, [interval[0], midpoint], tolerance / 2.0, depth=newdepth, coarse_sum=fine_sum_left) + recursion_function(math_funct, [midpoint, interval[1]], tolerance / 2.0, depth=newdepth, coarse_sum=fine_sum_right)

def x_sq(real):
	return real**2

def x_qd(real):
	return real**4

def sqrt_pi(real):
	return math.sqrt(math.sin(math.pi * (real / 2.0)))

def func_4(real):
	if real < (1.0 / 3.0):
		return 1
	else:
		return 0

def safe_log(real):
	try:
		return math.log(real)
	except:
		return np.inf

def safe_div(real1, real2):
	try: 
		return 1 / abs(real1 - real2)
	except:
		return np.inf

def eval_function(funct, interval, title, plot_to_graph):
	'''evaluates the funciton and plots the log of 1 over the tolerance vs log of function evaluations in red
	evaluates and plots the log of the inverse tolerance vs log of the inverse of the error in green
	and evaluates and plots the negative log of good interval length vs x in blue 
	'''
	global func_calls
	global max_depth
	global good_interval_length
	actual_value = recursion_function(funct, interval, 10.0**-10)
	tolerances = [10.0**(-2), 3*(10.0**(-3)), 10.0**(-3), 3*(10.0**(-4)), 10.0**(-4), 3*(10.0**(-5)), 10.0**(-5), 3*(10.0**(-6)), 10.0**(-6)]
	res_list = []
	for tolerance in tolerances:
		good_interval_length = []
		max_depth = 0
		func_calls = 0
		res_list.append([recursion_function(funct, interval, tolerance), max_depth, func_calls])
	log_inv_tol = [safe_log(1.0 / x) for x in tolerances]
	log_funct_eval = [safe_log(x[2]) for x in res_list]
	log_error_array = [safe_log(safe_div(actual_value, x[0])) for x in res_list]
	x_values = [x[0] for x in good_interval_length]
	int_lengths = [x[1] for x in good_interval_length]
	plt.title(title)
	if plot_to_graph == 'log_funct_eval':
		plt.plot(log_inv_tol, log_funct_eval, 'r')
		plt.xlabel('Log of inverse tolerance')
		plt.ylabel('Log of function evaluations')
	elif plot_to_graph == 'log_error_array':
		plt.plot(log_inv_tol, log_error_array, 'g')
		plt.xlabel('Log of inverse tolerance')
		plt.ylabel('Log of the inverse of actual error')
	else:
		plt.plot(x_values, int_lengths, 'b')
		plt.xlabel('x values')
		plt.ylabel('good interval length')
	plt.show()

#eval_function(x_sq, [0, 2], 'x^2 on [0,2]', 'log_funct_eval')
#eval_function(sqrt_pi, [0,1], 'sqrt( sin( pi x/2)) on [0,1]', 'log_funct_eval')
#eval_function(x_qd, [-1, 1], 'x^4 on [-1, 1]', 'log_funct_eval')
#eval_function(func_4, [0, 1], '1 if x<1.0/3.0 and 0 if x>=1.0/3.0 on [0,1]', 'log_funct_eval')
#eval_function(x_sq, [0, 2], 'x^2 on [0,2]', 'log_error_array')
#eval_function(sqrt_pi, [0,1], 'sqrt( sin( pi x/2)) on [0,1]', 'log_error_array')
#eval_function(x_qd, [-1, 1], 'x^4 on [-1, 1]', 'log_error_array')
#eval_function(func_4, [0, 1], '1 if x<1.0/3.0 and 0 if x>=1.0/3.0 on [0,1]', 'log_error_array')
#eval_function(x_sq, [0, 2], 'x^2 on [0,2]', 'good_interval_length')
#eval_function(sqrt_pi, [0,1], 'sqrt( sin( pi x/2)) on [0,1]', 'good_interval_length')
#eval_function(x_qd, [-1, 1], 'x^4 on [-1, 1]', 'good_interval_length')
#eval_function(func_4, [0, 1], '1 if x<1.0/3.0 and 0 if x>=1.0/3.0 on [0,1]', 'good_interval_length')