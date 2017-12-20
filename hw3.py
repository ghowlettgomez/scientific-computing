import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy, pylab

def array_multiply(array, scalar):
	return [i * scalar for i in array]

def array_add(array1, array2):
	return [i + j for i, j in zip(array1, array2)]

def array_subtract(array1, array2):
	return [i - j for i,j in zip(array1, array2)]

def norm(array):
	res = 0
	for i in array:
		res += i ** 2.0
	return math.sqrt(res)

def step_function(f, Ynow, t_now, dt):
	#takes in f, ynow, tnow and dt and return the new y, new f, and error indicator
	f_now = f(t_now, Ynow)
	k1 = f(t_now,Ynow)
	k2 = f(t_now+(dt/2), array_add(Ynow, array_multiply(k1, dt/2)))
	k3 = f(t_now+(dt/2), array_add(Ynow, array_multiply(k2, dt/2)))
	k4 = f(t_now+dt, array_add(Ynow, array_multiply(k3, dt)))
	dY = array_multiply(array_add(array_add(array_add(k1, array_multiply(k2, 2)), array_multiply(k3, 2)), k4), dt/6.0)
	Ynew = array_add(Ynow, dY)
	f_new = f(t_now + dt, Ynew)
	ei = error_indicator(f, Ynow, Ynew, f_now, f_new, t_now, dt)
	return Ynew, f_new, ei

def hermite_cubic_solve(dt, Ynow, Ynew, f_now, f_new, t):
	#solves the hermite cubic function
	h_00 = (2.0 * (t ** 3.0)) - (3.0 * (t ** 2.0)) + 1.0
	h_10 = (t ** 3.0) - (2.0 * (t ** 2.0)) + t 
	h_01 = (-2.0 * (t ** 3.0)) + (3.0 * (t ** 2.0))
	h_11 = (t ** 3.0) - (t ** 2.0)
	return array_add(array_multiply(Ynow, h_00), array_add(array_multiply(f_now, h_01 * dt), array_add(array_multiply(Ynew, h_01), array_multiply(f_new, h_11 * dt))))

def error_indicator(f, Ynow, Ynew, f_now, f_new, t_now, dt): 
	t_new = t_now + dt
	c1 = 0.4 - math.sqrt(0.06)
	c2 = 0.4 + math.sqrt(0.06)
	t_c1 = t_now + (dt * c1)
	t_c2 = t_now + (dt * c2) 
	Y_c1 = hermite_cubic_solve(dt, Ynow, Ynew, f_now, f_new, t_c1)
	Y_c2 = hermite_cubic_solve(dt, Ynow, Ynew, f_now, f_new, t_c2)
	f_c1 = f(t_c1, Y_c1)
	f_c2 = f(t_c2, Y_c2)
	w1 = (3 * c2 -1) / (6 * (c1 - c2) * (c1-1))
	w2 = (3 * c1 -1) / (6 * (c2 - c1) * (c2-1))
	wnew = 1 - w1 - w2
	dY_alt = array_multiply(array_multiply(array_add(f_c1, array_multiply(array_add(f_c2, array_multiply(f_new, wnew)), w2)), w1), dt)
	dY = array_subtract(Ynew, Ynow)
	return norm(array_subtract(dY, dY_alt))

class simTime():
	def __init__(self, time, dt, tol, agrow, ashrink, dtmax, dtmin, endTime):
		self.time = float(time) 
		self.dt = dtmin 
		self.tol = float(tol)
		self.agrow = float(agrow) 
		self.ashrink = float(ashrink)
		self.dtmin = float(dtmin) 
		self.dtmax = float(dtmax) 
		self.endTime = float(endTime)
		self.stepsSinceRejection = 0
		self.stepsRejected = 0
		self.stepsAccepted = 0

	def advance(self, current_sol, f):
		resfile = open('temp.txt', 'w')
		resfile.write('time, dt, log10dt, ei, Y_now\n')
		res = current_sol
		x_axis = []
		y_axis = []
		while self.time < self.endTime:
			x_axis.append(self.time)
			y_axis.append(norm(res[2:4]))
			while True:
				step = step_function(f, res, self.time, self.dt)
				resfile.write((','.join([str(self.time), str(self.dt), str(math.log(self.dt)), str(step[2]), str(res)])) + '\n')
				if step[2] > self.tol and self.dt > self.dtmin:
					self.dt /= 2.0
					self.stepsSinceRejection = 0
					self.stepsRejected += 1
				else:
					res = step[0]
					self.stepsSinceRejection += 1
					self.stepsAccepted += 1
					break
			if step[2] < self.tol / 4.0 and self.stepsSinceRejection >= 5:
				self.dt *= self.agrow
			elif step[2] < 3 * self.tol / 4.0:
				self.dt *= self.ashrink
			if self.dt < self.dtmin:
				self.dt = self.dtmin
			if self.dt > self.dtmax:
				self.dt = self.dtmax
			if self.time + self.dt > self.endTime:
				self.dt /= 2.0
				if self.time + self.dt > self.endTime:
					self.dt = self.endTime - self.time
			self.time += self.dt
		plt.title('Comet time vs velocity')
		print x_axis 
		print y_axis
		plt.xlabel('time')
		plt.xlabel('velocity')
		plt.plot(x_axis, y_axis)
		plt.show()
		#return res

def test_func(scalar, array): 
	return [scalar * i for i in array]

def apow(scalar, array):
	return [scalar ** i for i in array]

def comet_function(scalar, array):
	#the first two items are the position, the final two are the velocity, both are in x,y form
	print array
	res = [i for i in array]
	factor = norm(array[0:2]) ** 3.0
	res[2] = (-array[0] / factor)
	res[3] = (-array[1] / factor)
	res[0:2] = array[2:4]
	return res

def energy_func(scalar, array):
	return (0.5 * math.sqrt(norm(array[0:2])) - 1) / norm(array[0:2])

input_arr = [1.0, 0.0, 0.0, 0.3]
#input_arr = [2.0]
t_now = 0.0
dt = 0.1
length_time = 3.0 * 2.0 * math.pi / ((2.0 - math.sqrt(0.3)) ** 1.5)
current_sol = test_func(t_now, input_arr)
#print step_function(test_func, input_arr, t_now, dt)
sim = simTime(t_now, .05, .1, 1.2, .8, .1,  .001, length_time)
print sim.advance(input_arr,comet_function)


