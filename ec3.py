import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy, pylab

#Gregory Howlett-Gomez
#CMSC 28510

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
	f_now = f(0, Ynow)
	k1 = f(0,Ynow)
	k2 = f((dt/2), array_add(Ynow, array_multiply(k1, dt/2)))
	k3 = f((dt/2), array_add(Ynow, array_multiply(k2, dt/2)))
	k4 = f(dt, array_add(Ynow, array_multiply(k3, dt)))
	dY = array_multiply(array_add(array_add(array_add(k1, array_multiply(k2, 2)), array_multiply(k3, 2)), k4), dt/6.0)
	Ynew = array_add(Ynow, dY)
	f_new = f(dt, Ynew)
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
	def __init__(self, time, dt, tol, agrow, ashrink, dtmax, dtmin, endTime, p_prime):
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
		self.p_prime = p_prime

	def advance(self, current_sol, f):
		res = current_sol
		times = []
		p0xs = []
		p0ys = []
		p1xs = []
		p1ys = []
		normds = []
		while self.time < self.endTime:
			times.append(self.time)
			p0xs.append(res[0])
			p0ys.append(res[1])
			p1xs.append(res[2])
			p1ys.append(res[3])
			normds.append(norm(array_subtract(self.p_prime[2:], self.p_prime[0:2])))
			while True:
				step = step_function(f, res, self.time, self.dt)
				if step[2] > self.tol and self.dt > self.dtmin:
					self.dt /= 2.0
					self.stepsSinceRejection = 0
					self.stepsRejected += 1
				else:
					res = step[0]
					self.p_prime = step[1]
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
		plt.title('p0_x vs time')
		plt.xlabel('time')
		plt.ylabel('p0_x')
		plt.plot(times, p0xs)
		plt.show()
		plt.title('p0_y vs time')
		plt.ylabel('p0_y')
		plt.plot(times, p0ys)
		plt.show()
		plt.title('p1_x vs time')
		plt.ylabel('p1_x')
		plt.plot(times, p1xs)
		plt.show()
		plt.title('p1_y vs time')
		plt.ylabel('p1_y')
		plt.plot(times,p1ys)
		plt.show()
		plt.title('norm vs time')
		plt.ylabel('norm')
		plt.plot(times, normds)
		plt.show()
		return res

	def damped_oscillator(self, dt, array):
		d = array_subtract(array[2:], array[0:2])
		normd = norm(d)
		unitd = array_multiply(d, 1.0 / normd)
		fm = 10.0 * (1.0 - normd) - normd
		force = array_multiply(unitd, fm * dt)
		pp0 = array_subtract(self.p_prime[0:2], force)
		pp1 = array_add(self.p_prime[2:], force)
		return pp0 + pp1

	def moon(self, dt, array):
		d = array_subtract(array[2:], array[0:2])
		normd = norm(d)
		unitd = array_multiply(d, 1.0 / normd)
		fm = 10.0 * (1.0 - normd) - normd
		force = array_multiply(unitd, fm * dt)
		gravp0 = 1 / (norm(array[0:2]) ** 2)
		gravp1 = 1 / (norm(array[2:]) ** 2)
		pp0 = array_subtract(self.p_prime[0:2], force)
		pp1 = array_add(self.p_prime[2:], force)
		for i in range(len(pp0)):
			pp0[i] = pp0[i] - gravp0 if pp0[i] >= 0 else pp0[i] + gravp0
			pp1[i] = pp1[i] - gravp1 if pp1[i] >= 0 else pp1[i] + gravp1
		return pp0 + pp1

	def spinning_commet(self, dt, array):
		d = array_subtract(array[2:], array[0:2])
		normd = norm(d)
		unitd = array_multiply(d, 1.0 / normd)
		fm = 10.0 * (1.0 - normd) - normd
		force = array_multiply(unitd, fm * dt)
		gravp0 = 1 / (norm(array[0:2]) ** 2)
		gravp1 = 1 / (norm(array[2:]) ** 2)
		pp0 = array_subtract(self.p_prime[0:2], force)
		pp1 = array_add(self.p_prime[2:], force)
		for i in range(len(pp0)):
			pp0[i] = pp0[i] - gravp0 if pp0[i] >= 0 else pp0[i] + gravp0
			pp1[i] = pp1[i] - gravp1 if pp1[i] >= 0 else pp1[i] + gravp1
		return pp0 + pp1

#for damped oscillator
#p_start = [0.07, 0.0, -0.07, 0.0]
#p_prime = [0.0, 0.5, 0.0, -0.5]
#p_prime = [0.1, 1.5, 0.1, 0.5]
#sim = simTime(0.0, .05, 1e-4, 1.2, .8, .1,  .001, 10.0, p_prime)
#sim.advance(p_start,sim.damped_oscillator)

#for moon
#p_start = [1.0, .05, 1.0,-0.05]
#p_prime = [0.4, 1.0, -0.4, 1.0]
#sim = simTime(0.0, .05, 1e-4, 1.2, .8, .1,  .001, 12.0 * math.pi, p_prime)
#sim.advance(p_start, sim.moon)

#for spinning commet
#p_start = [1.0, .05, 1.0,-0.05]
#p_prime = [0.4, 0.5, -0.4, 0.5]
#sim = simTime(0.0, .05, 1e-4, 1.2, .8, .1,  .001, 6 * math.pi, p_prime)
#sim.advance(p_start, sim.spinning_commet)
