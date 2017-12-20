import math
import numpy as np
import matplotlib.pyplot as plt

#Gregory Howlett-Gomez
#Project 5, CS28510

class Solver():
	def __init__(self, function, interval, init_value, norm, dt, tolerance, dtmin, dt_factor, max_res):
		self.function = function 
		self.interval = interval
		self.dt_vals = [[interval[0]],[init_value]]
		self.y = interval[0]
		self.norm = norm
		self.dt = dt
		self.tolerance = tolerance
		self.dt_min = dtmin
		self.steps_since_error = 0
		self.dt_factor = dt_factor
		self.steps_rejected = 0
		self.max_res = max_res

	def time_advance(self):
	#advances time and adds time and the y values to list
		time = self.interval[0]
		i = 0
		dy = np.array([1,1,1])
		while time < self.interval[1]:
			dmid, res1 = self.backwards_euler(self.dt / 2.0, dy, self.dt_vals[0][i], self.dt_vals[1][i], self.dt_vals[0][i-1], self.dt_vals[1][i-1])
			dnew, res2 = self.backwards_euler(self.dt / 2.0, dy, self.dt_vals[0][i] + self.dt/2.0, dmid, self.dt_vals[0][i], self.dt_vals[1][i])
			snew, res3 = self.backwards_euler(self.dt, dy, self.dt_vals[0][i], self.dt_vals[1][i], self.dt_vals[0][i-1], self.dt_vals[1][i-1])
			error = abs(self.norm(dnew - snew))
			if self.dt < self.dt_min or error < self.tolerance and all([res1 < self.max_res, res2 < self.max_res, res3 < self.max_res]):
				self.steps_since_error += 1
				time += self.dt
				ynew = 2.0 * dnew - snew
				dy = ynew - self.y 
				self.y = ynew
				self.dt_vals[0].append(time)
				self.dt_vals[1].append(ynew)
				i += 1
				if self.dt < self.dt_min:
					continue
				if error < 0.25 * self.tolerance:
					self.dt *= self.dt_factor
				if error > 0.95 * self.tolerance:
					self.dt /= self.dt_factor
			else:
				self.dt /= self.dt_factor
				self.steps_since_error = 0
				self.steps_rejected += 1
				continue
			if self.dt_vals[0][i] + self.dt > self.interval[1] and self.dt > self.dt_min:
				self.dt /= 2.0

	def backwards_euler(self, dt, dy, tnow, ynow, told, yold):
	#does backwards euler method
		r1, j1 = self.function(tnow, ynow + dy)
		residual0 = (dy / dt) - r1
		jacobian_inv = np.linalg.inv((np.identity(3) / dt) - j1)
		ddy = -np.dot(jacobian_inv, residual0)
		dy += np.array([ddy.item(i) for i in range(dy.size)])
		r2, _ = self.function(tnow, ynow + dy)
		residual1 = (dy / dt) - r1
		rounding_error = (10 ** -15) * self.norm((np.abs(r2) + (1 / dt) * np.abs(dy)))
		residualratio = self.norm(residual1) / (self.norm(residual0) + 1000 * rounding_error)
		return ynow+dy, residualratio

def bz_norm(y):
	return (abs(y.item(0)) / (1.25 * (10.0 ** 5))) + (abs(y.item(1)) / 1800.0) + (abs(y.item(2)) / 3.0 * (10 ** 4)) 

def bz_reaction(t, y):
	y0 = y.item(0)
	y1 = y.item(1)
	y2 = y.item(2)
	f_0 = 77.27 * (y1 - (y0 * y1) + y0 - (-8.375 * (10.0 ** -6.0) * (y0 ** 2.0)))
	f_1 = (1.0 / 77.27) * (-y1 - (y0 * y1) + y2)
	f_2 = .161 * (y0 - y2)
	reaction = np.array([f_0, f_1, f_2])

	y0 = y.item(0)
	y1 = y.item(1)
	y2 = y.item(2)
	f_00 = 77.27 * (1.0 - y1 - 2.0 * (8.375 * (10.0 ** -6.0)) * y0)
	f_01 = 77.27 * (1.0 - y0)
	f_02 = 0
	f_10 = (1.0 / 77.27) * (-y1)
	f_11 = (1.0 / 77.27) * (-1.0 - y0)
	f_12 = (1.0 / 77.27)
	f_20 = 0.161
	f_21 = 0
	f_22 = -0.161
	jacobian = np.matrix([[f_00, f_01, f_02], [f_10, f_11, f_12], [f_20, f_21, f_22]])
	return reaction, jacobian

init_value = np.array([4.0, 1.1, 4.0]).T
s = Solver(bz_reaction, [0.0, 700.0], init_value, bz_norm, 0.5, 0.01, 0.00001, 1.2, 1.0/100.0)
s.time_advance()
plt.title('Norm vs Time')
plt.xlabel('Time')
plt.ylabel('Norm')
for i in range(3):
	yi = [s.dt_vals[1][j][i] for j in range(len(s.dt_vals[1]))]
	plt.plot(s.dt_vals[0], yi)
	plt.show()

