import math

class Tri_Diag():
	def __init__(self,b,d,a):
		self.b = b
		self.d = d
		self.a = a

	def multiply(self, v):
		res = []
		lena = len(self.a)
		for i in range(len(v)):
			d_factor = v[i] * self.d[i]
			a_factor = 0
			b_factor = 0
			if i < lena:
				a_factor = v[i + 1] * self.a[i]
			if i > 0:
				b_factor = v[i-1] * self.b[i - 1] 
			res.append(d_factor + b_factor + a_factor)
		return res

	def factor(self):
		for i in range(len(self.a)):
			self.a[i] /= self.d[i]
			self.d[i+1] -= (self.a[i] * self.b[i])

	def solve(self, y):
		leny = len(y)
		for i in range(leny - 1):
			y[i] /= self.d[i]
			y[i+1] -= (y[i] * self.b[i]) 
		x = [0 for i in range(leny)]
		for i in range(1, leny + 1):
			x[-i] = y[-i] / d[-i] 
			if i > 1:
				x[-i] = (y[-i] - (x[-(i-1)] * a[-(i-1)])) 
		return x


d = [2.1 + 0.1*i for i in range(0, 5)] 
b = [0.01*(i*i) for i in range(0, 4)]
a = [0.03*i for i in range(0, 4)]
x = [float(i) for i in range(1, 6)]
A = Tri_Diag(b,d,a)
z = A.multiply(x)
A.factor()
solution =  A.solve(z)
diff = sum([abs(solution[i] - x[i]) for i in range(len(x))])
print "The solution is: " + str(solution)
print "The difference is: " + str(diff)
'''
It will print 
The solution is: [1.0, 2.0, 2.9999999999999996, 3.9999999999999996, 5.0]
The difference is: 8.881784197e-16
'''





