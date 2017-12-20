import numpy as np
import numpy.linalg as linalg
import math
from random import randint

class QR_Decomposition():
	def __init__(self):
		self.qs = []

	def compute_householder(self, x):
		#computes the householder matrix from self.matrix and the given x
		norm = linalg.norm(x)
		lenx = len(x)
		ae_0 = [0 for i in range(lenx)]
		ae_0[0] = 1.0
		e_0 = np.asarray(ae_0)
		v = x - (norm * e_0)
		if v[0] == 0:
			beta = 0
		else:
			beta = 2.0 / (np.dot(v, v))
		P = np.subtract(np.identity(x.size), (np.outer(v, v) * beta))
		if self.lenmatrix <= lenx:
			return P 
		else:
			topleft = np.identity(self.lenmatrix - lenx)
			topright = np.zeros((self.lenmatrix - lenx, lenx))
			bottomleft = np.zeros((lenx, self.lenmatrix - lenx))
			return np.bmat([[topleft, topright], [bottomleft, P]])

	def solve_QR(self, matrix, y):
		self.compute_RQ(matrix)
		qy = np.dot(self.Q, y.T)
		return self.upper_triangular_solve(qy)

	def compute_RQ(self, matrix):
		#computes both r and q
		self.matrix = matrix
		self.lenmatrix = len(matrix)
		self.qs = []
		self.compute_R()
		self.compute_Q()

	def compute_R(self):
		#computes r
		for i in range(self.lenmatrix):
			x = []
			for j in range(i, self.lenmatrix):
				x.append(self.matrix.item(j,i))
			x = np.array(x)
			q = self.compute_householder(x)
			self.matrix = np.dot(q, self.matrix)
			self.qs.append(q)

	def compute_Q(self):
		#computes q
		Q = np.identity(self.lenmatrix)
		for q in self.qs:
			Q = np.dot(Q, q.T)
		self.Q = Q.T

	def upper_triangular_solve(self, y):
		#solves self.matrix x = y for x
		x = [y.item(i) for i in range(self.lenmatrix)]
		for i in range(1, self.lenmatrix + 1):
			for j in range(1, i):
				x[-i] -= (self.matrix.item(self.lenmatrix - i, self.lenmatrix - j) * x[-j])
			x[-i] /= self.matrix.item(self.lenmatrix - i, self.lenmatrix - i)
		return np.array(x)


def makeAs():
	#make A and A prime
	I = np.identity(5)
	I_prime = np.matrix([[(1.0 if 4 - j == i else 0.0) for j in range(5)] for i in range(5)])
	c = np.matrix([[1.0 / randint(10, 100) for j in range(5)] for i in range(5)])
	return [I + c, I_prime + c]

A, A1 = makeAs()
Z = np.array([float(i) for i in range(1, 6)])
qrd = QR_Decomposition()
B = np.dot(A, Z.T)
print qrd.solve_QR(A, B)
#returns [ 1.  2.  3.  4.  5.]
B1 = np.dot(A1, Z.T)
print qrd.solve_QR(A1, B1)
#returns [ 1.  2.  3.  4.  5.]
