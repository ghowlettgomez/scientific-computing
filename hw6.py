import numpy 
import math
import matplotlib.pyplot as plt

class continuousPieceWiseQuadratic():

	def __init__(self, mesh):
		self.mesh = mesh

	def fill(self, funct):
		self.values = []
		for item in mesh:
			self.values.append(funct(item))
		self.bumps = []
		for i in range(len(mesh) - 1):
			m_i = (mesh[i] + mesh[i+1]) / 2.0
			bump = funct(m_i) - ((self.values[i] + self.values[i+1]) / 2.0) 
			self.bumps.append(bump)

def cpwq_inner(cpwq1, cpwq2):
	return sum([i * j for i,j in zip(cpwq1.values, cpwq2.values)] + [i * j for i,j in zip(cpwq1.bumps, cpwq2.bumps)])

def cpwq_sum(cpwq1, cpwq2):
	res = continuousPieceWiseQuadratic(cpwq1.mesh)
	res.values = []
	res.bumps = []
	for value1, value2, in zip(cpwq1.values, cpwq2.values):
		res.values.append(value1 + value2)
	for bump1, bump2 in zip(cpwq1.bumps, cpwq2.bumps):
		res.bumps.append(bump1 + bump2)
	return res

def cpwq_minus(cpwq1, cpwq2):
	res = continuousPieceWiseQuadratic(cpwq1.mesh)
	res.values = []
	res.bumps = []
	for value1, value2, in zip(cpwq1.values, cpwq2.values):
		res.values.append(value1 - value2)
	for bump1, bump2 in zip(cpwq1.bumps, cpwq2.bumps):
		res.bumps.append(bump1 - bump2)
	return res

def cpwq_times(scalar, cpwq):
	res = continuousPieceWiseQuadratic(cpwq.mesh)
	res.values = []
	res.bumps = []
	for i in range(len(cpwq.values)):
		res.values.append(cpwq.values[i] * scalar)
	for i in range(len(cpwq.bumps)):
		res.bumps.append(cpwq.bumps[i] * scalar)
	return res

def mass_matrix_multiply(cpwq, funct):
	res = continuousPieceWiseQuadratic(cpwq.mesh)
	res.values = [0.0 for i in range(len(cpwq.values))]
	res.bumps = []
	for i in range(len(cpwq.values) - 1):
		dx = cpwq.mesh[i+1] - cpwq.mesh[i]
		m = (cpwq.mesh[i+1] + cpwq.mesh[i]) / 2.0
		Zm = cpwq.bumps[i] + 0.5 * (cpwq.values[i]+cpwq.values[i+1])
		value1 = (dx / 6.0) * ((cpwq.values[i] * funct(cpwq.mesh[i])) + 2.0 * (funct(m) * Zm))
		bump = (dx / 6.0) * (4.0 * Zm * funct(m))
		value2 = (dx / 6.0) * (2.0 * (Zm * funct(m)) + (cpwq.values[i+1] * funct(cpwq.mesh[i+1])))
		res.values[i] += value1
		res.values[i+1] += value2
		res.bumps.append(bump)
	return res 

def stiffness_matrix_multiply(cpwq, funct):
	res = continuousPieceWiseQuadratic(cpwq.mesh)
	res.values = [0.0 for i in range(len(cpwq.values))]
	res.bumps = []
	for i in range(len(cpwq.values) - 1):
		dx = cpwq.mesh[i+1] - cpwq.mesh[i]
		m = (cpwq.mesh[i+1] + cpwq.mesh[i]) / 2
		zpl = (cpwq.values[i+1]-cpwq.values[i])/dx + ((4/dx)  * cpwq.bumps[i])
		zpm = (cpwq.values[i+1]-cpwq.values[i])/dx
		zpr = (cpwq.values[i+1]-cpwq.values[i])/dx - ((4/dx) * cpwq.bumps[i])
		value1 = (dx/6) * ( zpl * funct(cpwq.mesh[i]) * (-1/dx) + 4 * zpm * funct(m) * (-1/dx) + zpr * funct(cpwq.mesh[i+1]) * (-1/dx) )
		bump = (dx/6) * ( zpl * funct(cpwq.mesh[i]) * (4/dx) + 0 + zpr * funct(cpwq.mesh[i+1]) * (-4/dx))
		res.bumps.append(bump)
		res.values[i] += value1
		res.values[i+1] -= value1
	return res

class Neumann_Problem():

	def __init__(self, mesh,cond, lrate,dflux,lflux,rflux):
		self.mesh = mesh
		self.cond = cond 
		self.lrate = lrate
		self.dflux = dflux
		self.lflux = lflux
		self.rflux = rflux

def multiply_matrices(cpwq, funct):
	S = stiffness_matrix_multiply(cpwq, funct)
	M = mass_matrix_multiply(cpwq, funct)
	return cpwq_sum(S, M)

def conjugate_gradient(nprob):
	df = continuousPieceWiseQuadratic(nprob.mesh)
	df.fill(nprob.dflux)
	x = df
	g = mass_matrix_multiply(df, nprob.lrate)
	g.values[0] += nprob.lflux
	g.values[-1] -= nprob.rflux
	Ax = multiply_matrices(df, nprob.dflux)
	r = cpwq_minus(g, Ax)
	p = r
	ks = []
	normseq = []
	for i in range(23):
		Apk = multiply_matrices(p, nprob.dflux)
		try:
			a = cpwq_inner(r, r) / cpwq_inner(p, Apk)
		except:
			break
		nx = cpwq_sum(x, cpwq_times(a, p))
		nr = cpwq_minus(r, cpwq_times(a, Apk))
		normseq.append(cpwq_inner(nr, nr))
		ks.append(i)
		if normseq[-1] <= 1e-20:
			return nx, normseq, ks
		beta = cpwq_inner(nr,nr) / cpwq_inner(r,r)
		p = cpwq_sum(nr, cpwq_times(beta, p))
		r = nr 
		x = nx
	return x, normseq, ks

def test_funct(real):
	return real

def test_weight(real):
	return real ** 2

def test_weight2(real):
	return real ** 3

def cond1(real):
	return 1

def lrate1(real):
	return 1.0e-8

def dflux1(real):
	return 0

def cond2(real):
	return 1

def lrate2(real):
	return 1

def dflux2(real):
	return 1

def cond3(real):
	if real < 0.5:
		return 10
	else:
		return 1

def lrate3(real):
	return 1.0e-8

def dflux3(real):
	return 0

def plot_problem(nprob):
	res1, norms1, ks1 = conjugate_gradient(nprob)
	plt.title('x vs y')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.plot(res1.mesh, res1.values)
	plt.show()
	plt.title('iteration vs norm')
	plt.xlabel('iteration')
	plt.ylabel('norm')
	plt.plot(ks1, norms1)
	plt.show()

mesh = [0.0, 0.1, 0.2, 0.3, 0.4, 0.49, 0.5, 0.625, 0.7, 0.8, 0.9, 1.0]
#nprob = Neumann_Problem(mesh, test_weight, test_weight2, test_funct, 5.0, 6.0)
#res = conjugate_gradient(nprob)
nprob1 = Neumann_Problem(mesh, cond1,lrate1,dflux1, -1.0, -1.0)
nprob2 = Neumann_Problem(mesh, cond2,lrate2,dflux2, 0.0, -math.e + 1.0/math.e)
nprob3 = Neumann_Problem(mesh, cond3,lrate3,dflux3, -1.0, -1.0)
plot_problem(nprob1)
plot_problem(nprob2)
plot_problem(nprob3)





