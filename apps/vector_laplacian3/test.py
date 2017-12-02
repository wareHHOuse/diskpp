import numpy as np
from numpy import linalg as LA

a = np.array([[2, 0, 0], [0, -1, 0], [0, 0, -0.1]])

print(a)

print(LA.eigvalsh(a))

print(LA.cond(a, 'fro'))

print(LA.cond(a, 2))

print(LA.cond(a, -2))
