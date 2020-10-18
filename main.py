import numpy as np
import math
from find_roots import *

data = []
tmp = 0

#read data file
with open('input.txt', 'rt') as inputfile:
        for inputline in inputfile:
                if tmp == 0:
                        tmp = 1
                else:
                        linedata = inputline.split()
                        data.append(float(linedata[2]))

#calculate average
mean = np.average(data, axis=0)
#calculate average of the cubed data
meanCube = np.average(np.power(data, 3), axis=0)

#calculate the cumulative
cumulative = 0
for element in data:
        if element < mean:
                cumulative = cumulative + 1

cumulative = cumulative / len(data)

def f(x):
   return cumulative + math.exp(-(mean / ((meanCube / math.gamma(1 + 3 / x)) ** (1 / 3))) ** x) - 1

print("\nTesting `find_root_bisection` ...")
x= find_root_bisection(f, 0.1, 100.0, xtol=1e-6, ftol=1e-6, verbose=True)
print('%3f' % x)
#
# print("\nTesting `find_root_secant` ...")
# x= find_root_secant(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
# print('%6f' % x)
#
# print("\nTesting `find_root` with default value of `contraction_factor` ...")
# x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
# print('%6f' % x)

# print("\nTesting `find_root` with `contraction_factor= 1.0` ...")
# x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, contraction_factor=1.0, verbose=True)
# print('%6f' % x)

print(mean / math.gamma(1 + 1 / x))