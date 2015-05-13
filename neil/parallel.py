#!/usr/local/bin/python
import numpy
import time
from multiprocessing import Pool
from functools import partial

def test_prime(n):
    prime = True
    for i in range(2,n):
        if n % i == 0:
            prime = False
    numpy.sqrt(2)
    return prime

def parallel_function(f):
    def easy_parallize(f, sequence,pool_size=8):
        """ assumes f takes sequence as input, easy w/ Python's scope """
        from multiprocessing import Pool
        pool = Pool(processes=pool_size) # depends on available cores
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = numpy.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallize, f)

N = 10000

"""
start = time.time()
serial = []
for i in range(N):
    serial.append(test_prime(i))
serial = numpy.asarray(serial)
end = time.time()
print end-start
"""

for i in [16]:
    start = time.time()
    parallel = parallel_function(test_prime)
    parallel_result = parallel(range(N),pool_size=i)
    end = time.time()
    print i,end-start

def a(a,b=True,c=2):
    print a,b,c
    return a

def a_helper(args):
    return partial(a,**args)

x = [1,2,3]
arg_list = {'b':False,'c':5}
print map(a,x)
print map(a_helper(arg_list),x)
