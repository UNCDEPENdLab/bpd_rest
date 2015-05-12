#!/usr/local/bin/python
import numpy
import time

def test_prime(n):
    prime = True
    for i in range(2,n):
        if n % i == 0:
            prime = False
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

N = 80000

"""
start = time.time()
serial = []
for i in range(N):
    serial.append(test_prime(i))
serial = numpy.asarray(serial)
end = time.time()
print end-start
"""

for i in [16,24,32]:
    start = time.time()
    parallel = parallel_function(test_prime)
    parallel_result = parallel(range(N),pool_size=i)
    end = time.time()
    print i,end-start
