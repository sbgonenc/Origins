"1 1 1+3 1+3+3 1+3+3+3+9"

#Using a decorator
from functools import lru_cache
"""
@lru_cache(maxsize= 1000)
def fibonacci2(n):
	if n == 1:
		rv = 1
	elif n == 2:
		rv = 1
	elif n> 2:
		rv = fibonacci2(n-2) + fibonacci2(n-1)
	return rv

for i in range(1, 1001):
	#print(i, ":", fibonacci2(i))


#Memoization_ Using dictionary
fibonacci_cache = {}
def fibonacci(n):
	if n in fibonacci_cache:
		return fibonacci_cache[n]
	if n == 1:
		rv = 1
	elif n == 2:
		rv = 1
	elif n> 2:
		rv = fibonacci(n-2) + fibonacci(n-1)
		fibonacci_cache[n] = rv
	return rv

for i in range(1, 10001):
	print(i, ":", fibonacci(i))
"""

@lru_cache(maxsize = 1000)
def Wascally(months, pairs):
	if months == 1:
		rv = 1
	elif months == 2:
		rv = 1
	elif months > 2:
		rv = Wascally(months-2, pairs)*pairs + Wascally(months-1,pairs)
	return rv

for i in range(1,101):
	print(i, ":", Wascally(i,3))