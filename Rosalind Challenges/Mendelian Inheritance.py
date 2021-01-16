""" k, m, n
k -> k/(K+n+m) * ((k-1)+(m)+(n)/(k-1+n+m)
m -> m/(k+m+n) *

27 23 20
--only homozygote recessive (n)--> n/(k+m+n) * ((n-1)/(k+m-1+n) + 1/2 * m/
"""

k = 27
m = 23
n = 20

rv = 1- (n*(n-1) + (1/4)*m*(m-1) + 2 * (1/2)*m*n)/((k+m+n)*(k+m+n-1))

print(rv)