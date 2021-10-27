import numpy as np
import math
import cmath


def alpha_func(eps):
    return (np.sqrt(eps) + 1)/(np.sqrt(eps)-1)

def k_n_func(n, alpha, eps, a):
    return 1/(2*a*np.sqrt(eps))*complex((np.pi*n),math.log(alpha))

def AGF(a, q, alpha, z, z_p):
    i_q    = complex(0,q)
    exp_1  = cmath.exp(i_q*(a-z))
    exp_2  = cmath.exp(-i_q*(a-z))
    exp_3  = cmath.exp(i_q*(a+z_p))
    exp_4  = cmath.exp(-i_q*(a+z_p))
    exp_den_1  = cmath.exp(i_q*a)
    exp_den_2  = cmath.exp(-i_q*a)
    first_term  = (exp_1 + alpha*exp_2)/(exp_den_1 + alpha*exp_den_2)
    second_term = (exp_3 + alpha*exp_4)/(-exp_den_1 + alpha*exp_den_2)
    return  1/(2*i_q)*first_term*second_term
    
def eigen_funct(eps, b_n, k_n, z, n):
    delta = complex(0, np.sqrt(eps)*k_n*z)
    return(b_n*cmath.exp(delta) + (-1)**n*(cmath.exp(-delta)))

def PFGF(eigen_1, eigen_2, k, k_n, n):
    summe = 0
    for j in range(n):
        summe = summe + (eigen_1*eigen_2)/(2*k*(k-k_n))
        #print(summe)
    return(summe)

eps   = 2
alpha = alpha_func(eps)
z     = 0.2
z_p   = 0.8
n     = 10
a     = 1
k     = 2
q     = k*np.sqrt(eps)
k_n   = k_n_func(n, alpha, eps, a)
i     = complex(0,1)
b_n   = (-i)**n/(2*np.sqrt(a*eps))
eigen_1 = eigen_funct(eps, b_n, k_n, z, n)
eigen_2 = eigen_funct(eps, b_n, k_n, z_p, n)

print(AGF(a, q, alpha, z, z_p))
print(PFGF(eigen_1, eigen_2, k, k_n, n))
