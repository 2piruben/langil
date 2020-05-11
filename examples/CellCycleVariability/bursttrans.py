import numpy as np
from math import factorial

def f(p,k,d):
        return k[:-1]/(k[1:]+p*d[1:])   
    # f has one elements rest than the rest

def sumfact(p,beta,r,k,d):

    s = np.zeros(len(k)-1)
    for j in range(p):
        s = s + factorial(p)/factorial(j)*beta[1:]**(p-j-1)*factorial_moment(j,beta,r,k,d)[1:]
    return s

def gtilde(p,beta,r,k,d):

    # print('gtilde: ',sumfact(p,beta,r,k,d),'/',(k[1:]+p*d[1:]))
    return r[1:]*sumfact(p,beta,r,k,d)/(k[1:]+p*d[1:])

def g(p,beta,r,k,d):

    # print('gtilde: ',sumfact(p,beta,r,k,d),'/',(k[1:]+p*d[1:]))
    return r[1:]*p*factorial_moment(p-1,beta,r,k,d)[1:]/(k[1:]+p*d[1:])

def delta(p,k,d):
# delta has N-1 elements, same than k
    delta_array = f(p,k,d)
    for i in range(1,len(delta_array)):
        delta_array[i] = delta_array[i]*delta_array[i-1]
    return delta_array

def thetatilde(p,beta,r,k,d):

    frac = gtilde(p,beta,r,k,d)/delta(p,k,d)
    return delta(p,k,d)*np.cumsum(frac)

def theta(p,beta,r,k,d):

    frac = g(p,beta,r,k,d)/delta(p,k,d)
    return delta(p,k,d)*np.cumsum(frac)

def pi_pop_Erlang(i,N): # probability of finding the cell in state i if the states are exp distributed with same rate
	return (2**(1/N)-1)/(2**((i/N)-1))

def factorial_moment(p,beta,r,k,d,mode = 'constitutive'):
    # mode can be constitutive or bursty

    if (p==0):

        nf0 = np.zeros(len(k))
        nf10 = 0
        frac = k[1:]/k[:-1]
        for i in range(1,len(k)):
            p = 1
            for j in range(1,i):
                p = p*frac[j]
            nf10 += p

        nf0[0] = 1.0/(1 + nf10)
        
        for i in range(1,len(k)):
            # print('i',i)
            # print('nf0',nf0)
            nf0[i] = nf0[i-1]*frac[i-1]
        return nf0

    else:

        nfp = np.zeros(len(k))
        s = 0
        delt = delta(p,k,d)
        if mode == 'constitutive':
            theta_mode = theta(p,beta,r,k,d)
            s = p*factorial_moment(p-1,beta,r,k,d)[0]
        elif mode == 'bursty':
            theta_mode = thetatilde(p,beta,r,k,d)
            for j in range(p):
                s = s + factorial(p)/factorial(j)*beta[0]**(p-j-1)*factorial_moment(j,beta,r,k,d)[0]

        # print("sum for n1 {}\n".format(s))
        # print("first sum {}\n".format(2*r[0]*s))
        # print("second sum {}\n".format(k[-1]*0.5**(p-1)*thetatild[-1]))
        nf1p = 2*r[0]*s + k[-1]*0.5**(p-1)*theta_mode[-1]
        nf1p = nf1p/( 2*(d[0]*p + k[0]) - k[-1]*0.5**(p-1)*delt[-1])
        nfp[0] = nf1p
        nfp[1:] = delt*nf1p + theta_mode

        return nfp

def meantrajErlang(beta,r,k,d,W,N):
# mean for the Erlang distribution with parameters beta,r,k,d, with N phases, W of them premitotic
    r_vec = r*np.ones(N)
    k_vec = k*np.ones(N)
    r_vec[W:] = r_vec[W:]*2 
    d_vec = d*np.ones(N)
    beta_vec = beta*np.ones(N)
    return sum(factorial_moment(1,beta_vec,r_vec,k_vec,d_vec))


def vartrajErlang(beta,r,k,d,W,N, mode = 'bursty'):
# mean for the Erlang distribution with parameters beta,r,k,d, with N phases, W of them premitotic
    r_vec = r*np.ones(N)
    k_vec = k*np.ones(N)
    r_vec[W:] = r_vec[W:]*2 
    d_vec = d*np.ones(N)
    beta_vec = beta*np.ones(N)

    n1 = factorial_moment(1,beta_vec,r_vec,k_vec,d_vec, mode = mode)
    # print('n1',n1)
    n2 = factorial_moment(2,beta_vec,r_vec,k_vec,d_vec, mode = mode)
    # print('n2',n2)
    var = sum(n2) + sum(n1) - sum(n1)*sum(n1)
    return var


def meanpopErlang(beta,r,k,d,W,N, mode = 'bursty'):
# mean for the Erlang distribution with parameters beta,r,k,d, with N phases, W of them premitotic
    r_vec = r*np.ones(N)
    k_vec = k*np.ones(N)*2**(1/N) # this is changed to get pop measurements
    r_vec[W:] = r_vec[W:]*2 
    d_vec = d*np.ones(N)
    beta_vec = beta*np.ones(N)

    n1 = factorial_moment(1,beta_vec,r_vec,k_vec,d_vec,mode = mode)/factorial_moment(0,beta_vec,r_vec,k_vec,d_vec,mode = mode)
    return sum(n1*pi_pop_Erlang(np.arange(N)+1,N))


def varpopErlang(beta,r,k,d,W,N,mode = 'bursty'):
# mean for the Erlang distribution with parameters beta,r,k,d, with N phases, W of them premitotic
    r_vec = r*np.ones(N)
    k_vec = k*np.ones(N)*2**(1/N)
    r_vec[W:] = r_vec[W:]*2 
    d_vec = d*np.ones(N)
    beta_vec = beta*np.ones(N)

    n1 = factorial_moment(1,beta_vec,r_vec,k_vec,d_vec,mode = mode)/factorial_moment(0,beta_vec,r_vec,k_vec,d_vec,mode = mode)
    sum1 = np.sum(n1*pi_pop_Erlang(np.arange(N)+1,N))
    # print('n1',n1)
    n2 = factorial_moment(2,beta_vec,r_vec,k_vec,d_vec,mode = mode)/factorial_moment(0,beta_vec,r_vec,k_vec,d_vec,mode = mode)
    sum2 = np.sum(n2*pi_pop_Erlang(np.arange(N)+1,N))
    # print('n2',n2)
    var = sum2 + sum1 - sum1*sum1
    return var







