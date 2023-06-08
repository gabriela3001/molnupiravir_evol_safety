import numpy as np


# BASIC MODEL
### Euler equation update Eq. (1,2) of Methods

def merck_euler(population,dt,b,q,m,n,a):
    mutation = 1
    x,y = population
    dy      = np.empty(2)
    dy[0] = x + dt*(b*x*(q**(m)*(q**n))-a*x)
    dy[1] = y + dt*(mutation*(b*x*(q**m)*(1-q**n)+b*y*(q**m)-a*y))
    return(dy)
    
def simulation(init,b,a0,a1,q0,q1,m,n,tim,ttr):
    
    dt = 1e-4
    initial = [init[0],init[1]]
    population = initial.copy()
    pop_evol = []
    
    current_q = q0
    current_a = a0

    for ngen in np.arange(30/dt):
        if ngen % 10 == 0:
            pop_evol.append(population)
        
        population = merck_euler(population,dt,b,current_q,m,n,current_a)
        
        # immune system kicks in
        if ngen > tim/dt:
            current_a = a1
            
        # Merck treatment starts
        if ngen > ttr/dt:
            current_q = q1
            
    pop_evol = np.array(pop_evol)
    return(pop_evol)


### Functions giving analytical solution along time for wild type x(t), mutant y(t) and total viral load v(t)
##### Eq. (3) (5) (6) of Methods


def func_x(t, b, a, q, m, n, x0, y0):
    return(x0*np.exp((b*q**(m+n) - a)*t))

def func_v(t, b, a, q, m, n, x0, y0):
    return(np.exp((b*(q**m)-a)*t)*(x0+y0))

def func_y(t, b, a, q, m, n, x0, y0):
    v_ =  func_v(t, b, a, q, m, n, x0, y0)
    x_ =  func_x(t, b, a, q, m, n, x0, y0)
    return(v_-x_)


### Total viral load in the growth phase
##### Eq. (4)

def growth_total_X(T, b, a, q, m, n):
    return((1/(-a+b*(q**(m+n))))*np.exp((-a+b*(q**(m+n)))*T))

def growth_total_Y(T, b, a, q, m, n):
    X_ = growth_total_X(T, b, a, q, m, n)
    V_ = growth_total_V(T, b, a, q, m, n)
    return(V_-X_)
    
def growth_total_V(T, b, a, q, m, n):
    return((1/(-a+b*(q**m)))*np.exp((-a+b*(q**m))*T))


### Total viral load in the clearance phase
##### Eq. (7) (12) (13) (14)

def clearance_total_V(T, b, a0, a1, q, m, n):
    return((1/(a1-b*(q**(m))))*np.exp((-a0+b*(q**(m)))*T))

def clearance_total_X(T, b, a0, a1, q, m, n):
    return((1/(a1-b*(q**(m+n))))*np.exp((-a0+b*(q**(m+n)))*T))

def clearance_total_Y(T, b, a0, a1, q, m, n):
    X_ = clearance_total_X(T, b, a0, a1, q, m, n)
    V_ = clearance_total_V(T, b, a0, a1, q, m, n)
    return(V_-X_)


### Total viral load produced over the course of an infection
##### Eq. (8) (9) (10) (11)

def totalV_wholeinfection(T, b, a0, a1, q0, q1, m, n):
    V = ((1/(-a0+b*(q0**(m)))) + (1/(a1-b*(q1**(m)))))*np.exp((-a0+b*(q0**(m)))*T)
    return(V)

def totalX_wholeinfection(T, b, a0, a1, q0, q1, m, n):
    X = ((1/(-a0+b*(q0**(m+n)))) + (1/(a1-b*(q1**(m+n)))))*np.exp((-a0+b*(q0**(m+n)))*T)
    return(X)

def totalY_wholeinfection(T, b, a0, a1, q0, q1, m, n):
    X = totalX_wholeinfection(T, b, a0, a1, q0, q1, m, n)
    V = totalV_wholeinfection(T, b, a0, a1, q0, q1, m, n)
    Y = V-X
    return(Y)


### Peak mutation rate u* for treatment at t = T
def calc_ustart_peaktreatment(T, b, a0, a1, q0, m, n):
    vT = np.exp(b*q0**m - a0)
    xT = np.exp(b*q0**(m+n) - a0)
    yT = vT-xT
    eta = yT/xT
    u_star = ((a1-b)/(m*b))*((n-eta*m)/(n+eta*m))
    return(u_star)

### Peak mutation rate u* for treatment at t = 0
def calc_F(T, b, a0, a1, m, u):
    h = b*T
    k = (b*(2*b-a0-a1))/((b-a0)*(a1-b))
    mu = m*u
    F = h+k-mu*((h**2)+(k**2))-((mu**2)*h*k*(2*h+k))-(mu**3)*(h**2)*(k**2)
    return(F)

def calc_k(b,a0,a1):
    k = (b*(2*b-a0-a1))/((b-a0)*(a1-b))
    return(k)

def calc_h(b,T):
    h = b*T
    return(h)

def solve_numerically_F(T, b, a0, a1, m, u_vals):
    return(u_vals[np.searchsorted(-calc_F(Tval, b, a0, a1val, mval, u_vals), 0)])