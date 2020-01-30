import numpy as np
x_new = []
y_new = []

x_new_four = []
y_new_four = []

x_expected = []
y_expected = []

psi_new_1 = []
psi_new_2 = []
psi_new_3 = []
psi_new_4 = []


def eqn(x,y):
    u1 = -0.005*x*x + 0.01*y*y
    u2 = -0.01*x*x+0.005*y*y
    return u1,u2

def strain(x,y):
    psi_1 = -0.01*x
    psi_2 =  0.02*y
    psi_3 = -0.02*x
    psi_4 =  0.01*y
    return psi_1,psi_2,psi_3,psi_4

x1 = np.array([0,1,1,0,0.5,1,0.5,0,0.5])
y1 = np.array([0,0,1,1,0,0.5,1,0.5,0.5])

x_new,y_new = eqn(x1,y1)
print(x_new,y_new)

print()
x2 = np.array([0,0.25,0.5,0.75,1,1,1,1,1,0.75,0.5,0.25,0,0,0,0])
y2 = np.array([0,0,0,0,0,0.25,0.5,0.75,1,1,1,1,1,0.75,0.5,0.25])

x_new_four,y_new_four = eqn(x2,y2)
print(x_new_four,y_new_four)

x_middle = np.array([0.25,0.5,0.75,0.25,0.5,0.75,0.25,0.5,0.75])
y_middle = np.array([0.25,0.25,0.25,0.5,0.5,0.5,0.75,0.75,0.75])

print()
x_expected,y_expected = eqn(x_middle,y_middle)
print(x_expected,y_expected)

x_psi = np.array([0,0.5,1,1,1,0.5,0,0])
y_psi = np.array([0,0,0,0.5,1,1,1,0.5])

print()
psi_new_1,psi_new_2,psi_new_3,psi_new_4 = strain(x_psi,y_psi)
print(psi_new_1,psi_new_2,psi_new_3,psi_new_4)
print()

psi_expe_1,psi_expe_2,psi_expe_3,psi_expe_4 = strain(0.5,0.5)
print(psi_expe_1,psi_expe_2,psi_expe_3,psi_expe_4)
