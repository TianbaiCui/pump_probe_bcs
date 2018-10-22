import sys, os, shutil
import time, math
import multiprocessing

from numpy import *


L_kxy = 400
L_kx = 400
mu = - 1.18992
#A0 = linspace(0.1, 1.6, 16)# sqrt(number * Delta0)
A0 = 1.12 #0.4 / 3.07652#1.05/1.3 #0.2294 #0.6*0.6

tau = 10#0.82 #3.076#20.5101 #10.4616 #5 * 2 * pi #2 * pi # number divided by Delta0
sigma = tau/5.0 #/ 2 # number divided by Delta0
omega_p = 1.3*0.0843672 #0.120865

T1 = 5  #1.0 # number * tau
T2 = 0.5 #0.35 # number *tau

x_max = 4 #number *tau

os.system("./main %f %f %f %f %f %f %f %f %f %f" %(L_kxy, L_kx, mu, A0, tau, sigma, omega_p, T1, T2, x_max))

#def run(i):
#    return os.system("./main %f %f %f %f %f %f %f %f %f %f" %(L_kxy, L_kx, mu, A0[i], tau, sigma, omega_p, T1, T2, x_max))

#if __name__ == '__main__':
#    pool = multiprocessing.Pool()
#    pool.map(run, range(len(A0)))

#for i in range(len(A0)):
#    os.system("./main %f %f %f %f %f %f %f %f %f %f" %(L_kxy, L_kx, mu, A0[i], tau, sigma, omega_p, T1, T2, x_max))

#os.system("./main %f %f %f %f %f %f %f %f %f %f" %(L_kxy, L_kx, mu, A0, tau, sigma, omega_p, T1, T2, x_max))


#T1_table = [0., tau, 2.5*tau, 10*tau]

#print "Hello! I'm rank %d from %d running in total..." %(comm.rank, comm.size)
#print T1_table

#for idx in range(0, 4):
        #os.system("./bcsdyn %f %f %f" %(A0, T1_table[idx], T1_table[idx]))
