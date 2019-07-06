''' plot susceptibility chi(q, w) ;
input: kappa Gamma nq
to run type "python chi_qw_fig.py kappa Gamma nq"
kappa, Gamma: plasma parameters
nq: wave number index, should be an integer
    0<= nq <= 165
    q = (nq+1)*dq, dq = 0.1809

 Yongjun Choi
'''

import matplotlib.pyplot as plt
import sys
import numpy as np
import sub

# input

path0 ="/Users/yongjunchoi/Research/Yukawa/data/Sets/T800dt001/"
r_path0=path0+"/kappa_0/gamma_10/"
kappa = int(sys.argv[1])
gamma = int(sys.argv[2])
nq = int(sys.argv[3])
#set working path
(N, Neq, Npd, q_max, T, L, Nq, w1, t1, dq, dw, dt) = \
    sub.checkSanity(r_path0, kappa, gamma)


# Load chi_qw data;
chi_qw = np.load("chi_qw_k"+str(kappa)+"g"+str(gamma)+".npy")

w = np.real(chi_qw[0, 1:])  # frequencies
q = np.real(chi_qw[1:, 0])  # wavenumbers
plt.figure()
plt.plot(w, np.real(chi_qw[nq+1, 1:]), 'b', label="$\chi'(q, \omega)$")
plt.plot(w, np.imag(chi_qw[nq+1, 1:]), 'r', label="$\chi''(q, \omega)$")

plt.legend(loc = "best")
plt.xlabel("$\omega$")
plt.ylabel("$\chi(q, \omega)$")
plt.title("$\kappa =$"+str(kappa)+", $\Gamma =$"+str(gamma)+", $q = $"+"${0:.2f}$".format(q[nq]))
plt.show()

