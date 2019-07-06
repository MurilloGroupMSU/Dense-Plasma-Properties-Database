''' plot S(q, w) ;
input: kappa Gamma nq
to run type "python sqw_fig.py kappa Gamma nq"
kappa, Gamma: plasma parameters
nq: wave number index, should be an integer
    0<= nq <= 165
    q = (nq+1)*dq, dq = 0.181

 Yongjun Choi
'''

import matplotlib.pyplot as plt
import sys
import numpy as np

# input
kappa = int(sys.argv[1])
gamma = int(sys.argv[2])
nq = int(sys.argv[3])

# Load sqw data;
sqw = np.load("sqw_k"+str(kappa)+"g"+str(gamma)+".npy")

w = sqw[0, 1:]  # frequencies
q = sqw[1:, 0]  # wavenumbers

dq = q[1] - q[0]
wavenumber = dq*(nq+1)

plt.figure()
plt.plot(w, sqw[nq+1, 1:])
plt.xlabel("$\omega$")
plt.ylabel("$S(q, \omega)$")
plt.title("$\kappa =$"+str(kappa)+", $\Gamma =$"+str(gamma)+", $q = $"+"${0:.2f}$".format(wavenumber))
plt.show()
