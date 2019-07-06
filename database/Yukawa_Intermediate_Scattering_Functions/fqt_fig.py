''' plot F(q, t) ;
input: kappa Gamma nq
to run type "python fqt_fig.py kappa Gamma nq"
kappa, Gamma: plasma parameters
nq: wave number index, should be an integer
    0<= nq <= 165
    q = (nq+1)*dq, dq = 0.181

 Yongjun Choi
'''

import matplotlib.pyplot as plt
import numpy as np
import sys

# input
kappa = int(sys.argv[1])
gamma = int(sys.argv[2])
nq = int(sys.argv[3])

# Load F(q, t) data;
fqt = np.load("fqt_k"+str(kappa)+"G"+str(gamma)+".npy")

t = fqt[0, 1:]  # timestamps
q = fqt[1:, 0]  # wavenumbers

plt.figure()
plt.plot(t, fqt[nq+1, 1:], 'b')
plt.xlim(0, 100)
plt.xlabel("$t$")
plt.ylabel("$F(q, t)$")
plt.title("$\kappa =$"+str(kappa)+", $\Gamma =$"+str(gamma)+", $q = $"+"${0:.2f}$".format(q[nq]))
plt.show()

