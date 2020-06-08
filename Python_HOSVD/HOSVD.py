from sktensor.tucker import hooi, hosvd
from sktensor.core import ttm
import numpy as np
from numpy import linalg as la
from sktensor import dtensor
from scipy.io.matlab import loadmat
import matplotlib.pyplot as plt
import time

start_time = time.time()
#Load MAT file of microseismic data
tensor_xyz = loadmat('xyztensor')
tensor = loadmat('tensor_test')

#Convert to tensor
T = dtensor(tensor["tensor"])
#T = dtensor(tensor_xyz["tensor_xyz"])
T_arr=np.array(T)

##Tolerance
#eps_0=[9, 2, 1] #0.2
#eps_1=[32, 8, 5] #0.02
#eps_2=[84, 11, 12] #0.002
#eps_3=[181, 12, 31] #0.0002
eps_3=[181, 34, 11] #0.0002

## HOOI
core,U= hooi(T, eps_3, init='nvecs')
#Factor matrices and core tensor
core_S = np.array(core)
print(core_S.shape)

U1 = U[0]
U2 = U[1]
U3 = U[2]

Trec = ttm(core, U)
Trec_S = np.array(Trec)

rel_er=la.norm(T-Trec)/la.norm(T)  #Relative Error
print(rel_er)

#Compression rate
CR=(T.shape[0]*T.shape[1]*T.shape[2])/(core_S.shape[0]*core_S.shape[1]*core_S.shape[2])
print(100-CR)

#Frobenius norm
sv1=np.zeros(core_S.shape[0])
sv2=np.zeros(core_S.shape[1])
sv3=np.zeros(core_S.shape[2])
for i in range(core_S.shape[0]):
    sv1[i]=la.norm(core_S[i,:,:], 'fro')
for i in range(core_S.shape[1]):
    sv2[i]=la.norm(core_S[:,i,:], 'fro')
for i in range(core_S.shape[2]):
    sv3[i]=la.norm(core_S[:,:,i], 'fro')

#Normalized singular values
sv1_n=sv1/np.max(sv1)
sv2_n=sv2/np.max(sv2)
sv3_n=sv3/np.max(sv3)
sv=[sv1_n,sv2_n,sv3_n]

#Elapsed time
elapsed_time = time.time() - start_time
print(elapsed_time)

##Plotting
#Normalized singular values
fig1, (ax0, ax1, ax2) = plt.subplots(1, 3)
ax0.plot(sv1_n)
ax0.set_title('Mode 1',fontsize=10)

ax1.plot(sv2_n)
ax1.set_title('Mode 2',fontsize=10)

ax2.plot(sv3_n)
ax2.set_title('Mode 3',fontsize=10)

fig1.suptitle('Normalized singular values', fontsize=16)

#Comparison between Recovered and Original
k=19
a=Trec[:,:,k]
b=T[:,:,k]
print(type(a))
print(b.shape)

fig2, (ax0, ax1) = plt.subplots(1, 2)

c = ax0.pcolor(a)
ax0.set_title('Recovered')
ax0.invert_yaxis()

c = ax1.pcolor(b)
ax1.set_title('Original')
ax1.invert_yaxis()

fig2.suptitle('Comparison between Recovered and Original', fontsize=16)
fig2.set_size_inches(7,7)
plt.show()

#General Structure: Microseismic traces
fig3, axs = plt.subplots(3, 2)

ax0 = axs[0, 0]
ax0.plot(Trec[:,3,6:8])
ax0.set_title('Reduced structure of microseismic data',fontsize=10)
ax0.set_xticks([])
ax0.set_yticks([])
ax1 = axs[1, 0]
ax1.plot(Trec[:,3,6:8])
ax1.set_xticks([])
ax1.set_yticks([])
ax2 = axs[2, 0]
ax2.plot(Trec[:,3,6:8])
ax2.set_xticks([])
ax2.set_yticks([])

ax3 = axs[0, 1]
ax3.set_title('Structure of microseismic data',fontsize=10)
ax3.plot(T[:,3,6:8])
ax3.set_xticks([])
ax3.set_yticks([])
ax4 = axs[1, 1]
ax4.plot(T[:,3,6:8])
ax4.set_xticks([])
ax4.set_yticks([])
ax5 = axs[2, 1]
ax5.plot(T[:,3,6:8])
ax5.set_xticks([])
ax5.set_yticks([])

fig3.suptitle('General Structure: Microseismic traces', fontsize=16)
plt.show()

#Number of variables
fig4, (ax0, ax1) = plt.subplots(1,2)
p=[CR, 100-CR]
x_b = np.arange(2)
ba=[T.size, T.size*(CR/100)]

ax0.pie(p, autopct='%1.1f%%')
ax0.set_title('Variables',fontsize=10)

ax1.bar(x_b[0], T.size, label='Tensor')
ax1.bar(x_b[1], T.size*(CR/100), label='HOSVD')
ax1.set_title('Number of variables',fontsize=10)
ax1.legend()

fig4.suptitle('Number of variables', fontsize=16)
plt.show()
