import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy.linalg import eigvals
#from tqdm.auto import tqdm
import numpy.linalg as nlg

############################################################################
def routine_diag(H):
    eigv = nlg.eigvals(H)
    return eigv

############################################################################

def GUE_matrices(n,ens):
    eigenvalues_ensemble = []
    spacing_ratios = []
    for i in range(ens):
    # Generate a random Hermitian matrix
        H0 = np.random.normal(0, 1/np.sqrt(2), (n,n)) + 1j*np.random.normal(0, 1/np.sqrt(2), (n,n))
        G = (H0 + H0.conj().T) / 2  # Make it Hermitian
        # Calculate the eigenvalues
        eigvs = routine_diag(G)
        eigso = np.sort(eigvs) #sort the eigenvalues
        # Append the eigenvalues to the ensemble
        eigenvalues_ensemble.append(eigvs) 
        for j in range(len(eigvs)-2):
        #spacings = np.diff(eigs)
            sratio = (eigso[j+2]-eigso[j+1])/(eigso[j+1]-eigso[j])
            spacing_ratios.append(sratio)
    return eigenvalues_ensemble, spacing_ratios

n = 1000  # Matrix size
ens = 50  # Number of random matrices

eigv,SR = GUE_matrices(n,ens)

###########################################################################

e1 = np.array(eigv).flatten()
e2 = np.array(SR).flatten()

##############Formulas of spacing ratios###################################
x = np.linspace(0, 10, 1000)
prgoe = 27*x*(x+1)/(8*(x** 2 + x + 1)**(5/2)) 
prgue = 81*np.sqrt(3)*x**2*(x+1)**2/(4*np.pi*(x ** 2 + x + 1)**(4)) 
prpoi = 1/(x+1)**2

###################Plotting of the results#####################################
plt.hist(e1,range = [-,80],  bins=50, density=True, alpha=0.5, color='blue')
#plt.savefig('GOE-ev-clustertrailN1000.pdf')
################################################################################

plt.hist(e2,range = [0,8],  bins=60,density=True, alpha=0.3, color='blue')
#plt.hist(es,range = [0,8],  bins=40,density=True, alpha=0.3, color='orange')
plt.plot(x,prgoe,alpha=0.8, color='black')
plt.plot(x,prgue,alpha=0.8, color='red')
plt.plot(x,prpoi,alpha=0.8, color='blue')
plt.xlabel('r', fontsize = 18)
plt.ylabel('P(r)', fontsize = 16)
plt.title('Histogram of Spacing Ratios of Eigenvalues')
# Function add a legend  
plt.legend(["GOE", "GUE", "Poisson"], loc ="upper right", fontsize = 12)
plt.grid(False)
#plt.savefig('GOE-sprat-clustertrailN1000.pdf')

########################################################################