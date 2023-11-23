import sys
import matplotlib.pyplot as plt
import numpy as np
import csv

Ndata = []
Nndata = []

for arg in range(len(sys.argv)-1):

    filename = str(sys.argv[arg+1])

    E = []
    N = []
    Nnorm = []

    print('Opening.... ' + filename)

    with open(filename) as File:
        mat = next(File)
        size = next(File)
        nevents = next(File)
        next(File)
        plots = csv.reader(File, delimiter=',')

        for row in plots:
            E.append(float(row[1]))
            N.append(int(row[2]))
            Nnorm.append(int(row[2])/int(nevents))

    total = sum(N)
    
    n_arr = np.array(N)
    nnorm_arr = np.array(Nnorm)

    Edata = np.transpose(np.array(E))
    #Ndata = np.hstack((Ndata,np.transpose(n_arr)))
    #Nndata = np.hstack((Nndata,np.transpose(nnorm_arr)))

    #plt.semilogx(E, nnorm_arr, label=size.strip()+'cm; '+str(total)+'/'+nevents.strip()+' events')
    plt.semilogx(E, nnorm_arr, label=mat.strip()+'; '+str(total)+'/'+nevents.strip()+' events')


#plt.semilogx(E, n_arr, color='green')
#plt.xscale('log')
plt.xlabel('Energy (MeV)')
plt.xlim([1e-8,10])
plt.ylabel('Counts per primary event')
#plt.title('2 MeV neutron beam onto ' + mat.strip())
plt.title('2 MeV neutron beam onto 30 cm')
plt.legend()
plt.show()
