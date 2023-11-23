import sys
import matplotlib.pyplot as plt
import numpy as np
import csv

E = []
N = []
Enorm = []
Nnorm = []
N2 = []

filename = str(sys.argv[1])

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

plt.bar(E, Nnorm, widths, color='green',label='Total num.: '+str(total))
plt.xscale('log')
plt.xlabel('Energy (MeV)')
plt.xlim([1e-8,10])
plt.ylabel('Counts per primary event')
plt.title('2 MeV neutron beam onto '+ size.strip() + 'cm ' + mat.strip() + '\n('+ nevents.strip() + ' events)')
plt.legend()
plt.show()

#hist = np.histogram(N,E)
#Ea = np.array(E)
#widths = (Ea[1:] - Ea[:-1])
#hist_norm = hist[0]/widths

#plt.bar(Ea[:-1],hist_norm,widths)
#plt.xscale('log')
#plt.show()