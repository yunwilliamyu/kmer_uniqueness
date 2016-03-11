#!/usr/bin/env python
# kmer_analysis.py
# Generates figure from kmer_analysis.log type files

import fileinput
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ksize = []
totalloci = []
uniquekmers = []
nonuniqueloci = []
uniqueloci = []
hammingloci = []
unambigloci = []
for line in fileinput.input():
    if (line == "========\n"):
        pass
    elif (line.startswith("k = ")):
        ksize.append(int(line.replace("k = ","")))
    elif (line.strip().startswith("Total loci:")):
        totalloci.append(int(line.replace("Total loci:","")))
    elif (line.strip().startswith("Unique k-mers:")):
        uniquekmers.append(int(line.replace("Unique k-mers:","")))
    elif (line.strip().startswith("Non-unique loci:")):
        nonuniqueloci.append(int(line.replace("Non-unique loci:","")))
    elif (line.strip().startswith("Unique loci:")):
        uniqueloci.append(int(line.replace("Unique loci:","")))
    elif (line.strip().startswith("Hamming-ed locations:")):
        pass
        #hammingloci.append(int(line.replace("Hamming-ed locations:","")))
    elif (line.strip().startswith("Unambiguous locations:")):
        unambigloci.append(int(line.replace("Unambiguous locations:","")))

        
print ksize
print totalloci
print uniquekmers
print nonuniqueloci
print uniqueloci
#print hammingloci
print unambigloci

plt.ioff()
#plt.plot([1.6,2.7])
plt.plot(ksize, totalloci, 'ro-', label="Total loci")
plt.plot(ksize, uniqueloci, 'bs-', label="Unambig loci (exact match)")
plt.plot(ksize, unambigloci, 'g*-', label="Unambig loci (Hamming dist 1)")
plt.grid(True)
plt.xlabel('k-mer size')
plt.ylabel('Number of loci')
plt.title('Identifying loci from k-mers')
plt.legend(loc=4)


plt.savefig('foo.png')

