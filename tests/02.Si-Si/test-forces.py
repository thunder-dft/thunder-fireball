#!/usr/bin/env python
import re
import numpy as np

# reading files and storing energies and forces
etot = []
f_tot = []
file_etot = open('etot.dat', 'w')
file_ftot = open('f_total.dat', 'w')

ebs = []
f_ebs = []
f_pulay =[]
file_ebs = open('ebs.dat', 'w')
file_febs = open('f_ebs.dat', 'w')
file_fpulay = open('f_pulay.dat', 'w')

usr =[]
f_usr = []
file_usr = open('usr.dat', 'w')
file_fusr = open('f_usr.dat', 'w')

uxcdcc =[]
file_uxcdcc = open('uxcdcc.dat', 'w')

components = []

kinetic = []
f_kinetic = []
file_KE = open('kinetic.dat', 'w')
file_fKE = open('f_kinetic.dat', 'w')

vna = []
f_vna = []
file_VNA = open('vna.dat', 'w')
file_fVNA = open('f_vna.dat', 'w')

file_VXC = open('vxc.dat', 'w')
file_VNL = open('vnl.dat', 'w')
file_EWDSR = open('ewaldsr.dat', 'w')
file_EWDLR = open('ewaldlr.dat', 'w')

for n in range(1, 200):
        # write out energies to .txt files from *.log files 
        if n < 10:
            source = '00%d.log' % n
        if n >= 10 and n < 100:
            source = '0%d.log' % n
        if n >= 100:
            source = '%d.log' % n

        file_source = open(source)
        for line in file_source:
            line = line.rstrip()
            
            if re.search('ETOT =', line):
                data = line.split()
                file_etot.write('%d %f\n' % (n, float(data[2])))
                etot.append(data[2])
            if re.search('ebs =', line):
                data = line.split()
                file_ebs.write('%d %f\n' % (n, float(data[2])))
                ebs.append(data[2])
            if re.search('uii - uee =', line):
                data = line.split()
                file_usr.write('%d %f\n' % (n, float(data[4])))
                usr.append(data[4])
            if re.search('uxcdcc =', line):
                data = line.split()
                file_uxcdcc.write('%d %f\n' % (n, float(data[2])))
                uxcdcc.append(data[2])
                
            # write out force components to .txt files
            if re.search('f_ebs    1', line):
                data = line.split()
                file_febs.write('%d %f\n' % (n, float(data[5])))
                f_ebs.append(data[5])     
            if re.search('f_usr    1', line):
                data = line.split()
                file_fusr.write('%d %f\n' % (n, float(data[5])))
                f_usr.append(data[5])
            if re.search('f_pulay    1', line):
                data = line.split()
                file_fpulay.write('%d %f\n' % (n, float(data[5])))
                f_pulay.append(data[5])
            if re.search('f_kinetic    1', line):
                data = line.split()
                file_fKE.write('%d %f\n' % (n, float(data[5])))
                f_kinetic.append(data[5])
            if re.search('f_vna_2c    1', line):
               data = line.split()
               file_fVNA.write('%d %f\n' % (n, float(data[5])))
               f_vna.append(data[5])            
            if re.search('f_total =     1', line):
                data = line.split()
                file_ftot.write('%d %f\n' % (n, float(data[5])))
                f_tot.append(data[5])

        if n < 10:
            source = '00%d.ENERGIES' % n
        if n >= 10 and n < 100:
            source = '0%d.ENERGIES' % n
        if n >= 100:
            source = '%d.ENERGIES' % n

        file_source = open(source)
        components = np.loadtxt(file_source,skiprows=1)
        
        # write out the energy components to the different files and to their
        # respective arrays
        kinetic.append(components[1])        
        file_KE.write('%d %f\n' % (n, float(components[1]))) 
        
        vna.append(components[2])
        file_VNA.write('%d %f\n' % (n, float(components[2])))
        
        file_VXC.write('%d %f\n' % (n, float(components[3])))
        file_VNL.write('%d %f\n' % (n, float(components[4])))
        file_EWDSR.write('%d %f\n' % (n, float(components[5])))
        file_EWDLR.write('%d %f\n' % (n, float(components[6])))
# end reading files        
        
# Now calculate the deltaE contributions and compare these with the forces. 
#print(' ')
#print(f"deltaKE/deltaR{' '*10}forceKE")
#for i in range(1,len(kinetic)):
#        deltaKE = (float(kinetic[i]) - float(kinetic[i-1]))/0.01
#        forceKE = 0.5*(float(f_kinetic[i]) + float(f_kinetic[i-1]))
#        print("{:.4f}".format(deltaKE), ' '*15, "{:.4f}".format(forceKE))  

print(' ')
print(f"deltaVNA/deltaR{' '*10}forceVNA")
for i in range(1,len(vna)):
        deltaVNA = (float(vna[i]) - float(vna[i-1]))/0.01
        forceVNA = 0.5*(float(f_vna[i]) + float(f_vna[i-1]))
        print("{:.4f}".format(deltaVNA), ' '*15, "{:.4f}".format(forceVNA))  

print(' ')
print(f"deltaEBS/deltaR{' '*10}forceEBS")
for i in range(1,len(ebs)):
       deltaEBS = (float(ebs[i]) - float(ebs[i-1]))/0.01       
       forceEBS = 0.5*(float(f_ebs[i]) + float(f_ebs[i-1]))
       print("{:.4f}".format(deltaEBS), ' '*15, "{:.4f}".format(forceEBS))
        
print(' ')      
print(f"deltaUSR/deltaR{' '*10}forceUSR")
for i in range(1,len(usr)):
       deltaUSR = (float(usr[i]) + float(uxcdcc[i]) \
                   - float(usr[i-1]) - float(uxcdcc[i-1]))/0.01
       forceUSR = 0.5*(float(f_usr[i]) + float(f_usr[i-1]))
       print("{:.4f}".format(deltaUSR), ' '*15, "{:.4f}".format(forceUSR))
              
print(' ')
print(f"deltaETOT/deltaR{' '*10}forceTOT")
for i in range(1,len(etot)):
        deltaETOT = (float(etot[i]) - float(etot[i-1]))/0.01
        forceTOT = 0.5*(float(f_tot[i]) + float(f_tot[i-1]))
        print("{:.4f}".format(deltaETOT), ' '*15, "{:.4f}".format(forceTOT))        