#!/usr/bin/env python
import re
import numpy as np

# reading files and storing energies and forces
etot = []
fO_tot = []
fH1_tot = []
fH2_tot = []
file_etot = open('etot.dat', 'w')
file_fO_tot = open('fO_total.dat', 'w')
file_fH1_tot = open('fH1_total.dat', 'w')
file_fH2_tot = open('fH2_total.dat', 'w')

ebs = []
fO_ebs = []
fH1_ebs = []
fH2_ebs = []
file_ebs = open('ebs.dat', 'w')
file_fO_ebs = open('fO_ebs.dat', 'w')
file_fH1_ebs = open('fH1_ebs.dat', 'w')
file_fH2_ebs = open('fH2_ebs.dat', 'w')

usr = []
fO_usr = []
fH1_usr = []
fH2_usr = []
file_usr = open('usr.dat', 'w')
file_fO_usr = open('fO_usr.dat', 'w')
file_fH1_usr = open('fH1_usr.dat', 'w')
file_fH2_usr = open('fH2_usr.dat', 'w')

uxcdcc =[]
file_uxcdcc = open('uxcdcc.dat', 'w')

components = []

fO_pulay = []
fH1_pulay = []
fH2_pulay = []
file_fO_pulay = open('fO_pulay.dat', 'w')
file_fH1_pulay = open('fH1_pulay.dat', 'w')
file_fH2_pulay = open('fH2_pulay.dat', 'w')

kinetic = []
fO_kinetic = []
fH1_kinetic = []
fH2_kinetic = []
file_kinetic = open('kinetic.dat', 'w')
file_fO_kinetic = open('fO_kinetic.dat', 'w')
file_fH1_kinetic = open('fH1_kinetic.dat', 'w')
file_fH2_kinetic = open('fH2_kinetic.dat', 'w')

vna = []
fO_vna_2c = []
fH1_vna_2c = []
fH2_vna_2c = []
fO_vna_3c = []
fH1_vna_3c = []
fH2_vna_3c = []
file_vna = open('vna.dat', 'w')
file_fO_vna_2c = open('fO_vna_2c.dat', 'w')
file_fH1_vna_2c = open('fH1_vna_2c.dat', 'w')
file_fH2_vna_2c = open('fH2_vna_2c.dat', 'w')
file_fO_vna_3c = open('fO_vna_3c.dat', 'w')
file_fH1_vna_3c = open('fH1_vna_3c.dat', 'w')
file_fH2_vna_3c = open('fH2_vna_3c.dat', 'w')

file_VXC = open('vxc.dat', 'w')
file_VNL = open('vnl.dat', 'w')
file_EWDSR = open('ewaldsr.dat', 'w')
file_EWDLR = open('ewaldlr.dat', 'w')

for n in range(201, 300):
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
                file_fO_ebs.write('%d %f\n' % (n, float(data[3])))
                fO_ebs.append(data[3])   
            if re.search('f_ebs    2', line):
                data = line.split()
                file_fH1_ebs.write('%d %f\n' % (n, float(data[3])))
                fH1_ebs.append(data[3])
            if re.search('f_ebs    3', line):
                data = line.split()
                file_fH2_ebs.write('%d %f\n' % (n, float(data[3])))
                fH2_ebs.append(data[3])                     
                
            if re.search('f_usr    1', line):
                data = line.split()
                file_fO_usr.write('%d %f\n' % (n, float(data[3])))
                fO_usr.append(data[3])
            if re.search('f_usr    2', line):
                data = line.split()
                file_fH1_usr.write('%d %f\n' % (n, float(data[3])))
                fH1_usr.append(data[3])
            if re.search('f_usr    3', line):
                data = line.split()
                file_fH2_usr.write('%d %f\n' % (n, float(data[3])))
                fH2_usr.append(data[3])
                
            if re.search('f_pulay    1', line):
               data = line.split()
               file_fO_pulay.write('%d %f\n' % (n, float(data[3])))
               fO_pulay.append(data[3])
            if re.search('f_pulay    2', line):
               data = line.split()
               file_fH1_pulay.write('%d %f\n' % (n, float(data[3])))
               fH1_pulay.append(data[3])                
            if re.search('f_pulay    3', line):
               data = line.split()
               file_fH2_pulay.write('%d %f\n' % (n, float(data[3])))
               fH2_pulay.append(data[3])
                
            if re.search('f_kinetic    1', line):
               data = line.split()
               file_fO_kinetic.write('%d %f\n' % (n, float(data[3])))
               fO_kinetic.append(data[3])
            if re.search('f_kinetic    2', line):
               data = line.split()
               file_fH1_kinetic.write('%d %f\n' % (n, float(data[3])))
               fH1_kinetic.append(data[3])                
            if re.search('f_kinetic    3', line):
               data = line.split()
               file_fH2_kinetic.write('%d %f\n' % (n, float(data[3])))
               fH2_kinetic.append(data[3])

            if re.search('f_vna_2c    1', line):
               data = line.split()
               file_fO_vna_2c.write('%d %f\n' % (n, float(data[3])))
               fO_vna_2c.append(data[3])
            if re.search('f_vna_2c    2', line):
               data = line.split()
               file_fH1_vna_2c.write('%d %f\n' % (n, float(data[3])))
               fH1_vna_2c.append(data[3])                
            if re.search('f_vna_2c    3', line):
               data = line.split()
               file_fH2_vna_2c.write('%d %f\n' % (n, float(data[3])))
               fH2_vna_2c.append(data[3])
            if re.search('f_vna_3c    1', line):
               data = line.split()
               file_fO_vna_3c.write('%d %f\n' % (n, float(data[3])))
               fO_vna_3c.append(data[3])
            if re.search('f_vna_3c    2', line):
               data = line.split()
               file_fH1_vna_3c.write('%d %f\n' % (n, float(data[3])))
               fH1_vna_3c.append(data[3])                
            if re.search('f_vna_3c    3', line):
               data = line.split()
               file_fH2_vna_3c.write('%d %f\n' % (n, float(data[3])))
               fH2_vna_3c.append(data[3])               
                       
            if re.search('f_total =     1', line):
               data = line.split()
               file_fO_tot.write('%d %f\n' % (n, float(data[3])))
               fO_tot.append(data[3])
            if re.search('f_total =     2', line):
               data = line.split()
               file_fH1_tot.write('%d %f\n' % (n, float(data[3])))
               fH1_tot.append(data[3])
            if re.search('f_total =     3', line):
               data = line.split()
               file_fH2_tot.write('%d %f\n' % (n, float(data[3])))
               fH2_tot.append(data[3])                  

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
        file_kinetic.write('%d %f\n' % (n, float(components[1]))) 
        
        vna.append(components[2])
        file_vna.write('%d %f\n' % (n, float(components[2])))
        
        file_VXC.write('%d %f\n' % (n, float(components[3])))
        file_VNL.write('%d %f\n' % (n, float(components[4])))
        file_EWDSR.write('%d %f\n' % (n, float(components[5])))
        file_EWDLR.write('%d %f\n' % (n, float(components[6])))
# end reading files        
        
# Now calculate the deltaE contributions and compare these with the forces. 
print(' ')
print(f"delta_KE/deltaR{' '*10}force_KE")
for i in range(1,len(kinetic)):
        delta_KE = (float(kinetic[i]) - float(kinetic[i-1]))/0.01
        force_KE = 0.25*(float(fO_kinetic[i]) + float(fH2_kinetic[i]) - float(fH1_kinetic[i]) \
                         + float(fO_kinetic[i-1]) + float(fH2_kinetic[i-1]) - float(fH1_kinetic[i-1]))
        print("{:.4f}".format(delta_KE), ' '*15, "{:.4f}".format(force_KE)) 

print(' ')      
print(f"delta_VNA/deltaR{' '*10}force_VNA")
for i in range(1,len(vna)):
        delta_VNA = (float(vna[i]) - float(vna[i-1]))/0.01
        force_VNA = 0.25*(float(fO_vna_2c[i]) + float(fH2_vna_2c[i]) - float(fH1_vna_2c[i]) \
                          + float(fO_vna_3c[i]) + float(fH2_vna_3c[i]) - float(fH1_vna_3c[i]) \
                          + float(fO_vna_2c[i-1]) + float(fH2_vna_2c[i-1]) - float(fH1_vna_2c[i-1]) \
                          + float(fO_vna_3c[i-1]) + float(fH2_vna_3c[i-1]) - float(fH1_vna_3c[i-1]))
        print("{:.4f}".format(delta_VNA), ' '*15, "{:.4f}".format(force_VNA)) 

print(' ')
print(f"delta_EBS/deltaR{' '*10}force_EBS")
for i in range(1,len(ebs)):
        delta_EBS = (float(ebs[i]) - float(ebs[i-1]))/0.01       
        force_EBS = 0.25*(float(fO_ebs[i]) + float(fH2_ebs[i]) - float(fH1_ebs[i]) \
                          + float(fO_ebs[i-1]) + float(fH2_ebs[i-1]) - float(fH1_ebs[i-1]))
        force_Pulay = 0.25*(float(fO_pulay[i]) + float(fH2_pulay[i]) - float(fH1_pulay[i]) \
                          + float(fO_pulay[i-1]) + float(fH2_pulay[i-1]) - float(fH1_pulay[i-1]))
        print("{:.4f}".format(delta_EBS), ' '*15, "{:.4f}".format(force_EBS), ' '*15, "{:.4f}".format(force_Pulay))
        
print(' ')      
print(f"delta_USR/deltaR{' '*10}force_USR")
for i in range(1,len(usr)):
        delta_USR = (float(usr[i]) + float(uxcdcc[i]) \
                     - float(usr[i-1]) - float(uxcdcc[i-1]))/0.01
        force_USR = 0.25*(float(fO_usr[i]) + float(fH2_usr[i]) - float(fH1_usr[i]) \
                          + float(fO_usr[i-1]) + float(fH2_usr[i-1]) - float(fH1_usr[i-1]))
        print("{:.4f}".format(delta_USR), ' '*15, "{:.4f}".format(force_USR))
       
print(' ')
print(f"delta_ETOT/deltaR{' '*10}force_TOT")
for i in range(1,len(etot)):
        delta_ETOT = (float(etot[i]) - float(etot[i-1]))/0.01
        force_TOT = 0.25*(float(fO_tot[i]) + float(fH2_tot[i]) - float(fH1_tot[i]) \
                          + float(fO_tot[i-1]) + float(fH2_tot[i-1]) - float(fH1_tot[i-1]))
        force_Pulay = 0.25*(float(fO_pulay[i]) + float(fH2_pulay[i]) - float(fH1_pulay[i]) \
                          + float(fO_pulay[i-1]) + float(fH2_pulay[i-1]) - float(fH1_pulay[i-1]))
        print("{:.4f}".format(delta_ETOT), ' '*15, "{:.4f}".format(force_TOT), ' '*15, "{:.4f}".format(force_Pulay))        