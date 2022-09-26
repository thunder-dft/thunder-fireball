import linecache
import re
import math
pathetot= 'etot.txt'
docuetot= open(pathetot, 'w')
pathebs= 'ebs.txt'
docuebs= open(pathebs, 'w')
pathuii= 'uii-uee.txt'
docuuii= open(pathuii, 'w')
pathuxcdcc= 'uxcdcc.txt'
docuuxcdcc= open(pathuxcdcc, 'w')
pathforce1= 'force1.txt'
docuforce1= open(pathforce1, 'w')
pathforce2= 'force2.txt'
docuforce2= open(pathforce2, 'w')
pathforce3= 'force3.txt'
docuforce3= open(pathforce3, 'w')
pathforce3= 'force3.txt'
docuforce3= open(pathforce3, 'w')

pathforcek= 'forcek.txt'
docuforcek= open(pathforcek, 'w')

pathforcevna= 'forcevna.txt'
docuforcevna= open(pathforcevna, 'w')

pathforcevnl= 'forcevnl.txt'
docuforcevnl= open(pathforcevnl, 'w')

pathforceusr= 'forceusr.txt'
docuforceusr= open(pathforceusr, 'w')

pathforcevxc= 'forcevxc.txt'
docuforcevxc= open(pathforcevxc, 'w')

force2= 0.000
force3= 0.000
for n in range(1, 200):
	if n < 10:
        	pathSour= '00%d.log' % n
	if n >= 10 and n < 100:
		pathSour= '0%d.log' % n
	if n >= 100:
	        pathSour= '%d.log' % n

        docuSour= open(pathSour)
        l = 0
        for line in docuSour:
                l = l+1
                line= line.rstrip()
                if re.search('ETOT', line):
                    data= line.split()
                    etot = float(data[2])
                        #docuetot.write('%d %f\n' % (n, float(data[2])))
                
#                if re.search('ebs', line):
#                    data= line.split()
                        #docuebs.write('%d %f\n' % (n, float(data[2])))
#                if re.search('uii', line):
#                    data= line.split()
                        #docuuii.write('%d %f\n' % (n, float(data[4])))
#                if re.search('uxcdcc', line):
#                    data= line.split()
#                        #docuuxcdcc.write('%d %f\n' % (n, float(data[2])))
#                if re.search('kinetic', line):
#                    newline = linecache.getline(pathSour, l +5)
#                    data= newline.split()
#                    print 'k'
#                    print data
#                    docuforcek.write('%d  %f\n' % (n,math.sqrt(float(data[2])**2+ float(data[3])**2 + float(data[4])**2)))
#                if re.search('(vna)', line):
#                    newline = linecache.getline(pathSour, l +5)
#                    data= newline.split()
#                    print 'vna'
#                    print data
#                    docuforcevna.write('%d  %f\n' % (n,math.sqrt(float(data[2])**2+ float(data[3])**2 +\
#                                    float(data[4])**2)))                	
                        
#               if re.search('(vnl)', line):
#                    newline = linecache.getline(pathSour, l +5)
#                    data= newline.split()
#                    print 'vnl'
#                    print data
#                    docuforcevnl.write('%d  %f\n' % (n,math.sqrt(float(data[2])**2+ float(data[3])**2 + float(data[4])**2)))
#                if re.search('(usr_)', line):
#		    newline = linecache.getline(pathSour, l +5)
#		    data= newline.split()
#		    print 'usr'
#		    print data
#		    docuforceusr.write('%d  %f\n' % (n,math.sqrt(float(data[2])**2+ float(data[3])**2 + float(data[4])**2)))
#		if re.search('(vxc_)', line):
#		    newline = linecache.getline(pathSour, l+5)
#	            data= newline.split()
#		    print 'vxc'
#		    print data
#		    docuforcevxc.write('%d  %f\n' % (n,math.sqrt(float(data[2])**2+ float(data[3])**2 + float(data[4])**2)))
		if re.search('ftot =', line):
		    data= line.split()
		    if int( data[2]) == 1:
		        docuforce1.write('%d  %f\n' % (n,math.sqrt(float(data[3])**2+ float(data[4])**2 + float(data[5])**2)))
                    if int( data[2]) == 2:
                    	docuforce2.write('%d  %f\n' % (n,math.sqrt(float(data[3])**2+ float(data[4])**2 + float(data[5])**2)))
                    if int( data[2]) == 3:
                    	docuforce3.write('%d  %f\n' % (n,math.sqrt(float(data[3])**2+ float(data[4])**2 + float(data[5])**2)))
	docuetot.write('%d %f\n' % (n, etot))
