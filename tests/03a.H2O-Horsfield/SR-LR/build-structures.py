#!/usr/bin/env python

#!initialize some numbers for starting position
iH1 = 1
xH1 =  0.61
yH1 =  0.00
zH1 =  0.00
iH2 = 1
xH2 =  0.00
yH2 =  1.00
zH2 =  0.00
iO2 = 8
xO2 =  0.00
yO2 = -4.00
zO2 =  0.00
iH3 = 1
xH3 =  1.00
yH3 = -4.00
zH3 =  0.00
iH4 = 1
xH4 = -0.242867
yH4 = -4.939103
zH4 =  0.00

for n in range(1, 101):
        # write out energies to .txt files from *.log files 
        if n < 10:
            source = '00%d.inp' % n
        if n >= 10 and n < 100:
            source = '0%d.inp' % n
        if n >= 100:
            source = '%d.inp' % n
           
        file_in = open('000.inp', 'r')
        file_out = open(source, 'wt')

        i = 0
        while i < 7:
            line = file_in.readline()
            file_out.write(line)
            i = i + 1
        print(' '*2, "{}".format(iH1), ' '*2, "{:.6f}".format(xH1), ' ', \
                      "{:.6f}".format(yH1), ' '*2, "{:.6f}".format(zH1), file = file_out)
        print(' '*2, "{}".format(iH2), ' '*2, "{:.6f}".format(xH2), ' ', \
                      "{:.6f}".format(yH2), ' '*2, "{:.6f}".format(zH2), file = file_out)        
        print(' '*2, "{}".format(iO2), ' '*2, "{:.6f}".format(xO2), ' ', \
                      "{:.6f}".format(yO2), ' '*2, "{:.6f}".format(zO2), file = file_out) 
        print(' '*2, "{}".format(iH3), ' '*2, "{:.6f}".format(xH3), ' ', \
                      "{:.6f}".format(yH3), ' '*2, "{:.6f}".format(zH3), file = file_out)
        print(' '*2, "{}".format(iH3), ' ', "{:.6f}".format(xH4), ' ', \
                      "{:.6f}".format(yH4), ' '*2, "{:.6f}".format(zH4), file = file_out)        
        xH1 = xH1 + 0.01
        file_in.close()
        file_out.close()
