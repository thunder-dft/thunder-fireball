#!/usr/bin/env python

import re
import numpy as np
import subprocess
import json
import os

pathetot= 'etot.txt'
docuetot= open(pathetot, 'w')
pathebs= 'ebs.txt'
docuebs= open(pathebs, 'w')
pathuii= 'uii-uee.txt'
docuuii= open(pathuii, 'w')
pathuxcdcc= 'uxcdcc.txt'
docuuxcdcc= open(pathuxcdcc, 'w')

for n in range(1, 200):
        if n < 10:
                pathSour= '00%d.log' % n
        if n >= 10 and n < 100:
                pathSour= '0%d.log' % n
        if n >= 100:
                pathSour= '%d.log' % n

        docuSour= open(pathSour)
        for line in docuSour:
                line= line.rstrip()
                if re.search('ETOT', line):
                        data= line.split()
                        docuetot.write('%d %f\n' % (n, float(data[2])))
                if re.search('ebs', line):
                        data= line.split()
                        docuebs.write('%d %f\n' % (n, float(data[2])))
                if re.search('uii', line):
                        data= line.split()
                        docuuii.write('%d %f\n' % (n, float(data[4])))
                if re.search('uxcdcc', line):
                        data= line.split()
                        docuuxcdcc.write('%d %f\n' % (n, float(data[2])))

def tester():
        rf=open('reference.json')
        ref=json.load(rf)
        ret = True

        for i in [70, 71, 72, 73]:
                data=np.loadtxt('fort.%s' % i)
                print 'min: %7.3f max: %7.3f avg: %7.3f' % (np.min(data[:,1]), np.max(data[:,1]), np.average(data[:,1]))
                if np.abs(np.average(data[:,1])-ref[str(i)]['avg'])>0.1:
                        ret=False
        return ret

def clean():
        files = os.listdir('.')
        for i in sorted(files):
                if i[-3:]=='inp':
                        pass
                elif i in [ 'fireball.x', 'Fdata', 'tester.py', 'reference.json' ]:
                        pass
                else:
                        os.remove(i)

def run():
        sp=subprocess.call(['./fireball.x'])

if __name__ == "__main__":
        run()
        ret=tester()
        if ret:
                print 'Test success!'
                clean()
        
