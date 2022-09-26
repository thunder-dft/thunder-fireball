pathebs= 'ebs.txt'
docuebs= open(pathebs, 'w')

pathke= 'ke.txt'
docuke= open(pathke, 'w')

pathvna= 'vna.txt'
docuvna= open(pathvna, 'w')

pathvxc= 'vxc.txt'
docuvxc= open(pathvxc, 'w')

pathvnl= 'vnl.txt'
docuvnl= open(pathvnl, 'w')

for n in range(1, 200):
        if n < 10:
            pathSour1= '00%d.ENERGIES' % n
        if n >= 10 and n < 100:
            pathSour1= '0%d.ENERGIES' % n
        if n >= 100:
            pathSour1= '%d.ENERGIES' % n
        docuSour1 = open(pathSour1)
        l =0
        for line in docuSour1:
            l = l+1
        m =0
        docuSour1 = open(pathSour1)
        for line in docuSour1:
            m = m+1
            if m == l:
                data= line.split()
                docuebs.write('%d  %s\n' % (n, data[0]))
                docuke.write('%d  %s\n' % (n, data[1]))
                docuvna.write('%d  %s\n' % (n, data[2]))
                docuvxc.write('%d  %s\n' % (n, data[3]))
                docuvnl.write('%d  %s\n' % (n, data[4]))
                
