path_in= 'dmatxvxc.txt'
docu_in= open(path_in)
path_out= 'pipi.txt'
docu_out= open(path_out, 'w')
n=1
for line in docu_in:
        docu_out.write('%d %s\n' % (n, line))
        n = n+ 1
