data = open('rho.txt')
data_out= open('rho+.txt', 'w')
for i, line in enumerate(data):
        if i > 5000000 and i < 10000000:
                data_out.write('%s\n' % ( line))
        if i > 10000000:
                break
