inp_file='force3.txt'
inp_data= open(inp_file)
out_file='force3_r.txt'
out_file= open(out_file, 'w') 
n = 1.00
for line in inp_data:
	data= line.split()
	out_file.write('%f %s\n' % (0.800 + n*0.003, data[1])) 
	n = n + 1
