import numpy as np
#import matplotlib.pyplot as plt

pathforce1= 'force1.txt'
docuforce1= np.loadtxt(pathforce1)
print docuforce1[:,0]
pathforce2= 'force2.txt'
docuforce2= np.loadtxt(pathforce2)
pathforce3= 'force3.txt'
docuforce3= np.loadtxt(pathforce3)
pathforcevna= 'forcevna.txt'
docuforcevna= np.loadtxt(pathforcevna)

pathforcevnl= 'forcevnl.txt'
docuforcevnl= np.loadtxt(pathforcevnl)

pathforceusr= 'forceusr.txt'
docuforceusr= np.loadtxt(pathforceusr)

pathforcevxc= 'forcevxc.txt'
docuforcevxc= np.loadtxt(pathforcevxc)


#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ff = np.loadtxt('rho.txt')
#x =ff[:,0]
#y =ff[:,1]
#z =ff[:,2]



#ax.scatter(x, y, z, c='r', marker='.')

#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')

#plt.savefig('rho-1.png')
#plt.show()
