from math import *
import numpy as np
def x_rotation (angle, vect):
    rota= [0.0, 0.0, 0.0]
    rota[1]= vect[1]*cos(angle) + vect[2]*sin(angle)
    rota[2]= vect[2]*cos(angle) - vect[1]*sin(angle)
    return rota;
#6           0.000000     0.000000     0.000000
#1           0.629118     0.629118     0.629118
#1           -0.629118    -0.629118    0.629118
#1           0.629118     -0.629118    -0.629118
#1           -0.629118    0.629118     -0.629118

angle=0.0000
r_meas = sqrt(0.629118**2 + 0.629118**2 + 0.629118**2)
r_unit2= [0.629118/r_meas, 0.629118/r_meas, 0.629118/r_meas]
r_unit3= [-0.629118/r_meas, -0.629118/r_meas, 0.629118/r_meas]
r_unit4= [0.629118/r_meas, -0.629118/r_meas, -0.629118/r_meas]
r_unit5= [-0.629118/r_meas, 0.629118/r_meas, -0.629118/r_meas]
r_atom2= np.multiply(r_unit2, 0.55)
r_atom3= np.multiply(r_unit3, 0.55)
r_atom4= np.multiply(r_unit4, 0.55)
r_atom5= np.multiply(r_unit5, 0.55)
delta_r2= np.multiply(r_unit2, 0.01)
delta_r3= np.multiply(r_unit3, 0.01)
delta_r4= np.multiply(r_unit4, 0.01)
delta_r5= np.multiply(r_unit5, 0.01)
for n in range(1, 200):
    if n < 10:
        path="00%d.inp" % n
    
    if n >= 10 and n < 100:
        path="0%d.inp" % n
    
    if n >= 100:
        path="%d.inp" % n
    
    myfile= open(path, 'w')
    myfile.write('  5            1\n')
    myfile.write("  999.00000    0.000000    0.000000\n")
    myfile.write("    0.00000  999.000000    0.000000\n")
    myfile.write("    0.00000    0.000000  999.000000\n")
    myfile.write("  1\n")
    myfile.write("  0.000  0.000  0.000   1.000\n")
    myfile.write("  6    0.000000   0.000000    0.000000\n")
    r_atom2= np.add(r_atom2,delta_r2)
    r_atom3= np.add(r_atom3,delta_r3)
    r_atom4= np.add(r_atom4,delta_r4)
    r_atom5= np.add(r_atom5,delta_r5)
    myfile.write("  1   %f   %f    %f\n" % (r_atom2[0], r_atom2[1], r_atom2[2]))
    myfile.write("  1   %f   %f    %f\n" % (r_atom3[0], r_atom3[1], r_atom3[2]))
    myfile.write("  1   %f   %f    %f\n" % (r_atom4[0], r_atom4[1], r_atom4[2]))
    myfile.write("  1   %f   %f    %f\n" % (r_atom5[0], r_atom5[1], r_atom5[2]))
    

    
              
