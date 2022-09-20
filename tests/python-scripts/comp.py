import numpy as np
elements= ['etot.txt', 'uii-uee.txt', 'uxcdcc.txt', 'force1.txt', \
           'force2.txt', 'force3.txt', 'forcek.txt', 'forcevna.txt', \
           'forcevnl.txt', 'forceusr.txt', 'forcevxc.txt', 'ebs.txt', \
           'ke.txt', 'vna.txt', 'vxc.txt', 'vnl.txt']

elemets_ref=['reference/etot_ref.txt', 'reference/uii-uee_ref.txt', \
              'reference/uxcdcc_ref.txt', 'reference/force1_ref.txt', \
           'reference/force2_ref.txt', 'reference/force3_ref.txt',\
            'reference/forcek_ref.txt', 'reference/forcevna_ref.txt', \
           'reference/forcevnl_ref.txt', 'reference/forceusr_ref.txt', \
           'reference/forcevxc_ref.txt', 'reference/ebs_ref.txt', \
           'reference/ke_ref.txt', 'reference/vna_ref.txt', \
           'reference/vxc_ref.txt', 'reference/vnl_ref.txt']
n=0
for i in elements:
    data_test= np.loadtxt(i)
    data_reff=np.loadtxt(elemets_ref[n])
    n= n+1
    valu_test= 0.0
    valu_refe= 0.0
    print shape(data_reff)
    #for j range(elements.shape):
        #valu_refe = valu_refe + data_reff[1][j]
        #valu_test = valu_test + data_test[1][j]
     #   print data_reff[1][j], data_test[1][j]