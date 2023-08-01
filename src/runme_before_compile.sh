#!/bin/bash
MASTER_PATH=thunder-master

for i in include Makefile MACHINES a.GLOBAL b.FUNCTIONS c.SYSTEM e.FDATA f.MPI g.XC_FUNCTIONALS h.SOLVESH i.GRID j.ASSEMBLERS k.DASSEMBLERS l.SCF m.MD o.OUTPUT u.UTIL
do
    if [ -h $i ]
    then
	rm $i
    fi
    ln -s ../../${MASTER_PATH}/src/$i
done
