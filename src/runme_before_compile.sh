#!/bin/bash

MASTER_PATH=fireball-master

for i in a.GLOBAL b.FUNCTIONS c.SYSTEM d.FUNCTIONS_EXTRA e.FDATA g.XC_FUNCTIONALS h.SOLVESH j.ASSEMBLERS k.DASSEMBLERS l.SCF m.MD n.NAMD p.THEORY o.OUTPUT u.UTIL include libs Makefile.in MACHINES
do
    if [ -h $i ]
    then
	rm $i
    fi
    ln -s ../../${MASTER_PATH}/src/$i
done
