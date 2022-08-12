./clean.com
cp structures.forces-HARRIS.inp structures.inp
./fireball.x
python test-forces-HARRIS.py > test-forces.HARRIS.dat
cp structures.forces-DOGS.inp structures.inp
./fireball.x
python test-forces-DOGS.py > test-forces.DOGS.dat
