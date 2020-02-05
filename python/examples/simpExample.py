from liggghts import liggghts

# Try loading in liggghts
lmp = liggghts()

# Try reading in an input file
lmp.file("in.simpExample")


# Try performing a command
print("")
print("Attempting to run simulation another 500 steps")
lmp.command("run 500")

# Try extracting a global value
print("")
print("Attempting to get the number of atoms in simulation")
numAtoms = lmp.extract_global("natoms", 0)
print("natoms =", numAtoms)

# Try extracting atom's positions
print("")
print("Attempting to get the atom's positions")
pos = lmp.extract_atom("x",3)
for k in range(0,numAtoms):
    print("Pos[%i] = [%f, %f, %f]" % (k, pos[k][0], pos[k][1], pos[k][2]))

# Try extracting a compute
print("")
print("Attempting to get the rotational KE")
rke = lmp.extract_compute("rke",0,0)
print("rke = %f" % rke)

# Try extracting a fix
print("")
print("Attempting to get the pos_z before python ran another 500 steps")
z0 = lmp.extract_fix("z0",1,1)
print("z0[0] = %f" % z0[0])

# Try extracting a variable
print("")
print("Attempting to get the density variable used in simulation")
density = lmp.extract_variable("density","all",0)
print("density = %f" % density)

# Testing get_natoms
print("")
print("Testing get_natoms")
natoms = lmp.get_natoms()
print("natoms = %i" % natoms)

# Testing gather_atoms
print("")
print("Testing gather_atoms (NOTICE THAT IDS ARE DIFFERENT!!!)")
x_gather_atoms = lmp.gather_atoms("x",1,3)
for k in range(natoms):
    x = x_gather_atoms[3*k+0]
    y = x_gather_atoms[3*k+1]
    z = x_gather_atoms[3*k+2]
    print("x[%i] = [%f, %f, %f]" % (k, x, y, z))
