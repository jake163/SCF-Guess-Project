import os

def give_geometry(filename):
    cwd = os.getcwd()
    coords = []
    if filename.endswith(".xyz"):
        with open(os.path.join(cwd, "geometries", filename)) as f:
            num_atoms = 0
            counter = 0
            for line in f:
                counter += 1
                if counter == 0:
                    num_atoms = int(line.split()[0])
                elif counter == 2:
                    continue
                elif counter > 2:
                    element = line.split()
                    symbol = element[0]
                    x,y,z = float(element[1]), float(element[2]), float(element[3])
                    coords.append([symbol, x, y, z])
                elif counter == 2+num_atoms:
                      break
    if filename.endswith(".sdf"):
        with open(os.path.join(cwd, "geometries", filename)) as f:
            counter = 0
            natoms = 0
            for line in f:
                if counter == 3:
                    natoms = int(line.split()[0])
                elif counter == 4 + natoms:
                    break
                elif counter > 3:
                    atom = line.split()
                    at_sym = atom[3]
                    x,y,z = float(atom[0]), float(atom[1]), float(atom[2])
                    coords.append([at_sym, x, y, z])
                counter += 1
    return coords
