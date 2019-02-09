import glob
import os
import numpy as np

from .utils import Atom, Residue, ActiveSite


def read_active_sites(dir):
    """
    Read in all of the active sites from the given directory.

    Input: directory
    Output: list of ActiveSite instances
    """
    files = glob.glob(dir + '/*.pdb')

    active_sites = []
    # iterate over each .pdb file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.pdb")):

        active_sites.append(read_active_site(filepath))

    print("Read in %d active sites"%len(active_sites))

    return active_sites

def generate_distance_types():
        groups = ["GROUP-0","GROUP-1","GROUP-2","GROUP-3","GROUP-4"]
        atom_ids = ["CA","CB","CENTROID"]
        distance_types = []
        for group_1 in groups:
            for atom_1 in atom_ids:
                for group_2 in groups:
                    for atom_2 in atom_ids:
                        if group_1 != group_2:
                            group_ID = max(group_1,group_2)+min(group_1,group_2)
                        else:
                            group_ID = group_1+group_2
                        if atom_1 != atom_2:
                            atom_ID = max(atom_1,atom_2)+min(atom_1,atom_2)
                        else:
                            atom_ID = atom_1+atom_2
                        entry_ID = [group_ID + atom_ID]
                        distance_types.append(entry_ID)
        distance_types = np.unique(distance_types)
        return distance_types

def entry_id(group_1,group_2,atom_1,atom_2):
        if group_1 != group_2:
                group_ID = max(group_1,group_2)+min(group_1,group_2)
        else:
            group_ID = group_1+group_2
        if atom_1 != atom_2:
            atom_ID = max(atom_1,atom_2)+min(atom_1,atom_2)
        else:
            atom_ID = atom_1+atom_2
        entry_ID = group_ID + atom_ID
        return entry_ID
def distance(a,b):
        x_dist = (a[0]-b[0])**2
        y_dist = (a[1]-b[1])**2
        z_dist = (a[2]-b[2])**2
        dist = (x_dist + y_dist + z_dist)**0.5
        return dist
def coordinates(residue):
        if residue.type == 'GLY':
                            coordinates = residue.atoms[1].coords
        elif residue.type == 'ALA':
            coordinates = residue.atoms[4].coords
        else:
            coordinates = (0,0,0)
            #Calculate centroid coordinates
            #Add up all the coordinates
            for side_chain in residue.atoms[5:]:
                a = coordinates
                b = side_chain.coords
                coordinates = tuple([item1 + item2 for item1, item2 in zip(a,b)])
            #Divide by the number of side chain atoms to get the centroid
            n = len(residue.atoms[5:])
            coordinates = (coordinates[0]/n, coordinates[1]/n, coordinates[2]/n)
        return coordinates
def read_active_site(filepath):
        """
        Read in a single active site given a PDB file

        Input: PDB file path
        Output: ActiveSite instance
        
        """
        #Define atom groupings
        group_dictionary = {}
        group_dictionary['ALA'] = 'GROUP-0'
        group_dictionary['VAL'] = 'GROUP-0'
        group_dictionary['ILE'] = 'GROUP-0'
        group_dictionary['LEU'] = 'GROUP-0'
        group_dictionary['GLY'] = 'GROUP-0'
        group_dictionary['PRO'] = 'GROUP-0'
        group_dictionary['LYS'] = 'GROUP-1'
        group_dictionary['ARG'] = 'GROUP-1'
        group_dictionary['HIS'] = 'GROUP-1'
        group_dictionary['ASP'] = 'GROUP-2'
        group_dictionary['GLU'] = 'GROUP-2'
        group_dictionary['GLN'] = 'GROUP-2'
        group_dictionary['ASN'] = 'GROUP-2'
        group_dictionary['TYR'] = 'GROUP-3'
        group_dictionary['PHE'] = 'GROUP-3'
        group_dictionary['TRP'] = 'GROUP-3'
        group_dictionary['CYS'] = 'GROUP-4'
        group_dictionary['SER'] = 'GROUP-4'
        group_dictionary['THR'] = 'GROUP-4'
        #Generate list of unique distance types
        distance_types = generate_distance_types()
        
        basename = os.path.basename(filepath)
        name = os.path.splitext(basename)

        if name[1] != ".pdb":
            raise IOError("%s is not a PDB file"%filepath)

        active_site = ActiveSite(name[0])

        r_num = 0
        
        # open pdb file
        with open(filepath, "r") as f:
            # iterate over each line in the file
            for line in f:
                if line[0:3] != 'TER':
                    # read in an atom
                    atom_type = line[13:17].strip()
                    x_coord = float(line[30:38])
                    y_coord = float(line[38:46])
                    z_coord = float(line[46:54])
                    atom = Atom(atom_type)
                    atom.coords = (x_coord, y_coord, z_coord)

                    residue_type = line[17:20]
                    residue_number = int(line[23:26])
                    residue_group = group_dictionary[residue_type]
                    # make a new residue if needed
                    if residue_number != r_num:
                        residue = Residue(residue_type, residue_number,residue_group)
                        r_num = residue_number

                    # add the atom to the residue
                    residue.atoms.append(atom)

                else:  # I've reached a TER card
                    #Add coordinates for centroid
                    #Since GLY has no side chain, make the centroid coords CA
                    #check if residue has two conformations and keep just conformation A
                    if residue.atoms[0].type == "N  A":
                        only_odd = np.array([num for num in range(len(residue.atoms)) if num % 2 == 0])
                        atoms_to_keep = []
                        for index in only_odd:
                            atoms_to_keep.append(residue.atoms[index])
                        residue.atoms = atoms_to_keep
                        #residue.atoms = residue.atoms[only_odd] 
                        for atom in residue.atoms:
                            atom.type = atom.type.replace(" ", "")
                            atom.type = atom.type[:-1]
                    atom_type = 'CENTROID'
                    atom = Atom(atom_type)
                    atom.coords = coordinates(residue)
                    residue.atoms.append(atom)
                    active_site.residues.append(residue) 
        #Store distances between all relevant atoms in each residue
        #Initialize distance dictionary
        active_site.distances = {}
        for distance_type in distance_types:
            active_site.distances[distance_type] = []
        index = 0
        for residue_a in active_site.residues:
            if residue_a.type == "GLY":
                atoms_1 = [residue_a.atoms[1],residue_a.atoms[-1],residue_a.atoms[-1]]
            else:
                atoms_1 = [residue_a.atoms[1],residue_a.atoms[4],residue_a.atoms[-1]]
            group_1 = residue_a.group
            for residue_b in active_site.residues[(index+1):]:
                if residue_b.type == "GLY":
                    atoms_2 = [residue_b.atoms[1],residue_b.atoms[-1],residue_b.atoms[-1]]
                else:
                    atoms_2 = [residue_b.atoms[1],residue_b.atoms[4],residue_b.atoms[-1]]
                group_2 = residue_b.group
                for atom_1 in atoms_1:
                    for atom_2 in atoms_2:
                        atom_1_type = atom_1.type
                        atom_2_type = atom_2.type
                        pairing_id = entry_id(group_1,group_2,atom_1_type,atom_2_type)
                        pairing_dist = [distance(atom_1.coords,atom_2.coords)]
                        active_site.distances[pairing_id] = active_site.distances[pairing_id] + pairing_dist
            index = index+1
        #Sort each distance list in ascendning order
        for item in active_site.distances.keys():
            active_site.distances[item] = sorted(active_site.distances[item])
        return active_site



def write_clustering(filename, clusters):
    """
    Write the clustered ActiveSite instances out to a file.

    Input: a filename and a clustering of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusters)):
        out.write("\nCluster %d\n--------------\n" % i)
        for j in range(len(clusters[i])):
            out.write("%s\n" % clusters[i][j])

    out.close()


def write_mult_clusterings(filename, clusterings):
    """
    Write a series of clusterings of ActiveSite instances out to a file.

    Input: a filename and a list of clusterings of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusterings)):
        clusters = clusterings[i]

        for j in range(len(clusters)):
            out.write("\nCluster %d\n------------\n" % j)
            for k in range(len(clusters[j])):
                out.write("%s\n" % clusters[j][k])

    out.close()
