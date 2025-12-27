from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

#Define the list of interface residues on each chain 
def interacting_residues(pdbqt_file, threshold, output_file):
    '''
    Docstring for interacting_residues
    
    :param pdbqt_file: path to the pdbqt file to get the sequence
    :param threshold: distance we want as cutoff 
    :param output_file: file where we want to save the interface residues
    '''

    #Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdbqt_file)
    model = structure[0] #getting the model of the structure

    #Keep atoms for each chain
    atoms_per_chain = {}
    for chain in model:
        atoms = list(chain.get_atoms())
        if atoms:
            atoms_per_chain[chain.id] = atoms

    #Identify interface residues
    interface_residues = {}                   
    for chain_id in atoms_per_chain:
        interface_residues[chain_id] = set() #for {chain1 : set(residue_name, residue_id, chain_of_residue), chain2: ...}

    #Execute a neighbor search for each atom
    all_atoms = [atom for atoms in atoms_per_chain.values() for atom in atoms]

    neigbour_search = NeighborSearch(all_atoms)

    for chain_id, atoms in atoms_per_chain.items():
        for atom in atoms:
            #Find all atoms within a threshold (we decided 8A)
            neighbours = neigbour_search.search(atom.coord, threshold)

            for n_atom in neighbours:
                #check to which chain the atom belongs to
                other_chain = n_atom.get_parent().get_parent().id # the first get_parent() gets the residue, the second the chain, and we keep the id

                if other_chain != chain_id:
                    res = atom.get_parent()
                    interface_residues[chain_id].add(res)
                    #interface_residues[chain_id].add((res.resname, res.id[1], other_chain))
                    # 
    

    interface_atoms = set() #collect atoms from detected residues

    for residues in interface_residues.values(): #code to keep features of each atom for the output file
        for res in residues:
            for atom in res.get_atoms():
                interface_atoms.add(atom)
          

    interface_residue_keys = set() #(chain, residue ID)

    for chain_id, residues in interface_residues.items():
        for res in residues:
            interface_residue_keys.add((chain_id, res.id[1]))


    #OUTPUT THE INTERFACE RESIDUES IN A PDBQT FILE

    #saving output in a pdbqt file
    with open(output_file, 'w') as fout, open(pdbqt_file, 'r') as fin:
        for line in fin:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                resnum = int(line[22:26])
                if (chain_id, resnum) in interface_residue_keys:
                    fout.write(line)
            elif line.startswith(("TER", "END")):
                fout.write(line)
    return 
    

