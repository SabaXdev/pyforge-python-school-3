from rdkit import Chem


def substructure_search(mols, mol):
    # Convert the SMILES string to an RDKit molecule object
    molecule = Chem.MolFromSmiles(mol)
    return [smile for smile in mols if Chem.MolFromSmiles(smile).HasSubstructMatch(molecule)]
