from typing import List
from fastapi import FastAPI, HTTPException, Request, UploadFile, File, Query, Path
from models import Molecule
from rdkit import Chem
from os import getenv


def substructure_search(mols, mol):
    # Convert the SMILES string to an RDKit molecule object
    molecule = Chem.MolFromSmiles(mol)
    return [smile for smile in mols if Chem.MolFromSmiles(smile).HasSubstructMatch(molecule)]


# print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))

app = FastAPI()

# Create a list to hold instances of Molecule and add instances of Molecule to the list
molecules: List[Molecule] = [Molecule(mol_id=1, name="CCO"),    # Ethanol
                             Molecule(mol_id=2, name="c1ccccc1"),   # Benzene
                             Molecule(mol_id=3, name="CC(=O)O"),    # Acetic acid
                             Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")]  # Aspirin


# This method checks balancing
@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules/{mol_id}", status_code=201)
def add_molecule(molecule: Molecule):
    # Add Molecule object into the molecules
    molecules.append(molecule)
    return molecule


@app.get("/molecules/{mol_id}", response_model=Molecule)
def retrieve_molecule(mol_id: int):
    for molecule in molecules:
        if molecule.mol_id == mol_id:
            return molecule
    raise HTTPException(status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.put("/molecules/{mol_id}", response_model=Molecule)
def update_molecule(mol_id: int, updated_molecule: Molecule):
    for index, molecule in enumerate(molecules):
        if molecule.mol_id == mol_id:
            molecule.name = updated_molecule.name
            return updated_molecule
    raise HTTPException(status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.delete("/molecules/{mol_id}", response_model=Molecule)
def delete_molecule(mol_id: int):
    for index, molecule in enumerate(molecules):
        if molecule.mol_id == mol_id:
            deleted_molecule = molecules.pop(index)
            return deleted_molecule
    raise HTTPException(status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.get("/molecules")
def retrieve_molecules():
    return molecules


@app.get("/search", response_model=List[str])
def search_molecules_by_smile(substructure_smile: str = Query(..., description="SMILES string of the substructure.")):
    # Create a list of molecule names(SMILES)
    smiles = [mol.name for mol in molecules]

    # If the passed argument is unknown type or does not exist raise an exception
    substructure = Chem.MolFromSmiles(substructure_smile)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES.")

    # Otherwise return a list of molecules containing all elements that has substructure of the given Smile
    return substructure_search(smiles, substructure_smile)


@app.post("/add")
async def upload_file(file: UploadFile = File(...)):
    # Other tasks are running while waiting entire file to be read
    content = await file.read()
    # Convert the raw bytes into a string using "UTF-8" encoding
    lines = content.decode("utf-8").splitlines()

    # Iterate over each line and add fetched molecule data into molecules list
    for line in lines:
        mol_id, name = line.split()
        molecules.append(Molecule(mol_id=int(mol_id), name=name))

    # Return success message
    return {"content": "Molecule/Molecules added successfully."}
