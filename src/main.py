from typing import List
from fastapi import FastAPI, HTTPException, UploadFile, File, Query
from models import Molecule
from rdkit import Chem
from os import getenv


def substructure_search(mols, mol):
    # Convert the SMILES string to an RDKit molecule object
    molecule = Chem.MolFromSmiles(mol)
    if not molecule:
        return []
    return [smile for smile in mols
            if Chem.MolFromSmiles(smile).HasSubstructMatch(molecule)]


app = FastAPI()

# Create a list to hold instances of Molecule and add instances
# of Molecule to the list
molecules: List[Molecule] = [
    Molecule(mol_id=1, name="CCO"),  # Ethanol
    Molecule(mol_id=2, name="c1ccccc1"),  # Benzene
    Molecule(mol_id=3, name="CC(=O)O"),  # Acetic acid
    Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
]


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules/{mol_id}", status_code=201)
def add_molecule(molecule: Molecule):
    # Add Molecule object into the molecules
    print("POST request received at /molecules")
    molecules.append(molecule)
    return molecule


@app.get("/molecules/{mol_id}", response_model=Molecule)
def retrieve_molecule(mol_id: int):
    for molecule in molecules:
        if molecule.mol_id == mol_id:
            return molecule
    raise HTTPException(
        status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.put("/molecules/{mol_id}", response_model=Molecule)
def update_molecule(mol_id: int, updated_molecule: Molecule):
    for index, molecule in enumerate(molecules):
        if molecule.mol_id == mol_id:
            molecule.name = updated_molecule.name
            return updated_molecule
    raise HTTPException(
        status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.delete("/molecules/{mol_id}", response_model=Molecule, status_code=200)
def delete_molecule(mol_id: int):
    for index, molecule in enumerate(molecules):
        if molecule.mol_id == mol_id:
            deleted_molecule = molecules.pop(index)
            print("Deleted molecule:", deleted_molecule)
            print("Remaining molecules:", molecules)
            return deleted_molecule
    raise HTTPException(
        status_code=404, detail=f'Molecule with id {mol_id} not found.')


@app.get("/molecules")
def retrieve_molecules():
    return molecules


@app.get("/search", response_model=List[str])
def search_molecules_by_smile(
        substructure_smile:
        str =
        Query(..., description="SMILES string of the substructure.")):
    smiles = [mol.name for mol in molecules]

    substructure = Chem.MolFromSmiles(substructure_smile)
    if not substructure:
        raise HTTPException(
            status_code=400, detail="Invalid substructure SMILES.")

    return substructure_search(smiles, substructure_smile)


@app.post("/add")
async def upload_file(file: UploadFile = File(...)):
    content = await file.read()
    lines = content.decode("utf-8").splitlines()

    global molecules
    existing_ids = {mol.mol_id for mol in molecules}

    for line in lines:
        parts = line.split()
        if len(parts) != 2:
            return {"detail": "Invalid file format"}

        mol_id, name = parts
        mol_id = int(mol_id)

        if mol_id not in existing_ids:
            molecules.append(Molecule(mol_id=mol_id, name=name))
            existing_ids.add(mol_id)

    return {"content": "Molecule/Molecules added successfully."}


@app.post("/clear_molecules", status_code=200)
def clear_molecules():
    global molecules
    molecules = []  # Clear the list
    return {"detail": "Molecules cleared successfully"}
