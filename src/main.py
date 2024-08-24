from typing import List
from fastapi import FastAPI, HTTPException, UploadFile, File, Query, Depends
from sqlalchemy.orm import Session

import crud
import models
import schemas
from database import engine, SessionLocal
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


models.Base.metadata.create_all(bind=engine)

app = FastAPI()


# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


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
def add_molecule(molecule: schemas.MoleculeCreate,
                 db: Session = Depends(get_db)):
    # Add Molecule object into the molecules
    return crud.create_molecule(db=db, molecule=molecule)


@app.get("/molecules/{mol_id}", response_model=schemas.Molecule)
def retrieve_molecule(mol_id: int, db: Session = Depends(get_db)):
    db_molecule = crud.get_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@app.put("/molecules/{mol_id}", response_model=schemas.Molecule)
def update_molecule(mol_id: int, molecule: schemas.MoleculeCreate,
                    db: Session = Depends(get_db)):
    db_molecule = crud.update_molecule_by_id(db=db, molecule_id=mol_id,
                                             updated_molecule=molecule)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@app.delete("/molecules/{mol_id}", response_model=schemas.Molecule,
            status_code=200)
def delete_molecule(mol_id: int, db: Session = Depends(get_db)):
    db_molecule = crud.delete_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@app.get("/molecules")
def retrieve_molecules(skip: int = 0, limit: int = 100,
                       db: Session = Depends(get_db)):
    return crud.get_all_molecules(db=db, skip=skip, limit=limit)


@app.get("/search", response_model=List[str])
def search_molecules_by_smile(
        substructure_smile: str =
        Query(..., description="SMILES string of the substructure."),
        db: Session = Depends(get_db)):
    try:
        return crud.search_molecules(
            db=db, substructure_smile=substructure_smile)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/add")
async def upload_file(file: UploadFile = File(...),
                      db: Session = Depends(get_db)):
    content = await file.read()

    try:
        return crud.add_molecules_from_file(
            db=db, file_content=content.decode("utf-8"))
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/clear_molecules", status_code=200)
def clear_molecules(db: Session = Depends(get_db)):
    return crud.clear_all_molecules(db=db)
