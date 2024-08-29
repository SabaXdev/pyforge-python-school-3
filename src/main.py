from typing import List, Iterator
from fastapi import FastAPI, HTTPException, UploadFile, File, Query, Depends
from sqlalchemy.orm import Session
from fastapi.middleware.cors import CORSMiddleware

import logging
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

# Logging configuration
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),  # Logs to console
        logging.FileHandler("app.log")  # Logs to a file
    ]
)

logger = logging.getLogger(__name__)


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
    logger.info("Fetching server ID")
    return {"server_id": getenv("SERVER_ID", "1")}


@app.post("/molecules/{mol_id}", status_code=201)
def add_molecule(molecule: schemas.MoleculeCreate,
                 db: Session = Depends(get_db)):
    # Add Molecule object into the molecules
    return crud.create_molecule(db=db, molecule=molecule)


@app.get("/molecules/{mol_id}", response_model=schemas.Molecule)
def retrieve_molecule(mol_id: int, db: Session = Depends(get_db)):
    logger.info(f"Retrieving molecule with ID {mol_id}")
    db_molecule = crud.get_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@app.put("/molecules/{mol_id}", response_model=schemas.Molecule)
def update_molecule(mol_id: int, molecule: schemas.MoleculeCreate,
                    db: Session = Depends(get_db)):
    logger.info(f"Updating molecule with ID {mol_id}")
    db_molecule = crud.update_molecule_by_id(db=db, molecule_id=mol_id,
                                             updated_molecule=molecule)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found for update")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"Molecule with ID {mol_id} successfully updated")
    return db_molecule


@app.delete("/molecules/{mol_id}", response_model=schemas.Molecule,
            status_code=200)
def delete_molecule(mol_id: int, db: Session = Depends(get_db)):
    logger.info(f"Attempting to delete molecule with ID {mol_id}")
    db_molecule = crud.delete_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found for deletion")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"Molecule with ID {mol_id} successfully deleted")
    return db_molecule


@app.get("/molecules", response_model=List[schemas.Molecule])
def retrieve_molecules(skip: int = 0, limit: int = Query(default=100, gt=0),
                       db: Session = Depends(get_db)):
    return crud.get_all_molecules(db=db, skip=skip, limit=limit)


@app.get("/search", response_model=List[str])
def search_molecules_by_smile(
        substructure_smile: str =
        Query(..., description="SMILES string of the substructure."),
        db: Session = Depends(get_db)):
    logger.info(f"Searching for molecules with substructure: "
                f"{substructure_smile}")
    try:
        results = crud.search_molecules(db=db, substructure_smile=substructure_smile)
        logger.info(f"Search completed, found {len(results)} matches: {results}")
        return results
    except ValueError as e:
        logger.error(f"Search failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/add")
async def upload_file(file: UploadFile = File(...),
                      db: Session = Depends(get_db)):
    content = await file.read()

    logger.info("File uploaded for processing")
    try:
        return crud.add_molecules_from_file(
            db=db, file_content=content.decode("utf-8"))
    except ValueError as e:
        logger.error(f"File processing error: {e}")
        raise HTTPException(status_code=400, detail=str(e))


@app.post("/clear_molecules", status_code=200)
def clear_molecules(db: Session = Depends(get_db)):
    logger.info("Clearing all molecules from the database")
    return crud.clear_all_molecules(db=db)
