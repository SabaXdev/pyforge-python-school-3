# src/routers/molecules.py

from typing import List
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
import logging

from src import schemas, crud
from src.database import SessionLocal

router = APIRouter(
    prefix="/molecules",
    tags=["molecules"],
    responses={404: {"description": "Not found"}},
)

logger = logging.getLogger(__name__)


# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@router.post("/", response_model=schemas.Molecule, status_code=201)
def add_molecule(molecule: schemas.MoleculeCreate, db: Session = Depends(get_db)):
    return crud.create_molecule(db=db, molecule=molecule)


@router.get("/{mol_id}", response_model=schemas.Molecule)
def retrieve_molecule(mol_id: int, db: Session = Depends(get_db)):
    logger.info(f"Retrieving molecule with ID {mol_id}")
    db_molecule = crud.get_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    return db_molecule


@router.put("/{mol_id}", response_model=schemas.Molecule)
def update_molecule(mol_id: int, molecule: schemas.MoleculeCreate, db: Session = Depends(get_db)):
    logger.info(f"Updating molecule with ID {mol_id}")
    db_molecule = crud.update_molecule_by_id(db=db, molecule_id=mol_id, updated_molecule=molecule)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found for update")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"Molecule with ID {mol_id} successfully updated")
    return db_molecule


@router.delete("/{mol_id}", response_model=schemas.Molecule, status_code=200)
def delete_molecule(mol_id: int, db: Session = Depends(get_db)):
    logger.info(f"Attempting to delete molecule with ID {mol_id}")
    db_molecule = crud.delete_molecule_by_id(db=db, molecule_id=mol_id)
    if db_molecule is None:
        logger.warning(f"Molecule with ID {mol_id} not found for deletion")
        raise HTTPException(status_code=404, detail="Molecule not found")
    logger.info(f"Molecule with ID {mol_id} successfully deleted")
    return db_molecule


@router.get("/", response_model=List[schemas.Molecule])
def retrieve_molecules(skip: int = 0, limit: int = Query(default=100, gt=0), db: Session = Depends(get_db)):
    return crud.get_all_molecules(db=db, skip=skip, limit=limit)
