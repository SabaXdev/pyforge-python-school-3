from rdkit import Chem
from sqlalchemy.orm import Session
import models
import schemas
import logging

logger = logging.getLogger(__name__)


def create_molecule(db: Session, molecule: schemas.MoleculeCreate):
    logger.info("Creating molecule with ID %d and name %s", molecule.mol_id, molecule.name)
    db_molecule = models.Molecule(mol_id=molecule.mol_id, name=molecule.name, description=molecule.description)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    logger.info("Molecule created successfully: %s", db_molecule.name)
    return db_molecule


def get_molecule_by_id(db: Session, molecule_id: int):
    logger.info("Fetching molecule with ID %d", molecule_id)
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == molecule_id).first()
    if db_molecule:
        logger.info("Molecule found: %s", db_molecule.name)
    return db_molecule


def update_molecule_by_id(db: Session, molecule_id: int, updated_molecule: schemas.MoleculeCreate):
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == molecule_id).first()
    if db_molecule:
        db_molecule.name = updated_molecule.name
        db_molecule.description = updated_molecule.description
        db.commit()
        db.refresh(db_molecule)
        return db_molecule
    return None


def delete_molecule_by_id(db: Session, molecule_id: int):
    logger.info("Deleting molecule with ID %d", molecule_id)
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id == molecule_id).first()
    if db_molecule:
        db.delete(db_molecule)
        db.commit()
        return db_molecule
    return None


def get_all_molecules(db: Session, skip: int = 0, limit: int = 100):
    logger.info("Fetching all molecules with skip=%d and limit=%d", skip, limit)

    # Query the molecules, applying skip and limit.
    molecules = db.query(models.Molecule).offset(skip).limit(limit)

    # Yield turns the function into a generator(iterator)
    count = 0
    for molecule in molecules:
        logger.info("Yielding molecule with mol_id=%d", molecule.mol_id)
        yield molecule
        count += 1

    logger.info("Yielded a total of %d molecules", count)
    return molecules


def search_molecules(db: Session, substructure_smile: str):
    smiles = db.query(models.Molecule.name).all()
    smiles = [smile[0] for smile in smiles]

    substructure = Chem.MolFromSmiles(substructure_smile)
    if not substructure:
        logger.error("Invalid substructure SMILES: %s", substructure_smile)
        raise ValueError("Invalid substructure SMILES.")

    matched_smiles = [smile for smile in smiles if Chem.MolFromSmiles(smile).HasSubstructMatch(substructure)]
    return matched_smiles


def add_molecules_from_file(db: Session, file_content: str):
    lines = file_content.splitlines()
    existing_ids = {mol.mol_id for mol in db.query(models.Molecule.mol_id).all()}

    for line in lines:
        parts = line.split(maxsplit=2)  # maxsplit=2 allows for up to 3 parts: id, name, and optional description
        if len(parts) < 2:
            logger.error("Invalid file format: %s", line)
            raise ValueError("Invalid file format. Must contain at least mol_id and name.")

        mol_id, name = int(parts[0]), parts[1]
        description = parts[2] if len(parts) == 3 else None  # Optional description

        if mol_id not in existing_ids:
            db_molecule = models.Molecule(mol_id=mol_id, name=name, description=description)
            logger.info('Adding molecule %s with ID %d', name, mol_id)
            db.add(db_molecule)
            existing_ids.add(mol_id)

    db.commit()
    logger.info("Molecule/Molecules added successfully.")
    return {'content': "Molecule/Molecules added successfully."}


def clear_all_molecules(db: Session):
    logger.info("Clearing all molecules from the database")
    db.query(models.Molecule).delete()
    db.commit()
    logger.info("All molecules cleared successfully.")
    return {'content': "Molecule/Molecules added successfully."}
