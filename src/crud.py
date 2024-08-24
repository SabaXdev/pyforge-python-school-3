from rdkit import Chem
from sqlalchemy.orm import Session

import models, schemas


def create_molecule(db: Session, molecule: schemas.MoleculeCreate):
    db_molecule = models.Molecule(mol_id=molecule.mol_id, name=molecule.name)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    return db_molecule


def get_molecule_by_id(db: Session, molecule_id: int):
    return (db.query(models.Molecule).
            filter(models.Molecule.mol_id == molecule_id).first())


def update_molecule_by_id(db: Session, molecule_id: int,
                          updated_molecule: schemas.MoleculeCreate):
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id ==
                                                   molecule_id).first()
    if db_molecule:
        db_molecule.name = updated_molecule.name
        db.commit()
        db.refresh(db_molecule)
        return db_molecule
    return None


def delete_molecule_by_id(db: Session, molecule_id: int):
    db_molecule = db.query(models.Molecule).filter(models.Molecule.mol_id ==
                                                   molecule_id).first()
    if db_molecule:
        db.delete(db_molecule)
        db.commit()
        return db_molecule
    return None


def get_all_molecules(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Molecule).offset(skip).limit(limit).all()


def search_molecules(db: Session, substructure_smile: str):
    smiles = db.query(models.Molecule.name).all()
    smiles = [smile[0] for smile in smiles]

    substructure = Chem.MolFromSmiles(substructure_smile)
    if not substructure:
        raise ValueError("Invalid substructure SMILES.")

    return [smile for smile in smiles if Chem.MolFromSmiles(smile).
    HasSubstructMatch(substructure)]


def add_molecules_from_file(db: Session, file_content: str):
    lines = file_content.splitlines()
    existing_ids = {mol.mol_id for mol in db.query(models.Molecule.mol_id).
    all()}

    for line in lines:
        parts = line.split()
        if len(parts) != 2:
            raise ValueError("Invalid file format")

        mol_id, name = parts
        mol_id = int(mol_id)

        if mol_id not in existing_ids:
            db_molecule = models.Molecule(mol_id=mol_id, name=name)
            db.add(db_molecule)
            existing_ids.add(mol_id)

    db.commit()
    return {"content": "Molecule/Molecules added successfully."}


def clear_all_molecules(db: Session):
    db.query(models.Molecule).delete()
    db.commit()
    return {"detail": "Molecules cleared successfully"}