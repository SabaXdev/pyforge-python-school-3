import logging
from typing import Optional
from pydantic import BaseModel, Field, field_validator, ConfigDict
from rdkit import Chem


class MoleculeBase(BaseModel):
    mol_id: int = Field(..., gt=0,
                        description="Unique ID for the molecules",
                        json_schema_extra={
                            "value_error.number.not_gt":
                                "ID must be greater than 0"
                        }
                        )
    name: str = Field(
        ..., min_length=1, max_length=100,
        description="Structure of chemical molecules(SMILES)",
        json_schema_extra={
            "value_error.any_str.min_length":
                "SMILES must be at least 1 character long"
        }
    )
    description: Optional[str] = (
        Field(None, max_length=255,
              description="Description of the molecule"))


class MoleculeCreate(MoleculeBase):
    # SMILES validation only during creation
    @field_validator("name")
    def validate_smiles(cls, value):
        logging.info(f"Validating SMILES: {value}")
        # Use RDKit to validate the SMILES string
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            raise ValueError("Invalid SMILES structure")
        return value


class Molecule(MoleculeBase):
    model_config = ConfigDict(from_attributes=True)
