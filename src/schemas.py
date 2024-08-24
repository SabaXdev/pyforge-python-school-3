import re
from typing import Optional
from pydantic import BaseModel, Field, field_validator, ConfigDict


class MoleculeBase(BaseModel):
    mol_id: int = Field(..., gt=0,
                        description="Unique ID for the molecules",
                        json_schema_extra=
                        {"value_error.number.not_gt"
                         : "ID must be greater than 0"})
    name: str = Field(
        ..., min_length=1, max_length=100,
        description="Structure of chemical molecules(SMILES)",
        json_schema_extra={"value_error.any_str.min_length"
                           : "SMILES must be at least 1 character long"}
    )
    description: Optional[str] = (
        Field(None, max_length=255,
              description="Description of the molecule"))

    @field_validator("name")
    def validate_smiles(cls, value):
        if not re.match(
                r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?"
                r"(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
                r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|"
                r"[().\[\]])*((?:[1-9]\d*)?)$",
                value,
        ):
            raise ValueError("Invalid SMILES structure")
        return value


class MoleculeCreate(MoleculeBase):
    pass


class Molecule(MoleculeBase):
    model_config = ConfigDict(from_attributes=True)
