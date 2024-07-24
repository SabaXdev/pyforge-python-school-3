from pydantic import BaseModel


class Molecule(BaseModel):
    mol_id: int
    name: str
