import pytest
from fastapi.testclient import TestClient
from main import Molecule, app

client = TestClient(app)


@pytest.fixture
def molecules():
    return {
        1: Molecule(mol_id=1, name="CCO"),
        2: Molecule(mol_id=2, name="c1ccccc1"),
        3: Molecule(mol_id=3, name="CC(=O)O"),
        4: Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")
    }


# Test for adding molecules
@pytest.mark.parametrize("molecule", [
    {"mol_id": 5, "name": "CH4"},
    {"mol_id": 6, "name": "CO"}
])
def test_add_molecule(molecule, molecules):
    response = client.post(f"/molecules/{molecule["mol_id"]}", json=molecule)
    assert response.status_code == 201
    assert response.json() == molecule

    # Verify the molecule was added
    get_response = client.get(f"/molecules/{molecule["mol_id"]}")
    assert get_response.status_code == 200
    assert get_response.json() == molecule
