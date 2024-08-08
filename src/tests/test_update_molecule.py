import pytest
from rdkit import Chem
from fastapi.testclient import TestClient
from main import substructure_search, Molecule, app


@pytest.fixture
def client():
    with TestClient(app) as test_client:
        yield test_client


@pytest.fixture
def molecules():
    return {
        1: Molecule(mol_id=1, name="CCO"),
        2: Molecule(mol_id=2, name="c1ccccc1"),
        3: Molecule(mol_id=3, name="CC(=O)O"),
        4: Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")
    }


@pytest.mark.parametrize("mol_id, updated_molecule_data", [
    (1, {"mol_id": 1, "name": "CH3OH"}),
    (2, {"mol_id": 2, "name": "C6H6"}),
    (3, {"mol_id": 3, "name": "H2O"}),
    (4, {"mol_id": 4, "name": "CO2"}),
    (5, {"mol_id": 5, "name": "CH4"})
])
# Test for adding molecules
def test_update_molecule(client, molecules, mol_id, updated_molecule_data):
    response = client.put(f"/molecules/{mol_id}", json=updated_molecule_data)

    if response.status_code == 200:

        updated_molecule = response.json()
        expected_molecule = molecules[mol_id]
        expected_molecule.name = updated_molecule_data["name"]

        # Convert expected molecule to dictionary for comparison
        expected_data = {"mol_id": expected_molecule.mol_id,
                         "name": expected_molecule.name}

        assert updated_molecule == expected_data

    else:
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
