import pytest
from fastapi.testclient import TestClient
from main import Molecule, app


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


# Test for retrieving molecules
@pytest.mark.parametrize("mol_id", [1, 2, 3, 4, 5])
def test_retrieve_molecule(client, molecules, mol_id):
    response = client.get(f"/molecules/{mol_id}")

    if response.status_code == 200:

        expected_molecule = molecules[mol_id]
        assert expected_molecule is not None

        # Convert expected molecule to dictionary for comparison
        expected_data = {"mol_id": expected_molecule.mol_id,
                         "name": expected_molecule.name}

        assert response.json() == expected_data

    else:
        # When status is 404, check the error message
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
    