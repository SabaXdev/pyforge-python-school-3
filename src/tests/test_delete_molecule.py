import pytest
from fastapi.testclient import TestClient
from main import app, Molecule


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


@pytest.mark.parametrize("mol_id", [1, 2, 3, 4, 5])
def test_delete_molecule(client, molecules, mol_id):
    if mol_id in molecules:
        # Ensure the molecule exists before deletion
        response = client.get(f"/molecules/{mol_id}")
        assert response.status_code == 200
        assert response.json() == {"mol_id": mol_id, "name": molecules[mol_id].name}

        # Perform the delete operation
        response = client.delete(f"/molecules/{mol_id}")
        assert response.status_code == 200
        assert response.json() == {"mol_id": mol_id, "name": molecules[mol_id].name}

        # Verify the molecule has been deleted
        response = client.get(f"/molecules/{mol_id}")
        assert response.status_code == 404
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
    else:
        # If the molecule didn't exist before deletion, expect a 404
        response = client.delete(f"/molecules/{mol_id}")
        assert response.status_code == 404
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
