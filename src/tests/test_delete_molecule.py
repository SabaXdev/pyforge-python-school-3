import pytest
from fastapi.testclient import TestClient
from main import app, Molecule


@pytest.fixture
def client():
    with TestClient(app) as test_client:
        yield test_client


@pytest.fixture
def initial_molecules():
    return {
        1: Molecule(mol_id=1, name="CCO"),
        2: Molecule(mol_id=2, name="c1ccccc1"),
        3: Molecule(mol_id=3, name="CC(=O)O"),
        4: Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")
    }


@pytest.fixture(autouse=True)
def reset_molecules(client):
    # Clear the in-memory storage before each test
    client.post("/clear_molecules")
    yield
    client.post("/clear_molecules")  # Clean up after each test


@pytest.mark.parametrize("mol_id", [1, 2, 3, 4, 5])
def test_delete_molecule(client, initial_molecules, mol_id):
    # Setup initial data
    for mol in initial_molecules.values():
        client.post(f"/molecules/{mol.mol_id}", json=mol.dict())

    if mol_id in initial_molecules:
        # Ensure the molecule exists before deletion
        response = client.get(f"/molecules/{mol_id}")
        assert response.status_code == 200
        assert response.json() == {"mol_id": mol_id, "name": initial_molecules[mol_id].name}

        # Perform the delete operation
        response = client.delete(f"/molecules/{mol_id}")
        assert response.status_code == 200
        assert response.json() == {"mol_id": mol_id, "name": initial_molecules[mol_id].name}

        # Verify the molecule has been deleted
        response = client.get(f"/molecules/{mol_id}")
        assert response.status_code == 404
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
    else:
        # If the molecule didn't exist before deletion, expect a 404
        response = client.delete(f"/molecules/{mol_id}")
        assert response.status_code == 404
        assert response.json() == {"detail": f"Molecule with id {mol_id} not found."}
