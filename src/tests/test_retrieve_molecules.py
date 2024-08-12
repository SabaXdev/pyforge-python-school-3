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


@pytest.fixture(autouse=True)
def reset_molecules(client):
    # Clear the in-memory storage before each test
    client.post("/clear_molecules")
    yield
    client.post("/clear_molecules")  # Clean up after each test


# Test for retrieving molecules
def test_retrieve_molecules(client, molecules):
    # Setup initial data
    for mol in molecules.values():
        client.post(f"/molecules/{mol.mol_id}",
                    json=mol.model_dump())  # Use model_dump()

    # Retrieve all molecules
    response = client.get("/molecules")
    assert response.status_code == 200

    # Convert the response to a set of Molecule identifiers
    response_molecules = response.json()
    response_molecules_set = {(mol["mol_id"], mol["name"])
                              for mol in response_molecules}

    # Create a set of expected molecule identifiers
    expected_molecules_set = {
        (mol.mol_id, mol.name) for mol in molecules.values()
    }

    # Validate that expected and actual sets match
    assert expected_molecules_set == response_molecules_set
