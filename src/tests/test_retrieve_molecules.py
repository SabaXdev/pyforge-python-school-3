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


# Test for retrieving molecules
def test_retrieve_molecules(client, molecules):
    response = client.get(f"/molecules")
    assert response.status_code == 200

    response_molecules = response.json()
    assert len(response_molecules) == len(molecules)

    for i, molecule in enumerate(response.json()):
        assert molecule["mol_id"] == molecules[i+1].mol_id
        assert molecule["name"] == molecules[i+1].name
