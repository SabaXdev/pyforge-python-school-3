from fastapi.testclient import TestClient
from main import app, Molecule
import io
import pytest


@pytest.fixture
def client():
    with TestClient(app) as test_client:
        yield test_client




@pytest.fixture
def sample_file():
    # Create a sample file content with molecule data
    content = "5 CH4\n6 CO"  # Ensure this matches the expected format
    return io.BytesIO(content.encode("utf-8"))


@pytest.fixture
def initial_molecules():
    return [
        Molecule(mol_id=1, name="CCO"),
        Molecule(mol_id=2, name="c1ccccc1"),
        Molecule(mol_id=3, name="CC(=O)O"),
        Molecule(mol_id=4, name="CC(=O)Oc1ccccc1C(=O)O")
    ]



def setup_initial_data(client, molecules):
    for mol in molecules.values():
        client.post("/molecules", json=mol.dict())


def test_upload_file(client, initial_molecules, sample_file):
    # Send POST request with sample file
    response = client.post("/add", files={"file": ("sample.txt", sample_file, "text/plain")})
    assert response.status_code == 200
    assert response.json() == {"content": "Molecule/Molecules added successfully."}

    # Retrieve all molecules to verify the new ones are added
    response = client.get("/molecules")
    assert response.status_code == 200

    # Convert the response to a set of Molecule identifiers
    response_molecules = response.json()
    response_molecules_set = {(mol["mol_id"], mol["name"]) for mol in response_molecules}

    # Create a set of expected molecule identifiers
    expected_molecules_set = {
        (mol.mol_id, mol.name) for mol in initial_molecules
    }

    # Add newly added molecules to the set
    expected_molecules_set.update({
        (5, "CH4"),
        (6, "CO")
    })

    # Validate that expected and actual sets match
    assert expected_molecules_set == response_molecules_set



