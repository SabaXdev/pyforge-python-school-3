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


@pytest.fixture
def sample_file():
    # Create a file-like object with sample data
    from io import BytesIO
    file_content = "5 CH4\n6 CO\n".encode('utf-8')
    return BytesIO(file_content)


def test_upload_file(client, initial_molecules, sample_file):
    # Send POST request with sample file
    response = client.post("/add", files={"file": ("sample.txt", sample_file, "text/plain")})
    assert response.status_code == 200
    assert response.json() == {"content": "Molecule/Molecules added successfully."}

    # Retrieve all molecules to verify the new ones are added
    response = client.get("/molecules")
    assert response.status_code == 200

    # Convert the response to a dictionary of Molecules
    response_molecules = {mol["mol_id"]: Molecule(**mol) for mol in response.json()}

    # Create a set of expected molecule identifiers
    expected_molecules_set = {
        (mol.mol_id, mol.name) for mol in initial_molecules.values()
    }

    # Add newly added molecules to the set
    expected_molecules_set.update({
        (5, "CH4"),
        (6, "CO")
    })

    # Create a set of actual response molecule identifiers
    response_molecules_set = {
        (mol["mol_id"], mol["name"]) for mol in response.json()
    }

    # Validate that expected and actual sets match
    assert expected_molecules_set == response_molecules_set
