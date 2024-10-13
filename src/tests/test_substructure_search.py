import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base

# Set up the SQLite database URL (using a file for tests)
SQLALCHEMY_DATABASE_URL = "sqlite:///./test_molecules.db"

# Create the engine and session for the test
engine = create_engine(SQLALCHEMY_DATABASE_URL,
                       connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False,
                                   bind=engine)

# Create a new test client for FastAPI
client = TestClient(app)


# Override the get_db dependency to use the test database session
@pytest.fixture(scope="function")
def db_session():
    Base.metadata.create_all(bind=engine)  # Create the tables
    db = TestingSessionLocal()
    try:
        yield db  # Provide the db session to the test
    finally:
        db.close()
        Base.metadata.drop_all(bind=engine)  # Drop the tables after the test


# Dependency override for get_db
def override_get_db():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db


# Test for searching molecules by substructure SMILES
@pytest.mark.parametrize("substructure_smile, expected_result", [
    ("CC", ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("O", ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),
    ("P", []),  # No molecules containing Phosphorus
])
def test_search_molecules_by_smile(substructure_smile, expected_result,
                                   db_session):
    # Add initial molecules
    initial_molecules = [
        {"mol_id": 1, "name": "CCO"},
        {"mol_id": 2, "name": "c1ccccc1"},
        {"mol_id": 3, "name": "CC(=O)O"},
        {"mol_id": 4, "name": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
    for mol in initial_molecules:
        client.post(f"/molecules/{mol['mol_id']}", json=mol)

    # Perform the search for the given substructure SMILES
    response = client.get(f"/search?substructure_smile={substructure_smile}")
    assert response.status_code == 200

    # Validate that the search results match the expected list of molecules
    response_molecules = response.json()
    assert sorted(response_molecules) == sorted(expected_result)
