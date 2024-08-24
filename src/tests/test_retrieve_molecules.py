import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base
import schemas  # Make sure to import your schemas module

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


# Test for retrieving molecules
def test_retrieve_molecules(db_session):
    # Add initial molecules
    initial_molecules = [
        {"mol_id": 1, "name": "CCO"},
        {"mol_id": 2, "name": "c1ccccc1"},
        {"mol_id": 3, "name": "CC(=O)O"},
        {"mol_id": 4, "name": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
    for mol in initial_molecules:
        client.post(f"/molecules/{mol['mol_id']}", json=mol)

    # Retrieve all molecules
    response = client.get("/molecules")
    assert response.status_code == 200

    # Convert the response to a set of Molecule identifiers
    response_molecules = response.json()
    response_molecules_set = {(mol["mol_id"], mol["name"])
                              for mol in response_molecules}

    # Create a set of expected molecule identifiers
    expected_molecules_set = {
        (mol["mol_id"], mol["name"]) for mol in initial_molecules
    }

    # Validate that expected and actual sets match
    assert expected_molecules_set == response_molecules_set
