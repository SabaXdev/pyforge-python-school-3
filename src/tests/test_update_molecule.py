import logging

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")


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
    logger.info("Setting up the test database.")
    Base.metadata.create_all(bind=engine)  # Create the tables
    db = TestingSessionLocal()
    try:
        yield db  # Provide the db session to the test
    finally:
        db.close()
        Base.metadata.drop_all(bind=engine)  # Drop the tables after the test
        logger.info("Test database teardown completed.")


# Dependency override for get_db
def override_get_db():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db


# Test for updating molecules
@pytest.mark.parametrize(
    "mol_id, update_data, expected_status_code, expected_response", [
        (1, {"mol_id": 1, "name": "C"}, 200,
         {"mol_id": 1, "name": "C"}),
        (2, {"mol_id": 2, "name": "P"}, 200,
         {"mol_id": 2, "name": "P"}),
        (3, {"mol_id": 3, "name": "CCO"}, 200,
         {"mol_id": 3, "name": "CCO"}),
        (4, {"mol_id": 4, "name": "CO"}, 200,
         {"mol_id": 4, "name": "CO"}),
        (5, {"mol_id": 5, "name": "CN1CCC[C@H]1c2cccnc2"}, 404,
         {"detail": "Molecule not found"})
    ]
)
def test_update_molecule(db_session, mol_id, update_data,
                         expected_status_code, expected_response):
    logger.info("Starting test for updating molecules")
    # Add initial molecules
    initial_molecules = [
        {"mol_id": 1, "name": "CCO"},
        {"mol_id": 2, "name": "c1ccccc1"},
        {"mol_id": 3, "name": "CC(=O)O"},
        {"mol_id": 4, "name": "CC(=O)Oc1ccccc1C(=O)O"}
    ]
    for mol in initial_molecules:
        response = client.post(f"/molecules/{mol['mol_id']}", json=mol)
        logger.info(f"API POST Response: {response.json()}")

    response = client.put(f"/molecules/{mol_id}", json=update_data)
    logger.info(f"API PUT Response: {response.json()}")

    # Check the status code first
    assert response.status_code == expected_status_code

    # Check the response JSON based on status code
    if expected_status_code == 200:
        response_data = response.json()
        assert response_data["mol_id"] == expected_response["mol_id"]
        assert response_data["name"] == expected_response["name"]
    else:
        # For 404 errors, just compare the entire response JSON
        assert response.json() == expected_response
