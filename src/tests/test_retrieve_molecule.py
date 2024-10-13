import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine, text
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


# Test for retrieving molecules
@pytest.mark.parametrize("mol_id, expected_status_code, expected_response", [
    (1, 200, {"mol_id": 1, "name": "CCO"}),
    (2, 200, {"mol_id": 2, "name": "c1ccccc1"}),
    (3, 200, {"mol_id": 3, "name": "CC(=O)O"}),
    (4, 200, {"mol_id": 4, "name": "CC(=O)Oc1ccccc1C(=O)O"}),
    (5, 404, {"detail": "Molecule not found"}),
    (10, 404, {"detail": "Molecule not found"})
])
def test_retrieve_molecule(mol_id, expected_status_code, expected_response,
                           db_session):
    # Insert test data if needed
    if mol_id in [1, 2, 3, 4]:  # Only insert molecules that should be
        # retrievable
        db_session.execute(
            text("INSERT INTO molecules (mol_id, name) "
                 "VALUES (:mol_id, :name)"),
            {"mol_id": mol_id, "name": expected_response["name"]}
        )
        db_session.commit()

    response = client.get(f"/molecules/{mol_id}")

    assert response.status_code == expected_status_code

    # Check the response JSON based on status code
    if expected_status_code == 200:
        response_data = response.json()
        assert response_data["mol_id"] == expected_response["mol_id"]
        assert response_data["name"] == expected_response["name"]
    else:
        # For 404 errors, just compare the entire response JSON
        assert response.json() == expected_response
