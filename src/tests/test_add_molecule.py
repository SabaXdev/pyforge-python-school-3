import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base

# Set up the SQLite database URL (using a file for tests)
SQLALCHEMY_DATABASE_URL = "sqlite:///./test_molecules.db"

# Create the engine and session for the test
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

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


# Test for adding molecules
@pytest.mark.parametrize("molecule_data", [
    {"mol_id": 5, "name": "CH4"},
    {"mol_id": 6, "name": "CO"}
])
def test_add_molecule(molecule_data, db_session):
    response = client.post(f"/molecules/{molecule_data['mol_id']}", json=molecule_data)
    assert response.status_code == 201

    # Extract only the relevant fields from the response and compare them
    response_data = response.json()
    assert response_data["mol_id"] == molecule_data["mol_id"]
    assert response_data["name"] == molecule_data["name"]

    # Verify the molecule was added by fetching it from the API
    get_response = client.get(f"/molecules/{molecule_data['mol_id']}")
    assert get_response.status_code == 200

    get_response_data = get_response.json()
    assert get_response_data["mol_id"] == molecule_data["mol_id"]
    assert get_response_data["name"] == molecule_data["name"]
