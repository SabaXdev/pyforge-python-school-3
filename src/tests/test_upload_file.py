import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from main import app, get_db
from models import Base
import io

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


# Test for uploading a file and adding molecules
@pytest.mark.parametrize(
    "file_content, expected_status_code, expected_response", [
        ("5 CH4\n6 CO\n", 200,  # Valid ID and SMILES format
         {'content': 'Molecule/Molecules added successfully.'}
         ),
        ("invalid_id invalid_smiles", 400,  # Invalid ID and SMILES should
         # result in a 400 Bad Request
         {'detail': "invalid literal for int() with base 10: 'invalid_id'"})
    ]
)
def test_upload_file(file_content, expected_status_code, expected_response,
                     db_session):
    # Simulate a file upload using FastAPI's UploadFile mechanism
    file = io.BytesIO(file_content.encode('utf-8'))
    files = {'file': ('molecules.txt', file, 'text/plain')}

    # Post the file to the /add endpoint
    response = client.post("/add", files=files)

    # Check the response status code
    assert response.status_code == expected_status_code
    response_data = response.json()

    # Check the response body if the status is 200
    if response.status_code == 200:
        assert response_data == expected_response
    else:
        # Handle the expected error response
        assert response_data == expected_response


# Test for uploading a file with invalid content
def test_upload_file_invalid_content(db_session):
    # Simulate a file with entirely invalid SMILES strings
    file_content = "invalid_smiles1\ninvalid_smiles2\n"
    file = io.BytesIO(file_content.encode('utf-8'))
    files = {'file': ('molecules.txt', file, 'text/plain')}

    # Post the file to the /add endpoint
    response = client.post("/add", files=files)

    # Check that the response returns a 400 status code
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid file format"}
