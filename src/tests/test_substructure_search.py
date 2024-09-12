import pytest
import logging
from fastapi.testclient import TestClient
from main import app, get_db, celery
from models import Base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from unittest.mock import patch, MagicMock
from celery.result import AsyncResult

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Set up the SQLite database URL (for test)
SQLALCHEMY_DATABASE_URL = "sqlite:///./test_molecules.db"
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

client = TestClient(app)


# Override get_db for testing
def override_get_db():
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()


app.dependency_overrides[get_db] = override_get_db


@pytest.fixture(scope="function")
def db_session():
    Base.metadata.create_all(bind=engine)
    db = TestingSessionLocal()
    try:
        yield db
    finally:
        db.close()
        Base.metadata.drop_all(bind=engine)


# Test for /search endpoint (Celery task submission)
@pytest.mark.parametrize("substructure_smile", [
    "CC", "c1ccccc1", "O", "P"
])
def test_search_molecules_by_smile(substructure_smile, db_session):
    logger.info("Starting test for searching molecules by smile")

    # Mock Celery task to avoid actual async task processing in test
    with patch('main.search_molecules_task.delay') as mock_task:
        mock_task.return_value = MagicMock(id="123", status="PENDING")

        # Perform the search for the given substructure SMILES
        response = client.get(f"/search?substructure_smile={substructure_smile}")
        assert response.status_code == 200
        response_data = response.json()

        # Ensure task_id and status are returned
        assert "task_id" in response_data
        assert response_data["task_id"] == "123"
        assert response_data["status"] == "PENDING"

    logger.info(f"Completed test for searching substructure smile {substructure_smile}.")


# Test for /search/{task_id} endpoint (Task status retrieval)
def test_get_task_status(db_session):
    logger.info("Starting test for getting task status")

    # Mock AsyncResult to simulate task behavior
    with patch('main.AsyncResult') as mock_async_result:
        # Simulate a task that is still processing
        mock_async_result.return_value = MagicMock(state="PENDING")

        # Check task status when it's pending
        response = client.get("/search/123")
        assert response.status_code == 200
        response_data = response.json()
        assert response_data["status"] == "Task is still processing"
        assert response_data["task_id"] == "123"

        # Simulate a task that has completed successfully
        mock_async_result.return_value.state = "SUCCESS"
        mock_async_result.return_value.result = ["CCO", "c1ccccc1"]

        # Check task status when it's completed
        response = client.get("/search/123")
        assert response.status_code == 200
        response_data = response.json()
        assert response_data["status"] == "Task completed"
        assert response_data["result"] == ["CCO", "c1ccccc1"]

    logger.info("Completed test for getting task status.")
