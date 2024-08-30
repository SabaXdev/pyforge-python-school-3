import logging

import pytest
from fastapi.testclient import TestClient
from main import app
import os

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


@pytest.mark.parametrize("server_id, expected", [
    ("1", "1"),
    ("2", "2"),
    ("my-server", "my-server"),
])
def test_get_server_with_param(server_id, expected):
    logger.info(f"Setting SERVER_ID environment variable to: {server_id}")
    os.environ["SERVER_ID"] = server_id

    client = TestClient(app)

    logger.info(f"Sending GET request to / with SERVER_ID: {server_id}")
    response = client.get("/")

    logger.info(f"Received response: {response.json()} "
                f"with status code: {response.status_code}")
    assert response.status_code == 200
    assert response.json() == {"server_id": expected}
    logger.info("Test passed successfully")
