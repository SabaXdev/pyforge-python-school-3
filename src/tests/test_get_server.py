import pytest
from fastapi.testclient import TestClient
from main import app
import os


@pytest.mark.parametrize("server_id, expected", [
    ("1", "1"),
    ("2", "2"),
    ("my-server", "my-server"),
])
def test_get_server_with_param(server_id, expected):
    os.environ["SERVER_ID"] = server_id
    client = TestClient(app)
    response = client.get("/")
    assert response.status_code == 200
    assert response.json() == {"server_id": expected}
