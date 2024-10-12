# src/routers/utils.py

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
import logging
import redis

from src import crud
from src.database import SessionLocal
from typing import Dict
from os import getenv

router = APIRouter(
    prefix="/utils",
    tags=["utils"],
    responses={404: {"description": "Not found"}},
)

logger = logging.getLogger(__name__)


# Dependency
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# Initialize Redis client
redis_client = redis.Redis(host='localhost', port=6379, db=0)


@router.get("/server")
def get_server():
    logger.info("Fetching server ID")
    return {"server_id": getenv("SERVER_ID", "1")}


@router.post("/clear_molecules", status_code=200)
def clear_molecules(db: Session = Depends(get_db)):
    logger.info("Clearing all molecules from the database")
    crud.clear_all_molecules(db=db)
    redis_client.flushdb()
    logger.info("Cache cleared from Redis")
    return {"message": "All molecules and cache cleared"}
