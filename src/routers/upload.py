# src/routers/upload.py

from fastapi import APIRouter, UploadFile, File, Depends, HTTPException
from sqlalchemy.orm import Session
import logging

from src import crud
from src.database import SessionLocal
from typing import Dict

router = APIRouter(
    prefix="/upload",
    tags=["upload"],
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


@router.post("/", response_model=Dict[str, str], status_code=200)
async def upload_file(file: UploadFile = File(...), db: Session = Depends(get_db)):
    content = await file.read()
    logger.info("File uploaded for processing")
    try:
        result = crud.add_molecules_from_file(db=db, file_content=content.decode("utf-8"))
        return {"message": "File processed successfully", "details": result}
    except ValueError as e:
        logger.error(f"File processing error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
