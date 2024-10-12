# src/routers/search.py

from fastapi import APIRouter, Query
from celery.result import AsyncResult
from typing import Dict
import logging

from src.celery_worker import search_molecules_task, celery

router = APIRouter(
    prefix="/search",
    tags=["search"],
    responses={404: {"description": "Not found"}},
)

logger = logging.getLogger(__name__)


@router.get("/", response_model=Dict[str, str])
async def search_molecules_by_smile(substructure_smile: str = Query(..., description="SMILES string of the substructure.")):
    logger.info(f"Searching for molecules with substructure: {substructure_smile}")
    task = search_molecules_task.delay(substructure_smile)
    return {"task_id": task.id, "status": task.status}


@router.get("/{task_id}", response_model=Dict[str, str])
async def get_task_status(task_id: str):
    task = AsyncResult(task_id, app=celery)

    if task.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task.result}
    else:
        return {"task_id": task_id, "status": task.state}
