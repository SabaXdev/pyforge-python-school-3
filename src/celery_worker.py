import crud
import json
import logging
import models
import redis
from celery import Celery
from fastapi import HTTPException

from database import SessionLocal, engine

celery = Celery(
    'tasks',
    broker='redis://localhost:6379/0',
    backend='redis://localhost:6379/0',
    broker_connection_retry_on_startup=True,
)

celery.conf.update(
    worker_hijack_root_logger=False,
    worker_log_format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    worker_task_log_format='%(asctime)s - %(name)s - %(levelname)s - '
                           '%(message)s',
    worker_log_color=False,
    loglevel='DEBUG'
)


models.Base.metadata.create_all(bind=engine)


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@celery.task
def search_molecules_task(substructure_smile: str):
    db = SessionLocal()
    # Check Redis cache for the substructure search result
    cache_key = f"substructure_search:{substructure_smile}"
    cached_result = get_cached_result(cache_key)

    if cached_result:
        try:
            # Decode bytes to string and convert to dictionary
            cached_result_dict = json.loads(cached_result.decode('utf-8'))
            results = cached_result_dict.get("results", [])
            if isinstance(results, list):  # Ensure results is a list
                logger.info(f"Returning cached result for substructure: "
                            f"{substructure_smile}")
                return results
        except json.JSONDecodeError:
            # Handle JSON decode error if needed
            return {}

    # Perform the substructure search if not cached
    logger.info(f"Performing substructure search for: {substructure_smile}")
    try:
        results = crud.search_molecules(
            db=db, substructure_smile=substructure_smile)

        set_cache(cache_key, {"results": results}, expiration=300)
        logger.info(f"Search completed, found {len(results)} "
                    f"matches: {results}")
        return results
    except ValueError as e:
        logger.error(f"Search failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    finally:
        db.close()


logger = logging.getLogger(__name__)
# Connect to redis
redis_client = redis.Redis(host='localhost', port=6379, db=0)


def get_cached_result(key):
    try:
        return redis_client.get(key)
    except Exception as e:
        logger.error(f"Error retrieving cache for key: {key}, Exception: {e}")
        raise


def set_cache(key: str, value: dict, expiration: int = 60):
    # Implement this function or import it if already implemented elsewhere
    redis_client.setex(key, expiration, json.dumps(value))
    print(redis_client.get('key'))
