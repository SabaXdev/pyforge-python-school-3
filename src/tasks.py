from celery_worker import celery


@celery.task
def add_task(x: int, y: int):
    return x + y
