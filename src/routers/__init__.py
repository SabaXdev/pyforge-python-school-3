# src/routers/__init__.py

from src.routers.molecules import router as molecules_router
from src.routers.search import router as search_router
from src.routers.upload import router as upload_router
from src.routers.utils import router as utils_router

__all__ = ["molecules_router", "search_router", "upload_router", "utils_router"]
