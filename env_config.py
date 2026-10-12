import os

from dotenv import load_dotenv

# The production server keeps its .env outside the repository. Set ENV_FILE
# to use a different file, e.g. when testing on another workstation.
DEFAULT_ENV_FILE = "/home/administrator/.env"


def load_env() -> str:
    path = os.environ.get("ENV_FILE", DEFAULT_ENV_FILE)
    load_dotenv(path)
    return path
