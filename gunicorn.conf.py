import os

from dotenv import load_dotenv

load_dotenv()

bind = os.environ.get("BIND", f"{os.environ.get('HOST', '0.0.0.0')}:{os.environ.get('PORT', '5000')}")
worker_class = "gthread"
workers = int(os.environ.get("WEB_CONCURRENCY", "1"))
threads = int(os.environ.get("GUNICORN_THREADS", "4"))
timeout = int(os.environ.get("GUNICORN_TIMEOUT", "120"))
graceful_timeout = int(os.environ.get("GUNICORN_GRACEFUL_TIMEOUT", "30"))
accesslog = "-"
errorlog = "-"
capture_output = True
loglevel = os.environ.get("GUNICORN_LOG_LEVEL", "info")
