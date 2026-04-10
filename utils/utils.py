import sys
import subprocess
import logging
import glob
import os
import time
import psutil
import time
import pickle
#import functools
from functools import wraps
from contextlib import contextmanager
from pathlib import Path



logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def execute_cmd(command) -> None:
    logger.info(f"Executing: {command}")
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {command}, error: {e}")
        raise


def find_file(pattern: str) -> str:
    """
    Args: pattern of file 
    Returns: file path
    """
    files = glob.glob(pattern)
    if not files:
        raise PathError(f"file not find: {pattern}")
    return files[0]


def mkdir(dir) -> None:
    try:
        os.makedirs(dir, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")


def run_with_single_thread(command, **kwargs):
    '''
    Run the command in single-thread mode without affecting the main process.
    Enforce single-threading to avoid conflicts with OpenBLAS.
    '''
    logger.info(f"Executing: {command} with single thread")
    env = os.environ.copy()
    env.update({
        'OPENBLAS_NUM_THREADS': '1',
        'MKL_NUM_THREADS': '1',
        'OMP_NUM_THREADS': '1',
        'NUMEXPR_NUM_THREADS': '1',
    })

    try:
        subprocess.check_call(command, env=env, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {command}, error: {e}")
        raise

    #return subprocess.run(cmd, env=env, **kwargs)


def format_time(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    if h > 0:
        return f"{h:.0f}h {m:.0f}m {s:.2f}s"
    elif m > 0:
        return f"{m:.0f}m {s:.2f}s"
    else:
        return f"{s:.2f}s"


def timer(func):
    '''
    decorator: measure execution time
    '''
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info(f"[{func.__name__}] start...")
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        logger.info(f"[{func.__name__}] done. time: {format_time(elapsed)}")
        return result
    return wrapper


# Context manager: temporarily switch the working directory
@contextmanager
def tmp_chdir(path: Path):
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def save_pickle(x, filename):
    mkdir(filename)
    with open(filename, 'wb') as file:
        pickle.dump(x, file)
    print(filename)


def load_pickle(filename, verbose=True):
    with open(filename, 'rb') as file:
        x = pickle.load(file)
    if verbose:
        print(f'Pickle loaded from {filename}')
    return x