import argparse
import logging
import os
from functools import wraps
import time
from datetime import timedelta


def add_log(func):
    '''
    logging start and done.
    '''
    logFormatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    module = func.__module__
    name = func.__name__
    logger_name = f'{module}.{name}'
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        if args and hasattr(args[0], 'debug') and args[0].debug:
            logger.setLevel(10)  # debug

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper