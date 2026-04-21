import logging
import sys
import time
from functools import wraps
from datetime import timedelta

import loompy as lp
import numpy as np
import scanpy as sc
import argparse
import subprocess


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

    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):

        logger.info('start...')
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        return result

    wrapper.logger = logger
    return wrapper


def h5ad2loom(h5ad, outdir):
    adata = sc.read_h5ad(f'{h5ad}')
    row_attrs = { 
            "Gene": np.array(adata.var_names) 
        }
    col_attrs = {
            "CellID": np.array(adata.obs_names) ,
            "nGene": np.array( np.sum(adata.X.transpose()>0 , axis = 0)).flatten() ,
            "nUMI": np.array( np.sum(adata.X.transpose() , axis = 0)).flatten() ,
        }
    lp.create(f'{outdir}/pyscenic_input.loom', adata.X.transpose(), row_attrs, col_attrs)

def mtx2loom(mtx, outdir):
    '''
    Not recommended because visualization requires h5ad files.
    '''
    mtx = sc.read_csv(f'{mtx}')
    row_attrs = {"Gene": np.array(mtx.var_names),}
    col_attrs = {"CellID": np.array(mtx.obs_names)}
    lp.create(f'{outdir}/pyscenic_input.loom', mtx.X.transpose(), row_attrs, col_attrs)