#!/usr/bin/env python3

"""Run reptar calculations"""

from __future__ import annotations

import argparse
import datetime
import os
import time
from collections.abc import Iterable

from .. import File
from ..logger import ReptarLogger, set_log_level
from ..utils import _load_config, get_obj_from_string
from . import Data, Driver
from .utils import prep_xtb_input_lines

log = ReptarLogger(__name__)


def process_worker_kwargs(
    tasks: Iterable[str],
    config: dict[str, Any],
    worker_name: str,
) -> dict:
    r"""Prepare worker keyword arguments as they are dependent on task and worker.

    Parameters
    ----------
    tasks
        Requested computations.
    config
        Configuration file.
    worker_name
        Function name to identify the worker.
    """
    worker_kwargs = config["worker"]["kwargs"]

    # Handling worker arguments
    if "opt" in tasks:
        if worker_name == "xtb_worker":
            opt_block = config["worker"]["blocks"].get("opt", None)
            constraints = config["worker"].get("constrain", None)
            input_lines = prep_xtb_input_lines(
                charge=config["system"]["charge"],
                multiplicity=config["system"]["multiplicity"],
                opt_block=opt_block,
                constraints=constraints,
                save_traj=False,
            )
            worker_kwargs["input_lines"] = input_lines

    return worker_kwargs


def run_calcs(config: dict, ray_address: str = "") -> Data:
    log.info("Opening file")
    rfile_path = os.path.abspath(config["rfile"]["path"])
    rfile = File(rfile_path, mode="a")

    log.info("Preparing data for calculations")
    # Getting rfile information
    tasks = config["tasks"]

    source_key = config["rfile"]["source"]["key"]
    source_labels = config["rfile"]["source"]["labels"]
    dest_key = config["rfile"]["destination"]["key"]
    dest_labels = config["rfile"]["destination"]["labels"]

    data = Data()
    data.rfile = rfile
    data.prepare_tasks(tasks, source_key, source_labels, dest_key, dest_labels)

    driver_kwargs = config["driver"]["kwargs"]
    if ray_address != "":
        driver_kwargs["ray_address"] = ray_address
    driver = Driver(**driver_kwargs)
    start_slice = config["driver"]["start_slice"]
    end_slice = config["driver"]["end_slice"]

    worker_obj = get_obj_from_string(config["worker"]["path"])
    worker_name = worker_obj.__name__
    log.info("Selected worker: %s", worker_name)
    worker_kwargs = process_worker_kwargs(tasks, config, worker_name)

    log.info("Requested tasks: %s", ", ".join(tasks))

    n_todo = data.get_idxs_todo(tasks, start_slice, end_slice).shape[0]
    log.info("Running %r calculations", n_todo)
    if n_todo > 0:
        t_start = time.perf_counter()
        data = driver.run(
            worker_obj, worker_kwargs, data, tasks, start_slice, end_slice
        )
        log.info("Calculations are done")
        t_end = time.perf_counter()
        duration = str(datetime.timedelta(seconds=t_end - t_start))
        log.info("Finished in %s", duration)
    else:
        log.warning("Check your bounds if you are selecting a subset of structures")

    n_remaining = data.get_idxs_todo(tasks, start_slice, end_slice).shape[0]
    log.info("There are %r incomplete calculations remaining", n_remaining)

    return data


def main():
    parser = argparse.ArgumentParser(description="Run calculations using reptar")
    parser.add_argument(
        "config_path",
        type=str,
        nargs="?",
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--ray_address",
        type=str,
        nargs="?",
        default="",
        help="Desired ray address (will override config file)",
    )
    parser.add_argument(
        "--log_level",
        type=str,
        nargs="?",
        default="info",
        help="Desired logging level",
    )

    args = parser.parse_args()
    set_log_level(args.log_level.upper())
    config = _load_config(args.config_path)
    run_calcs(config, args.ray_address)


if __name__ == "__main__":
    main()
