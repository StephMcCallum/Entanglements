#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    parameters = OrderedDict()
    # Define some system related parameters:
    parameters["lengths"] = [20,30,50,80,100,150]
    parameters["chains"] = [50]
    parameters["density"] = [1.3]
    parameters["remove_hydrogens"] = [True]
    parameters["remove_charges"] = [True]
    parameters["pppm_resolution"] = [(16, 16, 16)]
    parameters["pppm_order"] = [4]
    parameters["remove_charges"] = [True]
    parameters["auto_scale"] = [True]
    # Define some simulation related parameters:
    parameters["kT"] = [4.2]
    parameters["n_steps"] = [5e7]
    parameters["shrink_kT"] = [8.0]
    parameters["shrink_n_steps"] = [5e5]
    parameters["shrink_period"] = [10000]
    parameters["r_cut"] = [2.5]
    parameters["dt"] = [0.0003]
    parameters["tau_kT"] = [100] # Used as a multiple of dt
    parameters["gsd_write_freq"] = [1e4]
    parameters["log_write_freq"] = [1e3]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project() # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate job
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("nvt_done", False)
        job.doc.setdefault("sample_done", False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
