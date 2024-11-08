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
    parameters["chains"] = [
        (50, 20),
        (50, 30),
        (50, 50),
        (50, 75),
        (50, 100),
        (50, 150),
        (50, 200),
        (50, 300),
    ]
    parameters["density"] = [1.32]
    parameters["harmonic_bonds"] = [True]
    parameters["kT"] = [4.2]
    parameters["n_equil_steps"] = [5e7]
    parameters["n_prod_steps"] = [1e8]
    parameters["shrink_kT"] = [7.0]
    parameters["n_shrink_steps"] = [5e7]
    parameters["shrink_period"] = [10000]
    parameters["dt"] = [0.005]
    parameters["tau_kT"] = [100]
    parameters["gsd_write_freq"] = [1e6]
    parameters["log_write_freq"] = [1e4]
    parameters["sim_seed"] = [42]
    # Get FF from the MSIBI Project
    parameters["msibi_project"] = [
        "/home/erjank_project/pps-entanglement-length/Entanglements/angle-flow-with-pairs"
    ]
    parameters["msibi_job"] = ["34c9e9f8fa7d942743adbf6835395671"]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project()
    param_names, param_combinations = get_parameters()
    # Create workspace of jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("equilibrated", False)
        job.doc.setdefault("sampled", False)
        job.doc.setdefault("runs", 0)
        job.doc.setdefault("num_mols", job.sp.chains[0])
        job.doc.setdefault("lengths", job.sp.chains[1])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
