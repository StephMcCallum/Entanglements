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

import logging
from collections import OrderedDict
from itertools import product

import numpy as np
import signac


def get_parameters(ordered_dict=OrderedDict()):
    '''Use the listed parameters below to set up
    your MSIBI instructions.

    '''
    parameters = ordered_dict

    # Optimizer parameters
    parameters["nlist"] = ["Cell"]
    parameters["thermostat_tau"] = [0.03]
    parameters["dt"] = [0.0003]
    parameters["n_steps"] = [[1e6]]
    parameters["state_alphas"] = [[0.6]]
    parameters["n_iterations"] = [[10]]

    # State parameters
    parameters["single_chain_path"] = [
        "/home/erjank_project/PPS-MSIBI/pps-msibi/training-runs/single-chains"
    ]
    parameters["single_chain_job_id"] = ["29a7f0d216700e7c8534b8c11140ba06"]
    parameters["states"] = [
        [
            {"name": "A",
             "remove_hydrogens": True,
             "alpha": 0.6,
             "n_frames": 100,
             "cg_file_name": "target_1monomer_per_bead.gsd"
             },
        ],

    ]

    # Pair parameters
    parameters["pair_project_path"] = [
        "/home/erjank_project/PPS-MSIBI/pps-msibi/msibi-flow/pair-flow"
    ]
    parameters["pair_job_id"] = ["b4289482ca51d698463f1c62717db7d8"]
    parameters["pairs"] = [
            {
                "type1": "A",
                 "type2": "A",
                 "file_path": "pair_pot.csv",
             },
    ]

    # Angle parameters
    # Get final potential from initial angle run
    parameters["angle_project_path"] = [
        "/home/erjank_project/PPS-MSIBI/pps-msibi/msibi-flow/angle-flow"
    ]
    parameters["angle_job_id"] = ["d23d8c3ba38040009adda789dca01050"]
    parameters["head_correction"] = ["linear"]
    parameters["angles_nbins"] = [100]
    parameters["angles"] = [
            {
                "type1": "A",
                 "type2": "A",
                 "type3": "A",
                 "file_path": "A-A-A_angle.csv",
             },
    ]
    parameters["smoothing_window"] = [9]

    # Bond parameters
    parameters["bond_project_path"] = [
            "/home/erjank_project/PPS-MSIBI/pps-msibi/msibi-flow/bond-flow-with-pairs"
    ]
    parameters["bond_job_id"] = ["74161d051db9fe68ea155f379284e8d9"]
    parameters["bonds_nbins"] = [60]
    parameters["bonds"] = [
            {
                "type1": "A",
                 "type2": "A",
                 "file_path": "A-A_bond.csv",
             },
    ]

    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project()  # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("done", False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
