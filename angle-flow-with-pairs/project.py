"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help

"""

import os

import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment


class AngleMSIBI(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition", default="gpu",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition", default="batch",
            help="Specify the partition to submit to."
        )


# Definition of project-related labels (classification)
@AngleMSIBI.label
def completed(job):
    return job.doc.get("done")


def get_file(job, file_name):
    return os.path.abspath(os.path.join(job.ws, "..", "..", file_name))



@AngleMSIBI.post(completed)
@AngleMSIBI.operation(
    directives={"ngpu": 1, "executable": "python -u"}, name="optimize"
)
def optimize(job):
    from msibi import MSIBI, State, Bond, Angle, Pair
    import hoomd
    import os
    import numpy as np

    with job:
        job.doc["done"] = False
        if os.path.exists(job.fn("states")):
            dir_path = job.fn("states")
            os.system(f"rm -r {dir_path}")

        # Open up the pair project being used:
        pair_project = signac.get_project(job.sp.pair_project_path)
        pair_job = pair_project.open_job(id=job.sp.pair_job_id)

        print("Setting up MSIBI optimizer...")
        opt = MSIBI(
            nlist=hoomd.md.nlist.Cell,
            integrator_method=hoomd.md.methods.ConstantVolume,
            method_kwargs={},
            thermostat=hoomd.md.methods.thermostats.MTTK,
            thermostat_kwargs={"tau": job.sp.thermostat_tau},
            dt=job.sp.dt,
            gsd_period=job.sp.n_steps[0] // 500,
            nlist_exclusions=pair_job.sp.nlist_exclusions,
        )

        print("Creating State objects...")
        single_chain_project = signac.get_project(job.sp.single_chain_path)
        single_chain_job = single_chain_project.open_job(
            id=job.sp.single_chain_job_id
        )
        job.doc.target_state_path = single_chain_job.path
        for idx, state in enumerate(job.sp.states):
            print("state: ", state)
            gsd_file = single_chain_job.fn(state["cg_file_name"])
            state = State(
                name=state["name"],
                kT=single_chain_job.sp.kT,
                traj_file=gsd_file,
                n_frames=state["n_frames"],
                alpha=job.sp.state_alphas[0],
            )
            opt.add_state(state)

        print("Creating bond objects...")
        bond_project = signac.get_project(job.sp.bond_project_path)
        bond_job = bond_project.open_job(id=job.sp.bond_job_id)
        AA_bond = Bond(
            type1=pair_job.sp.bonds["type1"],
            type2=pair_job.sp.bonds["type2"],
            optimize=False,
            nbins=bond_job.doc.bonds_nbins,
        )
        AA_bond.set_from_file(file_path=bond_job.fn(job.sp.bonds["file_path"]))
        opt.add_force(AA_bond)

        print("Creating angle objects...")
        angle_project = signac.get_project(pair_job.sp.angle_project_path)
        angle_job = angle_project.open_job(id=pair_job.sp.angle_job_id)
        job.doc.angles_nbins = angle_job.sp.angles_nbins
        AAA_angle = Angle(
            type1=pair_job.sp.angles["type1"],
            type2=pair_job.sp.angles["type2"],
            type3=pair_job.sp.angles["type3"],
            optimize=True,
            nbins=angle_job.sp.angles_nbins,
        )
        AAA_angle.set_from_file(file_path=angle_job.fn(pair_job.sp.angles["file_path"]))
        opt.add_force(AAA_angle)

        print("Creating pair objects...")
        AA_pair = Pair(
            type1=pair_job.sp.pairs["type1"],
            type2=pair_job.sp.pairs["type2"],
            optimize=False,
            nbins=pair_job.sp.pairs_nbins,
            r_cut=pair_job.sp.r_cut
        )
        AA_pair.set_from_file(file_path=pair_job.fn(job.sp.pairs["file_path"]))
        opt.add_force(AA_pair)

        print("Running Optimization...")
        for n_iterations, n_steps, alpha in zip(
                job.sp.n_iterations,
                job.sp.n_steps,
                job.sp.state_alphas
        ):
            state.alpha = alpha
            opt.run_optimization(
                    n_steps=n_steps,
                    n_iterations=n_iterations,
                    backup_trajectories=True
            )
            AAA_angle.smooth_potential()
        print("Optimization done")

        # save the optimized pairs to file
        AAA_angle.save_potential(job.fn(f"{AAA_angle.name}_angle.csv"))
        AAA_angle.save_potential_history(
            job.fn(f"{AAA_angle.name}_potential_history.npy")
        )
        # save plots to file
        for state in opt.states:
            AAA_angle.save_state_data(
                    state=state,
                    file_path=job.fn(
                        f"state_{state.name}_angle_{AAA_angle.name}_data.npz"
                    )
            )
        opt.pickle_forces(job.fn("pps-msibi.pickle"))
        job.doc["done"] = True


if __name__ == "__main__":
    AngleMSIBI(environment=Fry).main()
