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
    from msibi import MSIBI, State, Bond, Angle
    import hoomd
    import os
    import numpy as np

    with job:
        job.doc["done"] = False
        if os.path.exists(job.fn("states")):
            dir_path = job.fn("states")
            os.system(f"rm -r {dir_path}")
        print("Setting up MSIBI optimizer...")
        opt = MSIBI(
            nlist=job.sp.nlist,
            integrator_method=hoomd.md.methods.ConstantVolume,
            method_kwargs={},
            thermostat=hoomd.md.methods.thermostats.MTTK,
            thermostat_kwargs={"tau": job.sp.thermostat_tau},
            dt=job.sp.dt,
            gsd_period=job.sp.n_steps[0] // 500,
            nlist_exclusions=job.sp.nlist_exclusions,
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

        print("Creating Bond objects...")
        bond_project = signac.get_project(job.sp.bond_project_path)
        bond_job = bond_project.open_job(id=job.sp.bond_job_id)

        AA_bond = Bond(
            type1=job.sp.bonds["type1"],
            type2=job.sp.bonds["type2"],
            optimize=False,
            nbins=bond_job.sp.bonds_nbins
        )
        AA_bond.set_from_file(file_path=bond_job.fn(job.sp.bonds["file_path"]))
        opt.add_force(AA_bond)

        print("Creating Angle objects...")

        AAA_angle = Angle(
            type1=job.sp.angles["type1"],
            type2=job.sp.angles["type2"],
            type3=job.sp.angles["type3"],
            optimize=True,
            nbins=job.sp.angles_nbins,
        )
        AAA_angle.set_quadratic(
            x0=job.sp.angles["x0"],
            x_min=job.sp.angles["x_min"],
            x_max=job.sp.angles["x_max"],
            k2=job.sp.angles["k2"],
            k3=job.sp.angles["k3"],
            k4=job.sp.angles["k4"],
        )
        AAA_angle.smoothing_window = job.sp.smoothing_window
        opt.add_force(AAA_angle)

        print("Running Optimization...")
        for n_iterations, n_steps, alpha in zip(job.sp.n_iterations,
                                                job.sp.n_steps,
                                                job.sp.state_alphas):
            state.alpha = alpha
            opt.run_optimization(n_steps=n_steps, n_iterations=n_iterations,
                                 backup_trajectories=True)
            AAA_angle.smooth_potential()
        # save the optimized bonds to file
        AAA_angle.save_potential(job.fn(f"{AAA_angle.name}_angle.csv"))
        AAA_angle.save_potential_history(
            job.fn(f"{AAA_angle.name}_potential_history.npy"))
        AAA_angle.plot_potentials(
            file_path=job.fn(f"{AAA_angle.name}_potential.png"), xlim=(1, np.pi), ylim=(0,70))
        AAA_angle.plot_potential_history(
            file_path=job.fn(f"{AAA_angle.name}_potential_history.png"), xlim=(1, np.pi), ylim=(-5,40))

        # save plots to file
        for state in opt.states:
            AAA_angle.save_state_data(state=state, file_path=job.fn(
                f"state_{state.name}_angle_{AAA_angle.name}_data.npz"))
            AAA_angle.plot_fit_scores(state=state, file_path=job.fn(
                f"{state.name}_{AAA_angle.name}_fitscore.png"))
            AAA_angle.plot_target_distribution(state=state, file_path=job.fn(
                f"{state.name}_{AAA_angle.name}_target_dist.png"))
            AAA_angle.plot_distribution_comparison(state=state,
                                                 file_path=job.fn(
                                                     f"{state.name}_{AAA_angle.name}_dist_comparison.png"))

        print("Optimization done")
        job.doc["done"] = True


if __name__ == "__main__":
    AngleMSIBI(environment=Fry).main()
