"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os


class MyProject(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpu",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )

# Definition of project-related labels (classification)
@MyProject.label
def nvt_done(job):
    return job.doc.nvt_done


@MyProject.label
def sample_done(job):
    return job.doc.sample_done


@MyProject.post(nvt_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="nvt"
)
def run_nvt(job):
    import unyt as u
    from unyt import Unit
    import flowermd
    from flowermd.base.system import Pack
    from flowermd.base.simulation import Simulation
    from flowermd.library import PPS, BeadSpring
    from flowermd.utils import get_target_box_mass_density

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
    
        molecules = PPS(lengths=job.sp.lengths, num_mols=job.sp.chains)
        molecules.coarse_grain(beads={"A": "c1ccc(S)cc1"})
        ff = BeadSpring(
                r_cut=2.5,
                beads={
                    "A": dict(epsilon=1, sigma=0.2),
                    },
                bonds={
                    "A-A": dict(r0=0.64, k=500),
                    },
                angles={"A-A-A": dict(t0=2.8, k=50)},)
        system = Pack(molecules=molecules, packing_expand_factor=9, density=job.sp.density) 
        

        gsd_path = job.fn("trajectory.gsd")
        log_path = job.fn("log.txt")

        sim = Simulation.from_system(
                system=system,
                dt=job.sp.dt,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        sim.save_restart_gsd(job.fn("init.gsd"))
        # Store unit information in job doc
        tau_kT = sim.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        job.doc.ref_mass = sim.reference_mass.to("amu").value
        job.doc.ref_mass_units = "amu"
        job.doc.ref_energy = sim.reference_energy.to("kJ/mol").value
        job.doc.ref_energy_units = "kJ/mol"
        job.doc.ref_length = sim.reference_length.to("nm").value
        job.doc.ref_length_units = "nm"
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"
        # Set up stuff for shrinking volume step
        print("Running shrink step.")
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        target_box = get_target_box_mass_density(
                mass=system.mass.to("g"),
                density=job.sp.density * (Unit("g/cm**3"))
        )
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.shrink_n_steps,
                period=job.sp.shrink_period,
                tau_kt=tau_kT,
                kT=shrink_kT_ramp
        )
        print("Shrink step finished.")
        print("Running simulation.")
        sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.n_steps, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("restart.gsd"))
        job.doc.nvt_done = True
        print("Simulation finished.")


@MyProject.pre(nvt_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 0, "executable": "python -u"}, name="sample"
)
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here
        job.doc.sample_done = True


if __name__ == "__main__":
    MyProject(environment=Fry).main()
