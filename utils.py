import numpy as np
import signac
import hoomd
import gsd.hoomd

def relabel_center_beads(job, beads_from_center=2, start=0, stop=None, stride=1):
    """
    Function takes CG signac job as input and creates a new GSD file with copied trajectory and relabels center 5 beads as new type "B".

    Parameters
    ----------

    job : signac job

    beads_from_center : int, default 2
    
    """
    gsd_file = job.fn("production.gsd")
    mid_chunk = np.array(range(job.doc.lengths//2 - beads_from_center, job.doc.lengths//2 + (beads_from_center+1)))
    all_chain_mid_indices = np.array([mid_chunk + job.doc.lengths * i for i in range(job.doc.num_mols)])

    original_traj = gsd.hoomd.open(gsd_file)
    with gsd.hoomd.open(job.fn("new_traj.gsd"), "w") as new_traj:
        for frame in original_traj[start:stop:stride]:
            frame.particles.types.append("B")
            old_ids = np.copy(frame.particles.typeid)
            for chunk in all_chain_mid_indices:
                 old_ids[chunk] = 1
            frame.particles.typeid = old_ids
            new_traj.append(frame)