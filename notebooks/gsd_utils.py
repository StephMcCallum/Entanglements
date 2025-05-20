import numpy as np
import signac
import hoomd
import gsd.hoomd

def relabel_center_gsd(job, start=0, stop=None, stride=1):
    """
    length = length of molecules
    chains = number of chains in the system
    """
    gsd_file = job.fn("production-combined.gsd")
    mid_chunk = np.array(range(job.doc.lengths//2 - 2, job.doc.lengths//2 + 3))
    all_chain_mid_indices = np.array([mid_chunk + job.doc.lengths * i for i in range(job.doc.num_mols)])

    original_traj = gsd.hoomd.open(gsd_file)
    with gsd.hoomd.open(job.fn("production-combined-center.gsd"), "w") as new_traj:
        for frame in original_traj[start:stop:stride]:
            frame.particles.types.append("B")
            old_ids = np.copy(frame.particles.typeid)
            for chunk in all_chain_mid_indices:
                 old_ids[chunk] = 1
            frame.particles.typeid = old_ids
            new_traj.append(frame)

def bootstrap_msd(job,length, start=0, stop=None, stride=1):
    steps_per_frame = int(5e5)
    ts = job.doc.real_time_step * 1e-15
    ts_frame = steps_per_frame * ts
    rand_list = []
    msd_list = []
    diff_const = []
    
    n = 50
    for i in range(n):
        rand_list.append(random.randint(0,400))
    for ind in rand_list:
        msd = msd_from_gsd(
            gsdfile=job.fn("production-combined-center.gsd"),
            start=ind,
            stop=(ind+397),
            atom_types="B",
            msd_mode="direct")
        msd_results = np.copy(msd.msd)
        np.array(msd_list)
        msd_list.append(msd_results)
    msd_avg = np.zeros((len(msd_list[0])))
    for j in msd_list:
        msd_avg += j
    msd_avg = msd_avg/len(msd_list)
    conv_factor = job.doc.ref_length**2
    job.doc.msd_units = "nm**2"
    msd_avg *= conv_factor
    time_array = np.arange(0, len(msd_avg), 1) * ts_frame
    for j in msd_list:
        slope, intercept, r_value, p_value, std_err = linregress(time_array[-133:], j[-133:])
        diff_const.append(slope)
    slope, intercept, r_value, p_value, std_err = linregress(time_array[-133:], msd_avg[-133:])
    diff_mean = np.mean(diff_const)
    std_dev = np.std(diff_const)
    std_err = np.std(diff_const)/np.sqrt(len(diff_const))
    diff_dict= [slope,diff_mean,std_dev,std_err]
    print(diff_dict)
    np.save(file=job.fn("diff_bs_50_400.npy"), arr=diff_dict)

def combine_gsd(job):
    og_prod_gsd = job.fn(f"production.gsd")
    next_prod_gsd = job.fn(f"production2.gsd")
    with gsd.hoomd.open(job.fn("production-combined.gsd"), "w") as new_traj:
        with gsd.hoomd.open(og_prod_gsd) as prod_traj:
            for frame in prod_traj:
                new_traj.append(frame)
        with gsd.hoomd.open(next_prod_gsd) as traj:
            for frame in traj[1:]:
                new_traj.append(frame)

