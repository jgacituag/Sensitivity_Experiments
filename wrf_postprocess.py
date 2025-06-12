#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from wrf import getvar, ALL_TIMES, udhel, to_np
from metpy.calc import most_unstable_cape_cin, surface_based_cape_cin
from metpy.units import units
import csv
import warnings

warnings.filterwarnings("ignore")

# -------------------- Configuration --------------------
data_root = "/home/jorge.gacitua/datosmunin2/EXPERIMENTS_UNWEATHER/DATA/TEST_Multi_8_vars_996"
csv_path  = "/home/jorge.gacitua/datosmunin/Sensitivity_Experiments/sampling_batch_996.csv"

# -------------------- Utility Functions --------------------
def extract_profile(ncfile):
    print("Extracting profile data...")
    t = getvar(ncfile, "temp", units='degC', timeidx=0)[:, 0, 0]
    td = getvar(ncfile, "td", units='degC', timeidx=0)[:, 0, 0]
    z = getvar(ncfile, "height", timeidx=0)[:, 0, 0]
    p = getvar(ncfile, "pressure", timeidx=0)[:, 0, 0]
    theta = getvar(ncfile, "theta", units='K', timeidx=0)[:, 0, 0]
    thetae = getvar(ncfile, "theta_e", units='K', timeidx=0)[:, 0, 0]
    u = getvar(ncfile, "ua", units='ms-1', timeidx=0)[:, 0, 0]
    v = getvar(ncfile, "va", units='ms-1', timeidx=0)[:, 0, 0]
    qv = getvar(ncfile, "QVAPOR", timeidx=0)[:, 0, 0]
    return dict(p=p, z=z, t=t, td=td, theta=theta, thetae=thetae, u=u, v=v, qv=qv)

def extract_reflectivity_cube(ncfile):
    print("Extracting reflectivity...")
    dbz = getvar(ncfile, "dbz", timeidx=ALL_TIMES, method="cat")  # (time, z, y, x)
    dbz_maxz = np.nanmax(dbz, axis=1)  # max over z -> (time, y, x)
    return np.transpose(dbz_maxz, axes=(2, 1, 0))  # -> (x, y, time)

def extract_updraft_helicity(ncfile):
    print("Extracting updraft helicity...")
    u = getvar(ncfile, 'ua', timeidx=ALL_TIMES)
    v = getvar(ncfile, 'va', timeidx=ALL_TIMES)
    w = getvar(ncfile, 'W', timeidx=ALL_TIMES)
    z = getvar(ncfile, 'zstag', timeidx=ALL_TIMES)
    mapfct = np.ones_like(z[:,0,:,:])
    
    uh_1_1p5 = to_np(udhel(zstag=z, mapfct=mapfct, u=u, v=v, wstag=w,bottom=1000, top=1500, dx=2000, dy=2000)).max(axis=(1, 2))
    uh_1_6   = to_np(udhel(zstag=z, mapfct=mapfct, u=u, v=v, wstag=w,bottom=1000, top=6000, dx=2000, dy=2000)).max(axis=(1, 2))
    return uh_1_1p5, uh_1_6

def extract_cape_mucape(ncfile):
    print("Extracting CAPE and MUCAPE...")
    p = getvar(ncfile, "pressure", timeidx=0)[:, 0, 0] * units.hPa
    t = getvar(ncfile, "temp", units='K', timeidx=0)[:, 0, 0] * units.K
    td = getvar(ncfile, "td", units='K', timeidx=0)[:, 0, 0] * units.K
    # check if p decreases with height
    if not np.all(np.diff(p) < 0):
        prof_p = p[::-1]
        prof_t = t[::-1]
        prof_td = td[::-1]
    else:
        prof_p = p
        prof_t = t
        prof_td = td

    mucape,mucin = most_unstable_cape_cin(prof_p, prof_t, prof_td)
    sbcape,sbcin = surface_based_cape_cin(prof_p, prof_t, prof_td)
    return float(sbcape.magnitude), float(mucape.magnitude)

def extract_scalar_timeseries(ncfile):
    print("Extracting scalar time series...")
    times = to_np(getvar(ncfile, "XTIME"))
    w = to_np(getvar(ncfile, "wa", timeidx=ALL_TIMES, method="cat"))
    t = to_np(getvar(ncfile, "temp", timeidx=ALL_TIMES, method="cat"))
    r = to_np(getvar(ncfile, "dbz", timeidx=ALL_TIMES, method="cat"))
    p = to_np(getvar(ncfile, "RAINNC", timeidx=ALL_TIMES, method="cat"))

    stats = {
        'time': times,
        'w_max': np.nanmax(w, axis=(1, 2, 3)),
        'w_p90': np.nanpercentile(w, 90, axis=(1, 2, 3)),
        'w_mean_pos': np.nanmean(np.where(w > 0, w, np.nan), axis=(1, 2, 3)),
        't_max': np.nanmax(t, axis=(1, 2, 3)),
        't_min': np.nanmin(t, axis=(1, 2, 3)),
        'precip_max': np.nanmax(p, axis=(1, 2)),
        'precip_mean': np.nanmean(np.where(p > 0, p, np.nan), axis=(1, 2)),
        'ref_p90': np.nanpercentile(r, 90, axis=(1, 2, 3)),
        'ref_gt10': np.sum(r > 10, axis=(1, 2, 3)),
        'ref_gt50': np.sum(r > 50, axis=(1, 2, 3)),
    }
    return stats

# -------------------- Main Loop --------------------
with open(csv_path, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)
    for i, row in enumerate(reader):
        expid = f"{i:03d}"
        exp_dir = os.path.join(data_root, expid)
        wrf_file = os.path.join(exp_dir, f"wrfout_{expid}.nc")

        if not os.path.exists(wrf_file):
            print(f"[!] Missing file for experiment {expid}")
            continue

        try:
            print(f"[INFO] Processing experiment {expid}...")
            # Ensure the experiment directory exists
            os.makedirs(exp_dir, exist_ok=True)
            # Load the WRF output file
            ncfile = Dataset(wrf_file)
            profile = extract_profile(ncfile)
            ref_cube = extract_reflectivity_cube(ncfile)
            uh1, uh6 = extract_updraft_helicity(ncfile)
            cape, mucape = extract_cape_mucape(ncfile)
            stats = extract_scalar_timeseries(ncfile)

            stats['cape'] = cape
            stats['mucape'] = mucape
            stats['uh_1_1.5'] = uh1
            stats['uh_1_6'] = uh6

            # Save all outputs
            np.savez_compressed(os.path.join(exp_dir, f"summary_{expid}.npz"), **stats)
            np.savez_compressed(os.path.join(exp_dir, f"profile_{expid}.npz"), **profile)
            np.savez_compressed(os.path.join(exp_dir, f"reflectivity_{expid}.npz"), reflectivity=ref_cube)
            np.save(os.path.join(exp_dir, f"parameters_{expid}.npy"), np.array(row, dtype=float))

        except Exception as e:
            print(f"[ERROR] Experiment {expid}: {e}")
            continue
