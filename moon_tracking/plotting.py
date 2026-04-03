"""
Functions to aid with plotting
"""

import numpy as np
from math import radians, cos, sin
from moon_position import moon_position, calendar_to_T

def ecliptic_to_cartesian(Delta, lambda_, beta):
    """Convert ecliptic coordinates to Cartesian coordinates.

    Parameters:
        Delta (array): Distance to the Moon in kilometers.
        lambda_ (array): Ecliptic longitude in degrees.
        beta (array): Ecliptic latitude in degrees.

    Returns:
        tuple: (x, y, z) Cartesian coordinates in kilometers.
    """
    lambda_rad = np.radians(lambda_)
    beta_rad = np.radians(beta)

    x = Delta * np.cos(beta_rad) * np.cos(lambda_rad)
    y = Delta * np.cos(beta_rad) * np.sin(lambda_rad)
    z = Delta * np.sin(beta_rad)

    return x, y, z


def orbit_cartesian(start_year, start_month, start_day,
                    n_days=89, step_hours=6):
    """
    Sample the Moon's position every step_hours for n_days days,
    returning Cartesian coordinates in km.

    Parameters
    ----------
    start_year, start_month, start_day : int
    n_days      : float   total duration in days (~89 for 3 synodic months)
    step_hours  : float   time step in hours

    Returns
    -------
    x, y, z : np.ndarray   Cartesian coordinates in km (Earth at origin)
    t_days  : np.ndarray   elapsed days from start
    lam     : np.ndarray   ecliptic longitude in degrees
    beta    : np.ndarray   ecliptic latitude in degrees
    delta   : np.ndarray   distance in km
    """
    # Build list of sample times
    total_steps = int(n_days * 24 / step_hours)
    t_days = np.arange(total_steps) * (step_hours / 24.0)

    lam_arr   = np.zeros(total_steps)
    beta_arr  = np.zeros(total_steps)
    delta_arr = np.zeros(total_steps)

    from datetime import datetime, timedelta
    start_dt = datetime(start_year, start_month, start_day, 0, 0, 0)

    for i, d in enumerate(t_days):
        dt = start_dt + timedelta(days=float(d))
        pos = moon_position(dt.year, dt.month, dt.day,
                            dt.hour, dt.minute, dt.second)
        lam_arr[i]   = pos["lambda_"]
        beta_arr[i]  = pos["beta"]
        delta_arr[i] = pos["Delta"]

    # Convert spherical ecliptic -> Cartesian (km)
    x, y, z = ecliptic_to_cartesian(delta_arr, lam_arr, beta_arr)

    return x, y, z, t_days, lam_arr, beta_arr, delta_arr