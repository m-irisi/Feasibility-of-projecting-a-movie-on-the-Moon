# moon_table.py
"""
Lunar position using the full 60-term Meeus Table 47.A / 47.B.

Reference: Meeus, J. (1998). Astronomical Algorithms, 2nd ed.,
           Chapter 47. Willmann-Bell.
           Chapront-Touze & Chapront (1983). A&A 124, 50.
"""

import math


# ── Step 1: time conversion ───────────────────────────────────────────────────

def calendar_to_T(year, month, day, hour=0, minute=0, second=0):
    """
    Convert a UTC calendar datetime to T (Julian centuries since J2000.0)
    and JDE (Julian Ephemeris Date).

    Applies delta_T = 69 s to convert UTC -> TT, appropriate for 2026.
    Reference: Meeus (1998), Chapter 7.
    """
    day_dec = day + hour/24.0 + minute/1440.0 + second/86400.0
    if month <= 2:
        year -= 1
        month += 12
    A   = int(year / 100)
    B   = 2 - A + int(A / 4)
    JD  = (int(365.25 * (year + 4716))
         + int(30.6001 * (month + 1))
         + day_dec + B - 1524.5)
    JDE = JD + 69.0 / 86400.0          # UTC -> TT
    T   = (JDE - 2451545.0) / 36525.0 # Meeus (1998), eq. 22.1 pg 143
    return T, JDE


# ── Step 2: fundamental angles ────────────────────────────────────────────────

def fundamental_angles(T):
    """
    Compute L', D, M, M', F in degrees (mod 360).
    Reference: Meeus (1998), eq. 47.1-47.5.
    """
    L_prime = (218.3164477 + 481267.88123421*T
          - 0.0015786*T**2
          + T**3/538841
          - T**4/65194000) % 360

    D  = (297.8501921 + 445267.1114034*T
          - 0.0018819*T**2
          + T**3/545868
          - T**4/113065000) % 360

    M  = (357.5291092 + 35999.0502909*T
          - 0.0001536*T**2
          + T**3/24490000) % 360

    M_prime = (134.9633964 + 477198.8675055*T
          + 0.0087414*T**2
          + T**3/69699
          - T**4/14712000) % 360

    F  = (93.2720950 + 483202.0175233*T
          - 0.0036539*T**2
          - T**3/3526000
          + T**4/863310000) % 360

    # Additional arguments for the small additive corrections
    A1 = (119.75 + 131.849*T)   % 360   # Venus perturbation
    A2 = (53.09  + 479264.290*T) % 360  # Jupiter perturbation
    A3 = (313.45 + 481266.484*T) % 360  # latitude correction

    return L_prime, D, M, M_prime, F, A1, A2, A3


# ── Step 3: eccentricity correction factor ────────────────────────────────────

def E_factor(T):
    """
    Correction for terms involving M (Sun's mean anomaly).
    Applied because Earth's orbital eccentricity modulates the Sun's
    perturbing force on the Moon.
    |c_M| = 1 -> multiply term by E
    |c_M| = 2 -> multiply term by E^2
    Reference: Meeus (1998), p.338.
    """
    return 1.0 - 0.002516*T - 0.0000074*T**2


# ── Step 4: Table 47.A — longitude (sl) and distance (sr) ────────────────────
#
# Format: (c_D, c_M, c_Mp, c_F, sl, sr)
# Contribution per row:
#   theta = c_D*D + c_M*M + c_Mp*Mp + c_F*F   (degrees)
#   sigma_l += sl * sin(theta)   [after E correction]
#   sigma_r += sr * cos(theta)   [after E correction]
# Reference: Meeus (1998), Table 47.A, p.339-340.

TABLE_47A = [
    ( 0,  0,  1,  0,  6288774, -20905355),
    ( 2,  0, -1,  0,  1274027,  -3699111),
    ( 2,  0,  0,  0,   658314,  -2955968),
    ( 0,  0,  2,  0,   213618,   -569925),
    ( 0,  1,  0,  0,  -185116,     48888),
    ( 0,  0,  0,  2,  -114332,     -3149),
    ( 2,  0, -2,  0,    58793,    246158),
    ( 2, -1, -1,  0,    57066,   -152138),
    ( 2,  0,  1,  0,    53322,   -170733),
    ( 2, -1,  0,  0,    45758,   -204586),
    ( 0,  1, -1,  0,   -40923,   -129620),
    ( 1,  0,  0,  0,   -34720,    108743),
    ( 0,  1,  1,  0,   -30383,    104755),
    ( 2,  0,  0, -2,    15327,     10321),
    ( 0,  0,  1,  2,   -12528,         0),
    ( 0,  0,  1, -2,    10980,     79661),
    ( 4,  0, -1,  0,    10675,    -34782),
    ( 0,  0,  3,  0,    10034,    -23210),
    ( 4,  0, -2,  0,     8548,    -21636),
    ( 2,  1, -1,  0,    -7888,     24208),
    ( 2,  1,  0,  0,    -6766,     30824),
    ( 1,  0, -1,  0,    -5163,     -8379),
    ( 1,  1,  0,  0,     4987,    -16675),
    ( 2, -1,  1,  0,     4036,    -12831),
    ( 2,  0,  2,  0,     3994,    -10445),
    ( 4,  0,  0,  0,     3861,    -11650),
    ( 2,  0, -3,  0,     3665,     14403),
    ( 0,  1, -2,  0,    -2689,     -7003),
    ( 2,  0, -1,  2,    -2602,         0),
    ( 2, -1, -2,  0,     2390,     10056),
    ( 1,  0,  1,  0,    -2348,      6322),
    ( 2, -2,  0,  0,     2236,     -9884),
    ( 0,  1,  2,  0,    -2120,      5751),
    ( 0,  2,  0,  0,    -2069,         0),
    ( 2, -2, -1,  0,     2048,     -4950),
    ( 2,  0,  1, -2,    -1773,      4130),
    ( 2,  0,  0,  2,    -1595,         0),
    ( 4, -1, -1,  0,     1215,     -3958),
    ( 0,  0,  2,  2,    -1110,         0),
    ( 3,  0, -1,  0,     -892,      3258),
    ( 2,  1,  1,  0,     -810,      2616),
    ( 4, -1, -2,  0,      759,     -1897),
    ( 0,  2, -1,  0,     -713,     -2117),
    ( 2,  2, -1,  0,     -700,      2354),
    ( 2,  1, -2,  0,      691,         0),
    ( 2, -1,  0, -2,      596,         0),
    ( 4,  0,  1,  0,      549,     -1423),
    ( 0,  0,  4,  0,      537,     -1117),
    ( 4, -1,  0,  0,      520,     -1571),
    ( 1,  0, -2,  0,     -487,     -1739),
    ( 2,  1,  0, -2,     -399,         0),
    ( 0,  0,  2, -2,     -381,     -4421),
    ( 1,  1,  1,  0,      351,         0),
    ( 3,  0, -2,  0,     -340,         0),
    ( 4,  0, -3,  0,      330,         0),
    ( 2, -1,  2,  0,      327,         0),
    ( 0,  2,  1,  0,     -323,      1165),
    ( 1,  1, -1,  0,      299,         0),
    ( 2,  0,  3,  0,      294,         0),
    ( 2,  0, -1, -2,        0,      8752),
]


# ── Step 5: Table 47.B — latitude (sb) ───────────────────────────────────────
#
# Format: (c_D, c_M, c_Mp, c_F, sb)
# Contribution per row:
#   theta = c_D*D + c_M*M + c_Mp*Mp + c_F*F   (degrees)
#   sigma_b += sb * sin(theta)   [after E correction]
# Reference: Meeus (1998), Table 47.B, p.341.

TABLE_47B = [
    ( 0,  0,  0,  1,  5128122),
    ( 0,  0,  1,  1,   280602),
    ( 0,  0,  1, -1,   277693),
    ( 2,  0,  0, -1,   173237),
    ( 2,  0, -1,  1,    55413),
    ( 2,  0, -1, -1,    46271),
    ( 2,  0,  0,  1,    32573),
    ( 0,  0,  2,  1,    17198),
    ( 2,  0,  1, -1,     9266),
    ( 0,  0,  2, -1,     8822),
    ( 2, -1,  0, -1,     8216),
    ( 2,  0, -2, -1,     4324),
    ( 2,  0,  1,  1,     4200),
    ( 2,  1,  0, -1,    -3359),
    ( 2, -1, -1,  1,     2463),
    ( 2, -1,  0,  1,     2211),
    ( 2, -1, -1, -1,     2065),
    ( 0,  1, -1, -1,    -1870),
    ( 4,  0, -1, -1,     1828),
    ( 0,  1,  0,  1,    -1794),
    ( 0,  0,  0,  3,    -1749),
    ( 0,  1, -1,  1,    -1565),
    ( 1,  0,  0,  1,    -1491),
    ( 0,  1,  1,  1,    -1475),
    ( 0,  1,  1, -1,    -1410),
    ( 0,  1,  0, -1,    -1344),
    ( 1,  0,  0, -1,    -1335),
    ( 0,  0,  3,  1,     1107),
    ( 4,  0,  0, -1,     1021),
    ( 4,  0, -1,  1,      833),
    ( 0,  0,  1, -3,      777),
    ( 4,  0, -2,  1,      671),
    ( 2,  0,  0, -3,      607),
    ( 2,  0,  2, -1,      596),
    ( 2, -1,  1, -1,      491),
    ( 2,  0, -2,  1,     -451),
    ( 0,  0,  3, -1,      439),
    ( 2,  0,  2,  1,      422),
    ( 2,  0, -3, -1,      421),
    ( 2,  1, -1,  1,     -366),
    ( 2,  1,  0,  1,     -351),
    ( 4,  0,  0,  1,      331),
    ( 2, -1,  1,  1,      315),
    ( 2, -2,  0, -1,      302),
    ( 0,  0,  1,  3,     -283),
    ( 2,  1,  1, -1,     -229),
    ( 1,  1,  0, -1,      223),
    ( 1,  1,  0,  1,      223),
    ( 0,  1, -2, -1,     -220),
    ( 2,  1, -1, -1,     -220),
    ( 1,  0,  1,  1,     -185),
    ( 2, -1, -2, -1,      181),
    ( 0,  1,  2,  1,     -177),
    ( 4,  0, -2, -1,      176),
    ( 4, -1, -1, -1,      166),
    ( 1,  0,  1, -1,     -164),
    ( 4,  0,  1, -1,      132),
    ( 1,  0, -1, -1,     -119),
    ( 4, -1,  0, -1,      115),
    ( 2, -2,  0,  1,      107),
]


# ── Step 6: compute the three sums ───────────────────────────────────────────

def _sin(deg):
    return math.sin(math.radians(deg))

def _cos(deg):
    return math.cos(math.radians(deg))

def compute_sums(T, D, M, Mp, F, A1, A2, A3, Lp):
    """
    Compute sigma_l, sigma_b, sigma_r by iterating over the tables.

    Each row contributes:
      sigma_l += sl * sin(theta) * E_correction
      sigma_r += sr * cos(theta) * E_correction
      sigma_b += sb * sin(theta) * E_correction

    where theta = c_D*D + c_M*M + c_Mp*Mp + c_F*F  (degrees)
    and E_correction = E^|c_M|.

    After the main sums, three small additive corrections are added
    to sigma_l and four to sigma_b (Venus, Jupiter, and nodal terms).
    Reference: Meeus (1998), p.338 and p.342.
    """
    E  = E_factor(T)
    E2 = E * E

    sigma_l = 0.0
    sigma_r = 0.0
    sigma_b = 0.0

    # Table 47.A: longitude and distance
    for (cD, cM, cMp, cF, sl, sr) in TABLE_47A:
        theta = (cD*D + cM*M + cMp*Mp + cF*F) % 360

        # Eccentricity correction based on |c_M|
        abs_cM = abs(cM)
        if abs_cM == 0:
            ecorr = 1.0
        elif abs_cM == 1:
            ecorr = E
        else:                    # abs_cM == 2
            ecorr = E2

        sigma_l += sl * _sin(theta) * ecorr
        sigma_r += sr * _cos(theta) * ecorr

    # Table 47.B: latitude
    for (cD, cM, cMp, cF, sb) in TABLE_47B:
        theta = cD*D + cM*M + cMp*Mp + cF*F

        abs_cM = abs(cM)
        if abs_cM == 0:
            ecorr = 1.0
        elif abs_cM == 1:
            ecorr = E
        else:
            ecorr = E2

        sigma_b += sb * _sin(theta) * ecorr

    # Small additive corrections to sigma_l (Meeus 1998, p.342)
    sigma_l += 3958 * _sin(A1)
    sigma_l += 1962 * _sin(Lp - F)
    sigma_l +=  318 * _sin(A2)

    # Small additive corrections to sigma_b (Meeus 1998, p.342)
    sigma_b -= 2235 * _sin(Lp)
    sigma_b +=  382 * _sin(A3)
    sigma_b +=  175 * _sin(A1 - F)
    sigma_b +=  175 * _sin(A1 + F)
    sigma_b +=  127 * _sin(Lp - Mp)
    sigma_b -=  115 * _sin(Lp + Mp)

    return sigma_l, sigma_b, sigma_r


# ── Step 7: final coordinates ─────────────────────────────────────────────────

def moon_position(year, month, day, hour=0, minute=0, second=0):
    """
    Full pipeline: UTC calendar date -> ecliptic coordinates.

    Returns
    -------
    dict:
        T        : Julian centuries since J2000.0
        JDE      : Julian Ephemeris Date
        lambda_  : geocentric ecliptic longitude (degrees)
        beta     : geocentric ecliptic latitude (degrees)
        Delta    : Earth-Moon distance (km)
    """
    T, JDE = calendar_to_T(year, month, day, hour, minute, second)
    Lp, D, M, Mp, F, A1, A2, A3 = fundamental_angles(T)

    sl, sb, sr = compute_sums(T, D, M, Mp, F, A1, A2, A3, Lp)

    # Equations 47.1, 47.2, 47.3 from Meeus (1998), p.342
    lam   = (Lp + sl / 1e6) % 360   # degrees
    beta  = sb / 1e6                  # degrees
    Delta = 385000.56 + sr / 1000.0  # km

    return {
        "T":       T,
        "JDE":     JDE,
        "lambda_": lam,
        "beta":    beta,
        "Delta":   Delta,
    }

