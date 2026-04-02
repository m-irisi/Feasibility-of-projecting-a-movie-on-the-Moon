from moon_position import moon_position
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

# ── Validation against Meeus (1998) worked example, p.342 ────────────────────
if __name__ == "__main__":
    # Meeus example: 1992 April 12, 0h TT
    # Expected: lambda = 133.167°, beta = -3.229°, Delta = 368409 km
    pos = moon_position(1992, 4, 12, 0, 0, 0)
    print("Validation — Meeus (1998) p.342 worked example:")
    print(f"  lambda = {pos['lambda_']:.3f}°  (expected 133.167°)")
    print(f"  beta   = {pos['beta']:.3f}°  (expected -3.229°)")
    print(f"  Delta  = {pos['Delta']:.0f} km  (expected 368409 km)")
    print()
    pos2 = moon_position(2026, 6, 29, 0, 58, 0)
    print("Full Moon, 29 June 2026 00:58 BST (= 23:58 UTC 28 June):")
    print(f"  lambda = {pos2['lambda_']:.3f}°")
    print(f"  beta   = {pos2['beta']:.3f}°")
    print(f"  Delta  = {pos2['Delta']:.0f} km")