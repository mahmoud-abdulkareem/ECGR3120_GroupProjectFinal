#!/usr/bin/env python3

# Using Python to simulate the droplet printer project

from __future__ import annotations
import math
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Given parameters 
# ---------------------------
# These are all the values the project gives us (velocities, distances, etc.)
Vx_base = 20.0              # horizontal droplet speed (m/s)
D  = 3e-3                    # distance from gun to paper
L1_base = 0.5e-3            # capacitor length
L2_base = 1.25e-3           # distance from cap exit to paper
W_base  = 1e-3              # gap between capacitor plates
rho = 1000.0                # density (close to water)
d_base   = 84e-6            # droplet diameter
q_base   = -1.9e-10         # charge of droplet
dpi = 300                   # printer resolution
page_height_in = 11.0       # page height
page_width_in  = 8.5        # page width

inch_to_m = 0.0254          # converting inches to meters

# ---------------------------
# Helper functions
# ---------------------------

def droplet_mass(d_m: float, rho_kg_m3: float) -> float:
    """
    Calculates droplet mass assuming it's a sphere.
    """
    return (4.0/3.0) * math.pi * (d_m/2.0)**3 * rho_kg_m3

def compute_times(Vx_m_s: float, D_m: float, L1_m: float, L2_m: float):
    """
    Computes the 3 times we need:
    - time inside the capacitor 
    - drift time after leaving capacitor 
    - total time to reach the paper
    """
    T_inside = L1_m / Vx_m_s
    t2 = L2_m / Vx_m_s
    t_center = D_m / Vx_m_s
    return T_inside, t2, t_center

def voltage_gain_K(m_kg: float, q_C: float, W_m: float, T_inside_s: float, t2_s: float) -> float:

    return m_kg * W_m / (q_C * T_inside_s * t2_s)

def build_vertical_staircase(K: float, N_pixels: int, pixel_pitch_m: float,
                             T_inside: float, t_center: float):
    
    start_y = -(page_height_in/2.0) * inch_to_m
    ys = start_y + pixel_pitch_m * np.arange(N_pixels)
    Vs = K * ys
    ts = t_center + T_inside * np.arange(N_pixels)
    return ts, Vs

# ---------------------------
# Parts 1–3 calculations
# ---------------------------
m_base = droplet_mass(d_base, rho)
T_inside_base, t2_base, t_center_base = compute_times(Vx_base, D, L1_base, L2_base)

# For Part 1
Q1_seconds = t_center_base

# For Part 2: number of droplets = 3300 pixels (11 inches * 300 dpi)
N_pixels = int(page_height_in * dpi)
period_per_drop = T_inside_base
time_to_draw_I = N_pixels * period_per_drop
Q2_total_seconds = t_center_base + time_to_draw_I

# For Part 3: compute K and then the full staircase voltage
K_base = voltage_gain_K(m_base, q_base, W_base, T_inside_base, t2_base)
pixel_pitch_m = (1.0 / dpi) * inch_to_m
dV_per_pixel_base = K_base * pixel_pitch_m
ts_base, Vs_base = build_vertical_staircase(K_base, N_pixels, pixel_pitch_m,
                                            T_inside_base, t_center_base)

# Print results for the first 3 parts
print("\n=== RESULTS (Parts 1–3) ===")
print(f"Droplet mass: {m_base:.3e} kg")
print(f"T_inside (cap time): {T_inside_base*1e6:.2f} µs")
print(f"t2 (drift time): {t2_base*1e6:.2f} µs")
print(f"Time to center (no voltage): {Q1_seconds*1e6:.2f} µs")
print(f"Time to draw 11\" line: {Q2_total_seconds*1e3:.2f} ms")
print(f"Voltage step per pixel: {dV_per_pixel_base:.2f} V")
print(f"Voltage range for full 11\": [{Vs_base.min():.0f}, {Vs_base.max():.0f}] V")

# Plot staircase for Parts 1–3
plt.figure(figsize=(10, 4))
plt.step(ts_base, Vs_base, where="post")
plt.xlabel("Time (s)")
plt.ylabel("Voltage V(t)")
plt.title("Staircase V(t) to draw 11\" vertical line at 300 dpi (base case)")
plt.tight_layout()
plt.show()

# ---------------------------
# PART 4 — Changing parameters
# ---------------------------
print("\n=== PART 4: Parameter variations ===")

def simulate_vertical_I_scenario(label, Vx, L1, L2, W, d, q):
    """
    Runs the same simulation as before but with the modified parameter.
    Prints results so we can compare each case.
    """
    m = droplet_mass(d, rho)
    T_inside, t2, t_center = compute_times(Vx, D, L1, L2)
    K = voltage_gain_K(m, q, W, T_inside, t2)
    dV_pixel = K * pixel_pitch_m
    total_time = t_center + N_pixels * T_inside
    ts, Vs = build_vertical_staircase(K, N_pixels, pixel_pitch_m, T_inside, t_center)

    print(f"\nScenario {label}:")
    print(f"T_inside = {T_inside*1e6:.2f} µs")
    print(f"t2 = {t2*1e6:.2f} µs")
    print(f"Total draw time = {total_time*1e3:.2f} ms")
    print(f"dV per pixel = {dV_pixel:.2f} V")
    print(f"Voltage range = [{Vs.min():.0f}, {Vs.max():.0f}] V")

    return ts, Vs, label

# Running all 5 cases
scenarios = []
scenarios.append(simulate_vertical_I_scenario("4(a): L2 -> 3L2",
    Vx_base, L1_base, 3*L2_base, W_base, d_base, q_base))

scenarios.append(simulate_vertical_I_scenario("4(b): L1 -> 2L1",
    Vx_base, 2*L1_base, L2_base, W_base, d_base, q_base))

scenarios.append(simulate_vertical_I_scenario("4(c): d -> 10d",
    Vx_base, L1_base, L2_base, W_base, 10*d_base, q_base))

scenarios.append(simulate_vertical_I_scenario("4(d): Vx -> 2Vx",
    2*Vx_base, L1_base, L2_base, W_base, d_base, q_base))

scenarios.append(simulate_vertical_I_scenario("4(e): q -> 5q",
    Vx_base, L1_base, L2_base, W_base, d_base, 5*q_base))

# Plot all scenarios on one graph
plt.figure(figsize=(10, 5))
for ts, Vs, lbl in scenarios:
    plt.step(ts, Vs, where="post", label=lbl)

plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.title("Part 4: Staircase voltages for different parameter changes")
plt.legend()
plt.tight_layout()
plt.show()

# ---------------------------
# PART 5 — Drawing the letter H
# ---------------------------
print("\n=== PART 5: Drawing H with 2 capacitors ===")

# Using same K for x and y to keep things simple
K_y = K_base
K_x = K_base

# Making a pretty big H but still reasonable voltage
page_height_m = page_height_in * inch_to_m
page_width_m  = page_width_in * inch_to_m
H_height_m = 0.90 * page_height_m
H_width_m  = 0.60 * page_width_m

# Coordinates of the big H
y_top =  H_height_m / 2
y_bottom = -H_height_m / 2
x_left = -H_width_m / 2
x_right =  H_width_m / 2

Ny_H = int(H_height_m / pixel_pitch_m)
Nx_H = int(H_width_m  / pixel_pitch_m)

# Building the H shape droplet by droplet
points = []

# left vertical bar
for k in range(Ny_H):
    points.append((x_left, y_top - k*pixel_pitch_m))

# right vertical bar
for k in range(Ny_H):
    points.append((x_right, y_top - k*pixel_pitch_m))

# middle horizontal bar
for k in range(Nx_H):
    points.append((x_left + k*pixel_pitch_m, 0.0))

points = np.array(points)
xs_H = points[:,0]
ys_H = points[:,1]

# Convert (x, y) → needed Vx, Vy
Vy_t = K_y * ys_H
Vx_t = K_x * xs_H

N_droplets_H = len(points)
T_inside_H = T_inside_base
total_time_H = t_center_base + N_droplets_H * T_inside_H
ts_H = t_center_base + T_inside_H * np.arange(N_droplets_H)

print(f"H height = {H_height_m:.3f} m")
print(f"H width  = {H_width_m:.3f} m")
print(f"Total droplets for H: {N_droplets_H}")
print(f"Total draw time = {total_time_H*1e3:.2f} ms")
print(f"Vy range: [{Vy_t.min():.0f}, {Vy_t.max():.0f}] V")
print(f"Vx range: [{Vx_t.min():.0f}, {Vx_t.max():.0f}] V")

# Plot both voltages for the H
plt.figure(figsize=(10,5))
plt.step(ts_H, Vy_t, where="post", label="Vertical capacitor V_y(t)")
plt.step(ts_H, Vx_t, where="post", label="Horizontal capacitor V_x(t)")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.title("Part 5: Staircase voltages to draw big 'H' at 300 dpi")
plt.legend()
plt.tight_layout()
plt.show()
