
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dpi = 300
inch_to_m = 0.0254

page_height_in = 11.0
page_width_in  = 8.5

page_height_m = page_height_in * inch_to_m
page_width_m  = page_width_in  * inch_to_m

pixel_pitch_m = (1.0 / dpi) * inch_to_m

# Big H size
H_height_m = 0.90 * page_height_m
H_width_m  = 0.60 * page_width_m

y_top    =  H_height_m / 2.0
y_bottom = -H_height_m / 2.0
x_left   = -H_width_m  / 2.0
x_right  =  H_width_m  / 2.0

Ny_H = int(H_height_m / pixel_pitch_m)
Nx_H = int(H_width_m  / pixel_pitch_m)


# Left vertical
left_x = np.full(Ny_H, x_left)
left_y = np.linspace(y_top, y_bottom, Ny_H)

# Right vertical
right_x = np.full(Ny_H, x_right)
right_y = np.linspace(y_top, y_bottom, Ny_H)

# Middle horizontal
mid_x = np.linspace(x_left, x_right, Nx_H)
mid_y = np.zeros(Nx_H)

# Combine for animation order:
xs_H = np.concatenate([left_x, [np.nan], right_x, [np.nan], mid_x])
ys_H = np.concatenate([left_y, [np.nan], right_y, [np.nan], mid_y])

N_droplets_H = len(xs_H)
print(f"Total droplets for H: {N_droplets_H}")

# -------------------------
# Plot setup
# -------------------------
fig, ax = plt.subplots(figsize=(5, 8))

margin_x = 0.05 * H_width_m
margin_y = 0.05 * H_height_m

ax.set_xlim(x_left - margin_x, x_right + margin_x)
ax.set_ylim(y_bottom - margin_y, y_top + margin_y)

ax.set_title("Droplet drawing the letter 'H' at 300 dpi")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

dot,   = ax.plot([], [], "bo", markersize=3)
trail, = ax.plot([], [], "b-", linewidth=1)

def init():
    dot.set_data([], [])
    trail.set_data([], [])
    return dot, trail

def update(frame):
    x = xs_H[frame]
    y = ys_H[frame]

    if np.isnan(x) or np.isnan(y):
        dot.set_data([], [])
    else:
        dot.set_data([x], [y])

    trail.set_data(xs_H[:frame+1], ys_H[:frame+1])
    return dot, trail

ani = FuncAnimation(
    fig,
    update,
    init_func=init,
    frames=N_droplets_H,
    interval=1,
    blit=True,
    repeat=False
)

plt.show()
