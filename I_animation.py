
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dpi = 300
inch_to_m = 0.0254
page_height_in = 11.0

pixel_pitch_m = (1.0 / dpi) * inch_to_m            # meters per pixel
N_pixels = int(page_height_in * dpi)               

# y-positions from bottom to top of page 
start_y = -(page_height_in / 2.0) * inch_to_m      # center page at y=0
ys = start_y + pixel_pitch_m * np.arange(N_pixels)
xs = np.zeros_like(ys)                             

# -------------------------
# Set up the figure
# -------------------------
fig, ax = plt.subplots(figsize=(3, 8))             
ax.set_xlim(-0.002, 0.002)                         
ax.set_ylim(ys.min(), ys.max())
ax.set_title("Droplet drawing the vertical line 'I' at 300 dpi")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

dot, = ax.plot([], [], "bo", markersize=3)
trail, = ax.plot([], [], "b-", linewidth=1)

def init():
    dot.set_data([], [])
    trail.set_data([], [])
    return dot, trail

def update(frame):
    # current droplet position
    x = xs[frame]
    y = ys[frame]
    dot.set_data([x], [y])               

    # trail of all previous droplets
    trail.set_data(xs[:frame+1], ys[:frame+1])

    return dot, trail

ani = FuncAnimation(
    fig,
    update,
    init_func=init,
    frames=N_pixels,
    interval=1,          # visual speed (ms between frames)
    blit=True,
    repeat=False
)

plt.show()
