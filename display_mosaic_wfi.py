"""Display a mosaic of all 18 Roman WFI SCAs positioned by true focal-plane coordinates."""

import keyring
import keyring.backends.null
keyring.set_keyring(keyring.backends.null.Keyring())

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib.pyplot as plt
import roman_datamodels as rdm
from astropy.visualization import ZScaleInterval

from streaming_utils import get_MAST_token, group_file_urls, stream_file_group_to_buffer, close_buffer_streams

# --- Configuration -----------------------------------------------------------
visit_id   = '0012401001001002001'
exp_num    = 1
data_level = 2
output_png = 'mosaic_all_scas.png'

# --- Focal-plane geometry (mm) -----------------------------------------------
# x/y are SCA centre positions; r is rotation (180 = flipped relative to r=0).
SCA_INFO = {
    1 : {'x' : -22.1400, 'y' :  12.1500, 'r' : 180},
    2 : {'x' : -22.2900, 'y' : -37.0300, 'r' : 180},
    3 : {'x' : -22.4400, 'y' : -82.0600, 'r' :   0},
    4 : {'x' : -66.4200, 'y' :  20.9000, 'r' : 180},
    5 : {'x' : -66.9200, 'y' : -28.2800, 'r' : 180},
    6 : {'x' : -67.4200, 'y' : -73.0600, 'r' :   0},
    7 : {'x' :-110.7000, 'y' :  42.2000, 'r' : 180},
    8 : {'x' :-111.4800, 'y' :  -6.9800, 'r' : 180},
    9 : {'x' :-112.6400, 'y' : -51.0600, 'r' :   0},
   10 : {'x' :  22.1400, 'y' :  12.1500, 'r' : 180},
   11 : {'x' :  22.2900, 'y' : -37.0300, 'r' : 180},
   12 : {'x' :  22.4400, 'y' : -82.0600, 'r' :   0},
   13 : {'x' :  66.4200, 'y' :  20.9000, 'r' : 180},
   14 : {'x' :  66.9200, 'y' : -28.2800, 'r' : 180},
   15 : {'x' :  67.4200, 'y' : -73.0600, 'r' :   0},
   16 : {'x' : 110.7000, 'y' :  42.2000, 'r' : 180},
   17 : {'x' : 111.4800, 'y' :  -6.9800, 'r' : 180},
   18 : {'x' : 112.6400, 'y' : -51.0600, 'r' :   0},
}

# --- Stream data -------------------------------------------------------------
mast_token  = get_MAST_token()
urls        = group_file_urls(mast_token, visit_id, exp_num='', data_level=data_level)
buffer_dict = stream_file_group_to_buffer(urls, exp_num=exp_num, mast_token=mast_token)

# --- Derive figure layout from focal-plane coordinates -----------------------
# Each SCA is ~40.88 mm × 40.88 mm (4088 px × 4088 px, 10 µm/px).
SCA_SIZE_MM = 40.88

xs = np.array([SCA_INFO[s]['x'] for s in SCA_INFO])
ys = np.array([SCA_INFO[s]['y'] for s in SCA_INFO])

# Bounding box of SCA centres plus half a chip on each side
x_min = xs.min() - SCA_SIZE_MM / 2
x_max = xs.max() + SCA_SIZE_MM / 2
y_min = ys.min() - SCA_SIZE_MM / 2
y_max = ys.max() + SCA_SIZE_MM / 2

fp_width  = x_max - x_min   # mm
fp_height = y_max - y_min   # mm

# Figure size in inches; keep aspect ratio
fig_width  = 18.0
fig_height = fig_width * fp_height / fp_width

fig = plt.figure(figsize=(fig_width, fig_height), facecolor='black')

zscale = ZScaleInterval()

# --- Place each SCA as a manually positioned Axes ----------------------------
# fig.add_axes takes [left, bottom, width, height] in figure-fraction units.

def fp_to_fig(x_mm, y_mm):
    """Convert focal-plane mm coordinates to figure-fraction (left, bottom)."""
    left   = (x_mm - x_min) / fp_width
    bottom = (y_mm - y_min) / fp_height
    return left, bottom

ax_w = SCA_SIZE_MM / fp_width   # axes width  in figure fraction
ax_h = SCA_SIZE_MM / fp_height  # axes height in figure fraction

# Small gap between chips (≈2 % of chip size) so borders are visible
gap = 0.02
ax_w_inner = ax_w * (1 - gap)
ax_h_inner = ax_h * (1 - gap)

for scanum, info in SCA_INFO.items():
    # Centre of this SCA in figure fraction
    cx, cy = fp_to_fig(info['x'], info['y'])

    left   = cx - ax_w_inner / 2
    bottom = cy - ax_h_inner / 2

    ax = fig.add_axes([left, bottom, ax_w_inner, ax_h_inner])
    ax.set_facecolor('black')

    buf = buffer_dict.get(scanum)
    if buf is not None:
        dm   = rdm.open(buf)
        data = dm.data[...].astype(float)

        # Rotate 180° for SCAs where r=180
        # this is skipped because the data is already in the correct orientation for display
        #if info['r'] == 180:
        #    data = data[::-1, ::-1]

        vmin, vmax = zscale.get_limits(data)
        ax.imshow(data, origin='lower', cmap='gray',
                  vmin=vmin, vmax=vmax, interpolation='nearest')

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor('0.4')
        spine.set_linewidth(0.5)

    ax.text(0.5, 0.97, f'SCA {scanum:02d}',
            transform=ax.transAxes,
            ha='center', va='top',
            fontsize=7, color='white')

fig.text(0.5, 0.995,
         f'Roman WFI — All 18 SCAs  |  Visit {visit_id}, Exposure {exp_num}',
         ha='center', va='top', fontsize=13, color='white')

plt.savefig(output_png, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor())
print(f'Saved {output_png}')
plt.show()

# --- Clean up ----------------------------------------------------------------
close_buffer_streams(buffer_dict)
