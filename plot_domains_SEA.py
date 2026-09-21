#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_domains_SEA.py

Description: Plot WRF-Chem nested domains (Domain 1 and Domain 2)
             over Southeast Asia, including monitoring stations.
             
Author: Nittaya (Knit-nit)
"""

import os
import subprocess
import warnings

import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.lines import Line2D
import numpy as np
import pyproj

warnings.filterwarnings('ignore')

# ================================================================
# 1. Environment Setup (Fix PROJ_LIB)
# ================================================================
env_path = '/home/nittaya/miniconda3/envs/ncl_stable'
proj_paths = [
    f'{env_path}/share/proj',
    f'{env_path}/share/proj6',
    '/usr/share/proj',
    '/usr/local/share/proj'
]

proj_found = False
for proj_path in proj_paths:
    if os.path.exists(proj_path):
        os.environ['PROJ_LIB'] = proj_path
        print(f"Set PROJ_LIB to: {proj_path}")
        proj_found = True
        break

if not proj_found:
    try:
        result = subprocess.run(
            ['find', env_path, '-name', 'proj.db', '-type', 'f'],
            capture_output=True, text=True
        )
        if result.stdout:
            proj_db_path = os.path.dirname(result.stdout.strip().split('\n')[0])
            os.environ['PROJ_LIB'] = proj_db_path
            print(f"Found PROJ at: {proj_db_path}")
            proj_found = True
    except Exception as e:
        print(f"Warning: Could not set PROJ_LIB automatically. Error: {e}")

# ================================================================
# 2. Function: Read Monitoring Stations
# ================================================================
def read_monitor_stations(filename):
    """Read monitoring station coordinates from a text file."""
    stations = []
    try:
        possible_paths = [
            filename,
            os.path.join(os.getcwd(), filename),
            os.path.join(os.path.dirname(os.path.abspath(__file__)), filename),
            '/home/nittaya/Build_WRF/WPS-4.4/monitorstation.txt'
        ]
        
        file_found = False
        for path in possible_paths:
            if os.path.exists(path):
                filename = path
                file_found = True
                break
        
        if not file_found:
            print(f"Error: Could not find monitorstation.txt")
            return stations
            
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        lat = float(parts[0])
                        lon = float(parts[1])
                        stations.append((lat, lon))
        print(f"Successfully read {len(stations)} monitoring stations.")
    except Exception as e:
        print(f"Error reading stations: {e}")
    return stations

# ================================================================
# 3. Domain Configuration (Based on namelist.wps)
# ================================================================
ref_lat, ref_lon = 18.744, 99.723
dx_d01, dy_d01 = 12000, 12000
grid_ratio = 3
dx_d02, dy_d02 = dx_d01 / grid_ratio, dy_d01 / grid_ratio

e_we_d01, e_sn_d01 = 160, 160
e_we_d02, e_sn_d02 = 241, 241
i_start_d02, j_start_d02 = 40, 40
deg_per_meter = 1 / 111000

# --- Domain 1 Boundaries ---
width_d01 = dx_d01 * e_we_d01
height_d01 = dy_d01 * e_sn_d01
width_d01_deg = width_d01 * deg_per_meter
height_d01_deg = height_d01 * deg_per_meter
lon_min_d01 = ref_lon - width_d01_deg / 2
lon_max_d01 = ref_lon + width_d01_deg / 2
lat_min_d01 = ref_lat - height_d01_deg / 2
lat_max_d01 = ref_lat + height_d01_deg / 2

# --- Domain 2 Boundaries ---
x0_d02 = (i_start_d02 - 1) * dx_d01
y0_d02 = (j_start_d02 - 1) * dy_d01
width_d02 = dx_d02 * e_we_d02
height_d02 = dy_d02 * e_sn_d02
x_center_d02 = x0_d02 + width_d02 / 2
y_center_d02 = y0_d02 + height_d02 / 2
x_center_d02_deg = (x_center_d02 - width_d01 / 2) * deg_per_meter + ref_lon
y_center_d02_deg = (y_center_d02 - height_d01 / 2) * deg_per_meter + ref_lat
lon_min_d02 = x_center_d02_deg - (width_d02 * deg_per_meter) / 2
lon_max_d02 = x_center_d02_deg + (width_d02 * deg_per_meter) / 2
lat_min_d02 = y_center_d02_deg - (height_d02 * deg_per_meter) / 2
lat_max_d02 = y_center_d02_deg + (height_d02 * deg_per_meter) / 2

print(f"Domain 1 bounds: lon [{lon_min_d01:.2f}, {lon_max_d01:.2f}], lat [{lat_min_d01:.2f}, {lat_max_d01:.2f}]")
print(f"Domain 2 bounds: lon [{lon_min_d02:.2f}, {lon_max_d02:.2f}], lat [{lat_min_d02:.2f}, {lat_max_d02:.2f}]")

# ================================================================
# 4. Load Data
# ================================================================
stations = read_monitor_stations("monitorstation.txt")

shapefile_path = "/home/nittaya/Documents/MapsDomains/Asia/Asia.shp"
use_shapefile = False
gdf_clipped = None

try:
    print(f"Loading shapefile from: {shapefile_path}")
    pyproj.datadir.set_data_dir(os.environ.get('PROJ_LIB', ''))
    
    gdf = gpd.read_file(shapefile_path)
    print(f"Shapefile loaded. CRS: {gdf.crs}")
    
    if gdf.crs is None:
        gdf = gdf.set_crs('EPSG:4326')
    elif str(gdf.crs) != 'EPSG:4326':
        gdf = gdf.to_crs('EPSG:4326')
    
    gdf_clipped = gdf.cx[lon_min_d01:lon_max_d01, lat_min_d01:lat_max_d01]
    print(f"Clipped to {len(gdf_clipped)} features")
    use_shapefile = True
    
except Exception as e:
    print(f"Error loading shapefile: {e}")
    use_shapefile = False

# ================================================================
# 5. Plotting
# ================================================================
crs_latlon = ccrs.PlateCarree()
fig, ax = plt.subplots(figsize=(15, 12), subplot_kw={'projection': crs_latlon})

# Expand extent slightly
margin = 0.5
ax.set_extent([
    lon_min_d01 - margin, lon_max_d01 + margin,
    lat_min_d01 - margin, lat_max_d01 + margin
], crs=crs_latlon)

# --- Plot Shapefile ---
if use_shapefile and gdf_clipped is not None and not gdf_clipped.empty:
    try:
        gdf_clipped.plot(ax=ax, column='ADMIN', cmap='Pastel1', 
                         edgecolor='gray', linewidth=0.8, alpha=0.9, zorder=1)
        
        # Add country labels
        country_colors = {
            'Thailand': '#FFF2CC', 'Laos': '#E2F0D9', 'Vietnam': '#FCE4D6',
            'Cambodia': '#E4DFEC', 'Myanmar': '#DDEBF7', 'China': '#F2F2F2',
            'Bangladesh': '#F8CBAD', 'Bhutan': '#D9E1F2'
        }
        
        for _, row in gdf_clipped.iterrows():
            if not row.geometry.centroid.is_empty:
                country_name = row['ADMIN']
                bg_color = country_colors.get(country_name, 'white')
                
                ax.text(row.geometry.centroid.x, row.geometry.centroid.y, 
                        country_name,
                        fontsize=11, ha='center', va='center', 
                        transform=crs_latlon,
                        fontweight='bold',
                        bbox=dict(facecolor=bg_color, alpha=0.85, 
                                  edgecolor='gray', linewidth=0.5,
                                  boxstyle='round,pad=0.3'), 
                        zorder=3)
    except Exception as e:
        print(f"Error plotting shapefile: {e}")
        use_shapefile = False

# --- Fallback if shapefile fails ---
if not use_shapefile:
    ax.add_feature(cfeature.LAND, facecolor='#f0f0f0', alpha=0.5)
    ax.add_feature(cfeature.OCEAN, facecolor='#e8f4f8', alpha=0.3)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, linestyle='--', edgecolor='gray')
    
    sea_countries = {
        'Thailand': (100.5, 15.5), 'Laos': (102.6, 19.2),
        'Vietnam': (105.5, 21.3), 'Cambodia': (104.9, 12.6),
        'Myanmar': (96.1, 21.0), 'China': (100.0, 24.5),
        'Bangladesh': (90.5, 23.5), 'Bhutan': (90.5, 27.5),
    }
    for name, (lon, lat) in sea_countries.items():
        if (lon_min_d01 - 1 <= lon <= lon_max_d01 + 1 and 
            lat_min_d01 - 1 <= lat <= lat_max_d01 + 1):
            ax.text(lon, lat, name, fontsize=11, ha='center', va='center',
                    transform=crs_latlon, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='white', 
                              alpha=0.8, edgecolor='gray', linewidth=0.5))

# --- Gridlines ---
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                  linewidth=0.5, color='gray', alpha=0.4, linestyle=':')
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 11, 'color': 'black'}
gl.ylabel_style = {'size': 11, 'color': 'black'}

# --- Domain 1 Boundary (Dark Blue) ---
ax.plot([lon_min_d01, lon_max_d01, lon_max_d01, lon_min_d01, lon_min_d01],
        [lat_min_d01, lat_min_d01, lat_max_d01, lat_max_d01, lat_min_d01],
        color='darkblue', linewidth=2.5, linestyle='-', 
        label=f'Domain 1 ({e_we_d01}x{e_sn_d01}), {dx_d01/1000:.0f} km', 
        transform=crs_latlon, zorder=5)

# --- Domain 2 Boundary (Red) ---
ax.plot([lon_min_d02, lon_max_d02, lon_max_d02, lon_min_d02, lon_min_d02],
        [lat_min_d02, lat_min_d02, lat_max_d02, lat_max_d02, lat_min_d02],
        color='red', linewidth=2.5, linestyle='-', 
        label=f'Domain 2 ({e_we_d02}x{e_sn_d02}), {dx_d02/1000:.1f} km', 
        transform=crs_latlon, zorder=6)

# --- Reference Point ---
ax.plot(ref_lon, ref_lat, marker='o', color='black', markersize=7,
        label='Reference Point', transform=crs_latlon, zorder=7)

# --- Monitoring Stations ---
if stations:
    station_lats = [s[0] for s in stations]
    station_lons = [s[1] for s in stations]
    
    ax.scatter(station_lons, station_lats, 
               facecolor='#b3d9ff',
               edgecolor='#0066cc',
               s=60, marker='o', 
               linewidth=1.5,
               label=f'Monitoring Stations (n = {len(stations)})',
               transform=crs_latlon, zorder=8)

# --- Scale Bar (200 km) ---
scale_bar_km = 200
scale_bar_deg = scale_bar_km / 111
scale_bar_x = lon_min_d01 + 0.5
scale_bar_y = lat_min_d01 + 0.5

# Draw main line
ax.plot([scale_bar_x, scale_bar_x + scale_bar_deg], 
        [scale_bar_y, scale_bar_y], 
        color='black', linewidth=3, transform=crs_latlon, zorder=10,
        solid_capstyle='butt')

# Draw end ticks
tick_height = 0.10
ax.plot([scale_bar_x, scale_bar_x], 
        [scale_bar_y - tick_height, scale_bar_y + tick_height], 
        color='black', linewidth=2, transform=crs_latlon, zorder=10)
ax.plot([scale_bar_x + scale_bar_deg, scale_bar_x + scale_bar_deg], 
        [scale_bar_y - tick_height, scale_bar_y + tick_height], 
        color='black', linewidth=2, transform=crs_latlon, zorder=10)

# Draw half tick (100 km)
half_point = scale_bar_x + scale_bar_deg / 2
ax.plot([half_point, half_point], 
        [scale_bar_y - tick_height/2, scale_bar_y + tick_height/2], 
        color='black', linewidth=1.5, transform=crs_latlon, zorder=10)

# Add scale bar text
ax.text(scale_bar_x, scale_bar_y - 0.15, '0', 
        ha='center', va='top', fontsize=10, fontweight='bold',
        transform=crs_latlon, zorder=10)
ax.text(scale_bar_x + scale_bar_deg, scale_bar_y - 0.15, f'{scale_bar_km} km', 
        ha='center', va='top', fontsize=10, fontweight='bold',
        transform=crs_latlon, zorder=10)

# --- North Arrow ---
arrow_x = lon_max_d01 - 0.6
arrow_y = lat_max_d01 - 0.6
arrow_length = 0.8

ax.annotate('', xy=(arrow_x, arrow_y), xytext=(arrow_x, arrow_y - arrow_length),
            arrowprops=dict(arrowstyle='->', color='black', lw=5,
                            shrinkA=0, shrinkB=0,
                            connectionstyle='arc3,rad=0.0'),
            transform=crs_latlon, zorder=10)

ax.text(arrow_x, arrow_y + 0.05, 'N', 
        fontsize=18, fontweight='bold', ha='center', va='bottom',
        transform=crs_latlon, zorder=10)

# --- Legend ---
legend_elements = [
    Line2D([0], [0], color='darkblue', linewidth=2.5, 
           label=f'Domain 1 ({e_we_d01}x{e_sn_d01}), {dx_d01/1000:.0f} km'),
    Line2D([0], [0], color='red', linewidth=2.5, 
           label=f'Domain 2 ({e_we_d02}x{e_sn_d02}), {dx_d02/1000:.1f} km'),
    Line2D([0], [0], marker='o', color='black', markersize=7, 
           linestyle='None', label='Reference Point'),
    Line2D([0], [0], marker='o', color='#0066cc', markersize=7, 
           linestyle='None', markerfacecolor='#b3d9ff', 
           label=f'Monitoring Stations (n = {len(stations)})'),
]

legend = ax.legend(handles=legend_elements, 
                   loc='upper right', 
                   bbox_to_anchor=(0.93, 0.94),
                   fontsize=12, 
                   framealpha=0.95, 
                   edgecolor='black', 
                   facecolor='white', 
                   fancybox=True, 
                   shadow=True, 
                   borderpad=0.6,
                   handlelength=2)

legend.get_frame().set_linewidth(0.9)

# --- Title ---
ax.set_title('WRF-Chem Nested Domains', 
             fontsize=18, fontweight='bold', pad=15)

# --- Adjust layout and Save ---
plt.tight_layout()

output_path = "/home/nittaya/Documents/MapsDomains/wrf_domains_final.jpg"
plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved plot to: {output_path}")

print("\nSummary of plot enhancements:")
print("  - Added gridlines (lat/lon)")
print("  - Scale bar: 200 km, placed in bottom right")
print("  - Legend: placed in upper right (inside Domain 1)")
print("  - North Arrow: placed inside Domain 1")
print("  - Domain 1: Dark blue")
print("  - Domain 2: Red")

plt.show()