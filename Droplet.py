import streamlit as st
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from fpdf import FPDF # Import FPDF for PDF generation
import io # For handling in-memory image data
import pandas as pd # Import pandas for data table
import math # Import math for new calculations
import Flowregime as fr # Import the new module

#@st.cache_data
def create_plot_buffers(plot_data_original, plot_data_adjusted, plot_data_after_gravity, plot_data_after_mist_extractor):
    import matplotlib.pyplot as plt
    import io

    micron_unit_label = "µm"
    buf_original = io.BytesIO()
    buf_adjusted = io.BytesIO()
    buf_after_gravity = io.BytesIO()
    buf_after_me = io.BytesIO()

    # --- Original ---
    dp_values_microns_original = plot_data_original['dp_values_ft'] * FT_TO_MICRON
    fig_original, ax_original = plt.subplots(figsize=(10, 6))
    ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
   #ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
    ax_original.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
    ax_original.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
    ax_original.tick_params(axis='y', labelcolor='black')
    ax_original.set_ylim(0, 1.2)
    ax_original.set_xlim(0, max(dp_values_microns_original) * 1.0 if dp_values_microns_original.size > 0 else 1000)
    ax2_original = ax_original.twinx()
    ax2_original.plot(dp_values_microns_original, plot_data_original['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
    ax2_original.set_ylabel('Volume Fraction', color='black', fontsize=12)
    ax2_original.tick_params(axis='y', labelcolor='black')
    max_norm_fv_original = max(plot_data_original['volume_fraction']) if plot_data_original['volume_fraction'].size > 0 else 0.1
    ax2_original.set_ylim(0, max_norm_fv_original * 1.2)
    plt.title('Entrainment Droplet Size Distribution (Before Inlet Device)', fontsize=14)
    lines_original, labels_original = ax_original.get_legend_handles_labels()
    lines2_original, labels2_original = ax2_original.get_legend_handles_labels()
    ax2_original.legend(lines_original + lines2_original, labels_original + labels2_original, loc='upper left', fontsize=10)
    ax_original.axhline(y=1, color='b', linestyle='--')
    plt.grid(True, linestyle='--', alpha=0.7)
    fig_original.savefig(buf_original, format="png", dpi=300)
    buf_original.seek(0)
    plt.close(fig_original)

    # --- Adjusted ---
    dp_values_microns_adjusted = plot_data_adjusted['dp_values_ft'] * FT_TO_MICRON
    fig_adjusted, ax_adjusted = plt.subplots(figsize=(10, 6))
    ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
   #ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
    ax_adjusted.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
    ax_adjusted.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
    ax_adjusted.tick_params(axis='y', labelcolor='black')
    ax_adjusted.set_ylim(0, 1.2)
    ax_adjusted.set_xlim(0, max(dp_values_microns_adjusted) * 1.0 if dp_values_microns_adjusted.size > 0 else 1000)
    ax2_adjusted = ax_adjusted.twinx()
    ax2_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
    ax2_adjusted.set_ylabel('Volume Fraction', color='black', fontsize=12)
    ax2_adjusted.tick_params(axis='y', labelcolor='black')
    max_norm_fv_adjusted = max(plot_data_adjusted['volume_fraction']) if plot_data_adjusted['volume_fraction'].size > 0 else 0.1
    ax2_adjusted.set_ylim(0, max_norm_fv_adjusted * 1.2)
    plt.title('Entrainment Droplet Size Distribution (After Inlet Device)', fontsize=14)
    lines_adjusted, labels_adjusted = ax_adjusted.get_legend_handles_labels()
    lines2_adjusted, labels2_adjusted = ax2_adjusted.get_legend_handles_labels()
    ax2_adjusted.legend(lines_adjusted + lines2_adjusted, labels_adjusted + labels2_adjusted, loc='upper left', fontsize=10)
    ax_adjusted.axhline(y=1, color='b', linestyle='--')
    plt.grid(True, linestyle='--', alpha=0.7)
    fig_adjusted.savefig(buf_adjusted, format="png", dpi=300)
    buf_adjusted.seek(0)
    plt.close(fig_adjusted)

    # --- After Gravity ---
    dp_values_microns_after_gravity = plot_data_after_gravity['dp_values_ft'] * FT_TO_MICRON
    fig_after_gravity, ax_after_gravity = plt.subplots(figsize=(10, 6))
    ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
   #ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
    ax_after_gravity.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
    ax_after_gravity.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
    ax_after_gravity.tick_params(axis='y', labelcolor='black')
    ax_after_gravity.set_ylim(0, 1.2)
    ax_after_gravity.set_xlim(0, max(dp_values_microns_after_gravity) * 1.0 if dp_values_microns_after_gravity.size > 0 else 1000)
    ax2_after_gravity = ax_after_gravity.twinx()
    ax2_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
    ax2_after_gravity.set_ylabel('Volume Fraction', color='black', fontsize=12)
    ax2_after_gravity.tick_params(axis='y', labelcolor='black')
    max_norm_fv_after_gravity = max(plot_data_after_gravity['volume_fraction']) if plot_data_after_gravity['volume_fraction'].size > 0 else 0.1
    ax2_after_gravity.set_ylim(0, max_norm_fv_after_gravity * 1.2)
    plt.title('Entrainment Droplet Size Distribution (After Gravity Settling)', fontsize=14)
    lines_after_gravity, labels_after_gravity = ax_after_gravity.get_legend_handles_labels()
    lines2_after_gravity, labels2_after_gravity = ax2_after_gravity.get_legend_handles_labels()
    ax2_after_gravity.legend(lines_after_gravity + lines2_after_gravity, labels_after_gravity + labels2_after_gravity, loc='upper left', fontsize=10)
    ax_after_gravity.axhline(y=1, color='b', linestyle='--')
    plt.grid(True, linestyle='--', alpha=0.7)
    fig_after_gravity.savefig(buf_after_gravity, format="png", dpi=300)
    buf_after_gravity.seek(0)
    plt.close(fig_after_gravity)

    # --- After Mist Extractor ---
    dp_values_microns_after_me = plot_data_after_mist_extractor['dp_values_ft'] * FT_TO_MICRON
    fig_after_me, ax_after_me = plt.subplots(figsize=(10, 6))
    ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
   #ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
    ax_after_me.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
    ax_after_me.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
    ax_after_me.tick_params(axis='y', labelcolor='black')
    ax_after_me.set_ylim(0, 1.2)
    ax_after_me.set_xlim(0, max(dp_values_microns_after_me) * 1.0 if dp_values_microns_after_me.size > 0 else 1000)
    ax2_after_me = ax_after_me.twinx()
    ax2_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
    ax2_after_me.set_ylabel('Volume Fraction', color='black', fontsize=12)
    ax2_after_me.tick_params(axis='y', labelcolor='black')
    max_norm_fv_after_me = max(plot_data_after_mist_extractor['volume_fraction']) if plot_data_after_mist_extractor['volume_fraction'].size > 0 else 0.1
    ax2_after_me.set_ylim(0, max_norm_fv_after_me * 1.2)
    plt.title('Entrainment Droplet Size Distribution (After Mist Extractor)', fontsize=14)
    lines_after_me, labels_after_me = ax_after_me.get_legend_handles_labels()
    lines2_after_me, labels2_after_me = ax2_after_me.get_legend_handles_labels()
    ax2_after_me.legend(lines_after_me + lines2_after_me, labels_after_me + labels2_after_me, loc='upper left', fontsize=10)
    ax_after_me.axhline(y=1, color='b', linestyle='--')
    plt.grid(True, linestyle='--', alpha=0.7)
    fig_after_me.savefig(buf_after_me, format="png", dpi=300)
    buf_after_me.seek(0)
    plt.close(fig_after_me)

    return buf_original, buf_adjusted, buf_after_gravity, buf_after_me

FT_TO_M = 0.3048  # 1 ft = 0.3048 m

def plot_mesh_pad_efficiency_with_pressure(Q_gas_ft3_s, A_installed_ft2,
                                           rho_l_fps, rho_g_fps, mu_g_fps,
                                           mesh_pad_type_params_fps, K_factor, K_dp, results):
    # Actual face velocity in ft/s
    V_face_actual_ft_s = Q_gas_ft3_s / A_installed_ft2

    # Allowable velocity in ft/s
    V_allow_ft_s = K_factor * np.sqrt((rho_l_fps - rho_g_fps) / rho_g_fps)

    # Velocity range for plotting (ft/s)
    V_max_ft_s = max(V_face_actual_ft_s, V_allow_ft_s) * 2
    velocities_ft_s = np.linspace(0.05, V_max_ft_s, 500)

    dp_microns = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11, 13, 15,17,19,20])

    def calculate_single_wire_efficiency(Stk):
        if Stk <= 0:
            return 0.0
        numerator = -0.105 + 0.995 * (Stk**1.00493)
        denominator = 0.6261 + (Stk**1.00493)
        if denominator == 0:
            return 0.0
        Ew = numerator / denominator
        return max(0.0, min(1.0, Ew))

    def mesh_pad_efficiency_func(dp_fps, V_face_ft_s, rho_l_fps, rho_g_fps, mu_g_fps, mesh_pad_type_params_fps):
        if V_face_ft_s <= 0 or mu_g_fps <= 0 or dp_fps <= 0:
            return 0.0
        Dw_fps = mesh_pad_type_params_fps["wire_diameter_in"] / 12.0
        pad_thickness_fps = mesh_pad_type_params_fps["thickness_in"] / 12.0
        specific_surface_area_fps = mesh_pad_type_params_fps["specific_surface_area_ft2_ft3"]
        if Dw_fps == 0:
            return 0.0
        Stk = ((rho_l_fps - rho_g_fps) * (dp_fps**2) * V_face_ft_s) / (18 * mu_g_fps * Dw_fps)
        Ew = calculate_single_wire_efficiency(Stk)
        exponent = -0.238 * specific_surface_area_fps * pad_thickness_fps * Ew
        E_pad = 1 - np.exp(exponent)
        return max(0.0, min(1.0, E_pad))

    def mesh_pad_pressure_drop(V_face_ft_s, rho_g_fps, voidage_frac, K_dp, pad_thickness_ft, Dw_fps):
        # Pressure drop per unit thickness (Pa per ft)
        pressure_drop_per_ft = (((1 - voidage_frac)) * (rho_g_fps * V_face_ft_s**2) * K_dp * 47.8803) / Dw_fps
        return pressure_drop_per_ft * pad_thickness_ft  # Total pressure drop in Pa

    fig, ax1 = plt.subplots(figsize=(10, 6))

    Dw_fps = mesh_pad_type_params_fps["wire_diameter_in"] / 12.0
    pad_thickness_ft = mesh_pad_type_params_fps["thickness_in"] / 12.0
    voidage_frac = mesh_pad_type_params_fps.get("voidage_percent", 98.5) / 100.0

    for dp_um in dp_microns:
        # Convert microns to feet
        dp_ft = dp_um * 1e-6 / 0.3048
        efficiencies = [mesh_pad_efficiency_func(dp_ft, v, rho_l_fps, rho_g_fps, mu_g_fps, mesh_pad_type_params_fps) * 100
                        for v in velocities_ft_s]
        ax1.plot(velocities_ft_s, efficiencies, color=ax1._get_lines.get_next_color(), linestyle='-', linewidth=1)  # no label for line
        E_actual = mesh_pad_efficiency_func(dp_ft, V_face_actual_ft_s, rho_l_fps, rho_g_fps, mu_g_fps, mesh_pad_type_params_fps) * 100
        ax1.plot(V_face_actual_ft_s, E_actual, 'o', markersize=5, label=f"Droplet {dp_um} μm Capture: {E_actual:.2f}%")

    ax1.set_xlabel("Face Velocity through Mesh Pad (ft/s)") 
    ax1.set_ylabel("Capture Efficiency (%)")
    ax1.grid(True)

    ax1.axvline(V_allow_ft_s, color='red', linestyle='--')
    ax1.axvline(V_allow_ft_s * 0.3, color='red', linestyle='--')

    # Shade operating region between 0.3*V_allow_ft_s and V_allow_ft_s
    ax1.axvspan(V_allow_ft_s * 0.3, V_allow_ft_s, color='green', alpha=0.2, label='Operating Region')

    # Shade reentrainment region beyond V_allow_ft_s to max x-axis limit
    ax1.axvspan(V_allow_ft_s, velocities_ft_s[-1], color='red', alpha=0.1, label='Flooding Region')

    # Add text labels (adjust y position as needed)
    ylim = ax1.get_ylim()
    y_pos = ylim[1] * 0.5    
    ax1.text((V_allow_ft_s * 0.3 + V_allow_ft_s) / 2, y_pos, f'Mesh Pad Operating Region\n Ks = {K_factor} ft/s ', color='green', ha='center', fontsize=9)
    ax1.text((V_allow_ft_s + velocities_ft_s[-1]) / 2, y_pos*1.5, 'Flooding Region', color='red', ha='center', fontsize=9)


    ax2 = ax1.twinx()
    dp_pressures = [mesh_pad_pressure_drop(v, rho_g_fps, voidage_frac, K_dp, pad_thickness_ft, Dw_fps) for v in velocities_ft_s]
    ax2.plot(velocities_ft_s, dp_pressures, color='black', linestyle=':', label="Pressure Drop (Pa)")

    # Add pressure drop point for actual face velocity
    dp_actual = mesh_pad_pressure_drop(V_face_actual_ft_s, rho_g_fps, voidage_frac, K_dp, pad_thickness_ft, Dw_fps)
    dp_actual_kpa = dp_actual /1000
    ax2.plot(V_face_actual_ft_s, dp_actual, 's', color='black', markersize=5, label=f"Pressure Drop (ΔP) : {dp_actual_kpa:.2f} kPa")

    ax2.set_ylabel("Pressure Drop (Pa)")

    # Filter legend: only markers
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    lines = []
    labels = []
    for h, l in zip(handles1, labels1):
        if hasattr(h, 'get_marker') and h.get_marker() != 'None':
            lines.append(h)
            labels.append(l)
    for h, l in zip(handles2, labels2):
        if hasattr(h, 'get_marker') and h.get_marker() != 'None':
            lines.append(h)
            labels.append(l)

    ax1.legend(lines, labels, loc='lower right', fontsize=8)

    plt.title(f"Mesh Pad Capture Efficiency & Pressure Drop") #\nQ_gas={Q_gas_ft3_s} ft³/s, Area={A_installed_ft2} ft²")
    plt.tight_layout()
    #return fig

    dp_mesh_pad_actual = mesh_pad_pressure_drop(V_face_actual_ft_s, rho_g_fps, voidage_frac, K_dp, pad_thickness_ft, Dw_fps)
    print(f"Gas flow Q = {Q_gas_ft3_s:.3f} ft³/s")
    print(f"Installed frontal area A = {A_installed_ft2:.3f} ft²")
    print(f"Actual face velocity V_face_actual = {V_face_actual_ft_s:.3f} ft/s ({V_face_actual_ft_s * FT_TO_M:.3f} m/s)")
    print(f"Allowable velocity from K factor = {V_allow_ft_s:.3f} ft/s ({V_allow_ft_s * FT_TO_M:.3f} m/s)")
    print(f"Mesh pad ΔP = {dp_mesh_pad_actual:.2f} Pa ({dp_mesh_pad_actual/1000:.4f} kPa)")

    results["dp_mesh_pad_actual"] = dp_mesh_pad_actual
    results["V_face_actual_ft_s"] = V_face_actual_ft_s
    results["V_allow_ft_s"] = V_allow_ft_s

## This is for finding the carryover. The droplt size will increase by a factor of 2.

    # Define table data
    ratio_percent_points = np.array(
        [100, 103, 105, 107, 110, 112, 114, 116, 118, 120, 122, 123, 124, 125, 127, 128, 131, 133, 136]
    )
    carryover_points = np.array(
        [0, 0.04, 0.1, 2, 3, 5, 8, 11, 16, 22, 31, 38, 46, 56, 68, 76, 91, 97, 100]
    )

    # Calculate ratio in %
    ratio_percent = (V_face_actual_ft_s / V_allow_ft_s) * 100

    # Handle below-minimum and above-maximum cases
    if ratio_percent < ratio_percent_points[0]:
        carryover_percent = carryover_points[0]
    elif ratio_percent > 136:
        carryover_percent = 100
    else:
        carryover_percent = np.interp(ratio_percent, ratio_percent_points, carryover_points)

    print(f"Velocity ratio = {ratio_percent:.2f}%")
    print(f"Estimated carryover = {carryover_percent:.2f}%")
    results["carryover_percent"] = carryover_percent

    return fig


# Function to generate the Flow Regime plot
def create_flow_regime_plot_buffer(plot_data):
    """
    Generates a flow regime map and returns it as a PNG image buffer.
    """
    import matplotlib.pyplot as plt
    import io

    fig, ax = plt.subplots(figsize=(10, 8))
    ax2= ax.twinx()
    
    # Get data from the plot_data dictionary
    line_a = plot_data['line_a']
    line_b = plot_data['line_b']
    line_c_x = plot_data['line_c_x']
    line_c_k = plot_data['line_c_k']
    line_d = plot_data['line_d']
    X = plot_data['X']
    F = plot_data['F']
    K = plot_data['K']

    # Plot the lines on the primary and secondary axes
    ax.plot(line_a[0], line_a[1], 'k-', label='Stratified Transition (Line A-A)')
    ax.plot(line_b[0], line_b[1], 'k--', label='Dispersed Bubble Transition (Line B-B)')
    ax.plot(line_d[0], line_d[1], 'k-.', label='Annular Transition (Line D-D)')
    ax2.plot(line_c_x, line_c_k, 'b--', label='K-based Transition (Line C-C)')

    # Plot the computed point based on the condition
    if K < 6:
        ax2.plot(X, K, 'go', markersize=10, label=f"Computed Point (K={K:.2f})")
    else:
        ax.plot(X, F, 'ro', markersize=10, label=f"Computed Point (X={X:.2f}, F={F:.2f})")

    # Add text labels for regimes
    plt.text(0.01, 13, "STRATIFIED WAVY", fontsize=12, color='BLACK')
    plt.text(0.005, 4000, "ANNULAR - DISPERSED LIQUID", fontsize=12, color='BLACK')
    plt.text(100, 1000, "DISPERSED BUBBLE", fontsize=12, color='BLACK')
    plt.text(100, 50, "INTERMITTENT", fontsize=12, color='BLACK')
    plt.text(0.1, 2, "STRATIFIED SMOOTH", fontsize=12, color='BLACK')

    # Set labels and title
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Martinelli Parameter (X)', fontsize=12)
    ax.set_ylabel('Froude Number-based Parameter (F)', fontsize=12)
    ax2.set_yscale('log')
    ax2.set_ylabel('K Parameter', fontsize=12, color='b')
    ax.set_xlim(0.001, 10000)
    ax.set_ylim(0.001, 10)
    ax2.tick_params(axis='y', colors='blue')
    ax.set_title('Taitel & Dukler Flow Regime Map (1976)', fontsize=14)
    ax.grid(True, which="both", ls="--")
    ax2.set_ylim(1, 10000)

    # Combine legends from both axes
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='upper right')

    # Save to a buffer
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    plt.close(fig)
    return buf


# --- Constants and Look-up Data (derived from the article) ---

# Table 1: Typical Liquid Surface Tension Values (dyne/cm)
# These will be converted to N/m for SI input/display
SURFACE_TENSION_TABLE_DYNE_CM = {
    "Water/gas": 72,
    "Light crude oil/gas": 32,
    "Heavy crude oil/gas": 37,
    "Condensate/gas": 25,
    "Liquid petroleum gas/gas": 10,
    "Natural gas liquids (high C2)/gas": 5,
    "Triethylene glycol": 45,
    "Amine solutions": 52.5, # Average of 45-60 from the range provided
}
# Conversion: 1 dyne/cm = 0.0022 poundal/ft (from the article)
DYNE_CM_TO_POUNDAL_FT = 0.0022
# Conversion: 1 dyne/cm = 0.001 N/m
DYNE_CM_TO_NM = 0.001

# Typical values for Upper-Limit Log Normal Distribution (from the article)
A_DISTRIBUTION = 4.0
DELTA_DISTRIBUTION = 0.72

# Figure 9 Approximation Data for Droplet Size Distribution Shift Factor
# Data extracted from the provided image (Figure 9 table)
SHIFT_FACTOR_DATA = {
    "No inlet device": {
        "rho_v_squared": np.array([0, 100, 200, 300, 400, 500, 600, 650, 675]),
        "shift_factor": np.array([1.00, 0.98, 0.95, 0.90, 0.77, 0.50, 0.20, 0.08, 0.08]) # Using 0.08 as per image
    },
    "Diverter plate": {
        "rho_v_squared": np.array([0,31, 148, 269, 368, 574, 626, 660, 744, 775, 812, 870, 903, 941, 952]),
        "shift_factor": np.array([1,0.99, 0.95, 0.91, 0.87, 0.79, 0.75, 0.71, 0.62, 0.56, 0.52, 0.48, 0.43, 0.35, 0.39])
    },
    "Half-pipe": {
        "rho_v_squared": np.array([0,287, 579, 843, 1064, 1236, 1389, 1519, 1587, 1659, 1743, 1782, 1819, 1892, 1917]),
        "shift_factor": np.array([1,0.96, 0.94, 0.90, 0.87, 0.83, 0.80, 0.71, 0.67, 0.57, 0.43, 0.35, 0.25, 0.07, 0.04])
    },
    "Vane-type": {
        "rho_v_squared": np.array([0,1433, 2297, 3162, 4026, 4891, 5323, 5754, 6229, 6583, 6686, 6775, 6862, 6891, 6894, 6979, 7070]),
        "shift_factor": np.array([1,0.99, 0.97, 0.95, 0.92, 0.89, 0.87, 0.83, 0.78, 0.72, 0.66, 0.60, 0.54, 0.48, 0.43, 0.37, 0.31])
    },
    "Cyclonic": {
        "rho_v_squared": np.array([0, 553, 2294, 2716, 4014, 5312, 5745, 7043, 7908, 8340, 8772, 9205, 9637, 10069, 10501, 10932, 11364, 11794, 12169, 12467, 12716, 12923, 13123, 13323, 13509, 13688, 13891]),
        "shift_factor": np.array([1, 0.99, 0.98, 0.97, 0.97, 0.95, 0.94, 0.93, 0.92, 0.91, 0.90, 0.89, 0.87, 0.86, 0.84, 0.81, 0.77, 0.72, 0.67, 0.61, 0.56, 0.50, 0.45, 0.39, 0.34, 0.28, 0.23])
    }
}


def get_shift_factor(inlet_device, rho_v_squared):
    """
    Calculates the droplet size distribution shift factor using linear interpolation
    based on the inlet device and inlet momentum (rho_g * V_g^2) from Figure 9 data.
    Handles out-of-range values by clamping to the nearest boundary and provides a warning.
    """
    if inlet_device not in SHIFT_FACTOR_DATA:
        st.warning(f"Unknown inlet device: '{inlet_device}'. Defaulting shift factor to 1.0.")
        return 1.0

    data = SHIFT_FACTOR_DATA[inlet_device]
    x_values = data["rho_v_squared"]
    y_values = data["shift_factor"]

    # Check if rho_v_squared is outside the defined range
    if rho_v_squared < x_values.min():
        st.warning(f"Inlet momentum ({rho_v_squared:.2f} lb/ft-sec^2) is below the minimum defined for '{inlet_device}' ({x_values.min():.2f} lb/ft-sec^2). Using minimum shift factor: {y_values.min():.3f}.")
        return y_values.min() # Use the smallest shift factor if below range
    elif rho_v_squared > x_values.max():
        # Corrected: Use the minimum shift factor for values above the maximum defined range
        st.warning(f"Inlet momentum ({rho_v_squared:.2f} lb/ft-sec^2) is above the maximum defined for '{inlet_device}' ({x_values.max():.2f} lb/ft-sec^2). Using minimum shift factor: {y_values.min():.3f}.")
        return y_values.min() # Use the smallest shift factor if above range

    # Perform linear interpolation
    shift_factor = np.interp(rho_v_squared, x_values, y_values)
    return float(shift_factor)


# --- Unit Conversion Factors (for internal FPS calculation) ---
M_TO_FT = 3.28084 # 1 meter = 3.28084 feet
KG_TO_LB = 2.20462 # 1 kg = 2.20462 lb
MPS_TO_FTPS = 3.28084 # 1 m/s = 3.28084 ft/s
PAS_TO_LB_FT_S = 0.67197 # 1 Pa.s (kg/m.s) = 0.67197 lb/ft.s
KG_M3_TO_LB_FT3 = 0.0624279 # 1 kg/m^3 = 0.0624279 lb/ft^3
NM_TO_POUNDAL_FT = 2.2 # 1 N/m = 1000 dyne/cm; 1 dyne/cm = 0.0022 poundal/ft => 1 N/m = 2.2 poundal/ft
IN_TO_FT = 1/12 # 1 inch = 1/12 feet

MICRON_TO_FT = 1e-6 * M_TO_FT
FT_TO_MICRON = 1 / MICRON_TO_FT

def to_fps(value, unit_type):
    """Converts a value from SI to FPS units for internal calculation."""
    if unit_type == "length": # meters to feet
        return value * M_TO_FT
    elif unit_type == "velocity": # m/s to ft/s
        return value * MPS_TO_FTPS
    elif unit_type == "density": # kg/m^3 to lb/ft^3
        return value * KG_M3_TO_LB_FT3
    elif unit_type == "viscosity": # Pa.s to lb/ft.s
        return value * PAS_TO_LB_FT_S
    elif unit_type == "surface_tension": # N/m to poundal/ft
        return value * NM_TO_POUNDAL_FT
    elif unit_type == "pressure": # psig to psi (psig is already a unit of pressure)
        return value # No conversion needed for psig to psi, just use the value directly
    elif unit_type == "diameter_in": # inches to feet
        return value * IN_TO_FT
    return value

def from_fps(value, unit_type):
    """Converts a value from FPS to SI units for display."""
    if unit_type == "length": # feet to meters
        return value / M_TO_FT
    elif unit_type == "velocity": # ft/s to m/s
        return value / MPS_TO_FTPS
    elif unit_type == "density": # lb/ft^3 to kg/m^3
        return value / KG_M3_TO_LB_FT3
    elif unit_type == "viscosity": # lb/ft.s to Pa.s
        return value / PAS_TO_LB_FT_S
    elif unit_type == "momentum": # lb/ft-s^2 to Pa
        return value * 1.48816 # 1 lb/ft-s^2 = 1.48816 Pa
    return value

# --- Functions for E (Entrainment Fraction) Calculation ---
def calculate_e_value(rho_L, rho_G, mu_L, mu_G, sigma_fps, D, u_Gs, Q_liquid_mass_flow_rate):
        # Hardcoded A_2 as per user request
        a2 = 9e-8  #  Refer to Pan and Hanratty (2002)
        #a2 = 8.8e-5  # Large Gas velocity U_g
        
        # Get inputs from the user
        g = 9.81 # float(input("Enter acceleration due to gravity (g): "))
        rho_l = rho_L / KG_M3_TO_LB_FT3 #float(input("Enter liquid density (ρ_L): "))
        rho_g = rho_G / KG_M3_TO_LB_FT3#float(input("Enter gas density (ρ_G): "))
        mu_g = mu_G / PAS_TO_LB_FT_S#float(input("Enter gas viscosity (μ_G): "))
        mu_l = mu_L / PAS_TO_LB_FT_S #float(input("Enter liquid viscosity (μ_L): "))
        d_capital = D / M_TO_FT #float(input("Enter the value for D: "))
        ug = u_Gs #float(input("Enter the value for U_G: "))
        sigma = sigma_fps / NM_TO_POUNDAL_FT #0.0115 #float(input("Enter the value for sigma (σ): "))
        w_l = Q_liquid_mass_flow_rate #float(input("Enter the value for W_L: "))
        
               
        # --- Step 1: Calculate drop diameter 'd' using Eq. 25 ---
        if sigma <= 0 or rho_g <= 0 or ug <= 0:
             raise ValueError("Sigma (σ), Gas Density (ρ_G), and U_G must be positive values.")
        
        # Eq 25 is d = ((D * sigma * 0.0091) / (rho_g * ug**2))^0.5
        d = math.pow((d_capital * sigma * 0.0091) / (rho_g * ug**2), 0.5)
        
        print(f"Calculated drop diameter (d) using Equation 25: {d:.6f}")
        print("-" * 40)
        
        # --- Step 2: Iterative calculation to find a consistent Re_p and m ---
        max_iterations = 500
        tolerance = 1e-6
        cd_guess = 0.4
        
        print("\nStarting iterative process for m...")
        for i in range(max_iterations):
            if cd_guess <= 0:
                print("Error: Cd guess became non-positive. Aborting.")
                return
            
            # Corrected formula for u_t, including buoyancy (rho_l - rho_g)
            u_t = math.sqrt((4 * d * g * (rho_l)) / (3 * cd_guess * rho_g))
            re_p_new = (d * u_t * rho_g) / mu_g
            cd_new = find_cd_from_rep(re_p_new)
            
            if abs(cd_new - cd_guess) < tolerance:
                print(f"Converged after {i+1} iterations.")
                final_re_p = re_p_new
                final_cd = cd_new
                break
            cd_guess = cd_new
        else:
            print("Warning: Calculation did not converge within the maximum number of iterations.")
            final_re_p = re_p_new
            final_cd = cd_new

        print("-" * 40)
        #m_value = calculate_exponent_m(final_re_p)
        if final_re_p < 1.92:
            m_value = 1
        elif final_re_p < 500:
            m_value = 0.6
        else:
            m_value = 0
        print(f"Final Converged Reynolds number (Re_p): {final_re_p:.6f}")
        print(f"Final Converged Drag Coefficient (Cd): {final_cd:.6f}")
        print(f"The calculated exponent 'm' is: {m_value:.6f}")
        print("-" * 40)

        # --- Step 3: Calculate omega (ω), Re_LFC, and E_M ---
        if mu_g == 0 or rho_l == 0:
            raise ValueError("Gas viscosity (μ_G) and Liquid density (ρ_L) cannot be zero.")
            
        omega = (mu_l / mu_g) * math.sqrt(rho_g / rho_l)
        
        if omega <= 0:
            raise ValueError("Omega (ω) must be a positive value for log10 calculation.")
        relfc = 7.3 * (math.log10(omega)) ** 3 + 44.2 * (math.log10(omega)) ** 2 - 263 * (math.log10(omega)) + 439        
        gamma_c = relfc * mu_l / 4
        
        if w_l == 0:
            raise ValueError("W_L cannot be zero.")
        em = 1 - (math.pi * d_capital * gamma_c / w_l)
        
        print(f"Calculated value of omega (ω): {omega:.6f}")
        print(f"Calculated value of Re_LFC: {relfc:.6f}")
        print(f"Calculated value of E_M: {em:.6f}")
        print("-" * 40)

        # --- Step 4: Calculate E using the original equation ---
        term1 = math.sqrt(d_capital * math.pow(ug, 3) * rho_l * rho_g) / sigma
        
        numerator_term2 = math.pow(rho_g, 1 - m_value) * math.pow(mu_g, m_value)
        denominator_term2 = math.pow(d, 1 + m_value) * g * (rho_l - rho_g)
        
        if (2 - m_value) == 0:
            print("Error: The value of 'm' cannot be 2 for this calculation.")
            return

        exponent_term2 = 1 / (2 - m_value)
        
        if denominator_term2 <= 0:
            raise ValueError("Denominator of term 2 must be positive.")
        
        term2 = math.pow(numerator_term2 / denominator_term2, exponent_term2)
        
        rhs = a2 * term1 * term2
        
        if (1 + rhs) == 0:
            raise ValueError("Division by zero. Please check your inputs.")

        e_value = em * (rhs / (1 + rhs))
        return e_value
        

def find_cd_from_rep(re_p):

    if re_p <= 0:
        return float('inf') # Return a very large value to handle non-positive Re_p

    # Apply the Schiller-Naumann correlation
    cd = (24.0 / re_p) * (1.0 + 0.15 * math.pow(re_p, 0.687)) + (0.42 / (1 + (42500 / re_p**1.16)))
    return cd


# Figure 6: C_d vs Re_p data (digitized from plot)
CD_VS_REP_DATA = {
    "Re_p": np.array([
        0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0,
        2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0, 30.0, 40.0, 50.0,
        70.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0, 2000.0,
        5000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0, 500000.0, 1000000.0
    ]),
    "Cd": np.array([
        25000, 12500, 5000, 2500, 1250, 500, 250, 125, 50, 25,
        12.5, 8.3, 6.25, 5.0, 3.5, 2.5, 1.5, 1.0, 0.8, 0.7,
        0.6, 0.5, 0.45, 0.42, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
        0.4, 0.4, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65
    ])
}

def calculate_terminal_velocity(dp_fps, rho_l_fps, rho_g_fps, mu_g_fps, g_fps=32.174):
    """
    Calculates the terminal settling velocity of a droplet using an iterative approach
    based on Eq. 2, Eq. 3, and Figure 6 (Cd vs Re_p).
    All inputs and outputs are in FPS units.
    g_fps: acceleration due to gravity in ft/s^2
    Returns: Vt_new, Cd, Re_p
    """
    if rho_g_fps == 0 or mu_g_fps == 0 or (rho_l_fps - rho_g_fps) <= 0:
        return 0.0, 0.0, 0.0 # Prevent division by zero or non-physical density difference

    # Initial guess for Vt (e.g., using Stokes' Law for small droplets)
    # This initial guess helps the iteration converge faster for typical values.
    # If dp_fps is very small, Stokes' Law is a good start.
    if dp_fps > 0:
        Vt_guess = (g_fps * dp_fps**2 * (rho_l_fps - rho_g_fps)) / (18 * mu_g_fps)
    else:
        Vt_guess = 0.0

    Vt_current = Vt_guess
    tolerance = 1e-6
    max_iterations = 100

    Re_p = 0.0
    Cd = 0.0

    for _ in range(max_iterations):
        if Vt_current <= 0 or dp_fps <= 0: # Handle cases where velocity or diameter is zero/negative
            Re_p = 0.0
        else:
            Re_p = (dp_fps * Vt_current * rho_g_fps) / mu_g_fps

        # Get Cd from Re_p using interpolation from Figure 6 data
        Cd = np.interp(Re_p, CD_VS_REP_DATA["Re_p"], CD_VS_REP_DATA["Cd"])

        # Ensure Cd is not zero or negative
        if Cd <= 0:
            Cd = 0.01 # Small positive value to avoid division by zero, or handle as error

        # Calculate new Vt using Eq. 2
        # Ensure the argument inside sqrt is non-negative
        arg_sqrt = (4 * g_fps * dp_fps * (rho_l_fps - rho_g_fps)) / (3 * Cd * rho_g_fps)
        if arg_sqrt < 0:
            Vt_new = 0.0 # Cannot have imaginary velocity
        else:
            Vt_new = arg_sqrt**0.5

        if abs(Vt_new - Vt_current) < tolerance:
            return Vt_new, Cd, Re_p

        Vt_current = Vt_new

    # If max_iterations reached without convergence, return the last calculated value and warn
    st.warning(f"Terminal velocity calculation did not converge for dp={dp_fps*FT_TO_MICRON:.2f} um after {max_iterations} iterations. Returning last value: {Vt_current:.6f} ft/s.")
    return Vt_current, Cd, Re_p


# Figure 2: F, Actual Velocity/Average (Plug Flow) Velocity vs. L/Di (digitized from plot)
F_FACTOR_DATA = {
    "No inlet device": {
        "L_over_Di": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "F_value": np.array([3.0, 2.5, 2.1, 1.8, 1.5, 1.3, 1.2, 1.1, 1.05, 1.02, 1.0])
    },
    "Diverter plate": {
        "L_over_Di": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "F_value": np.array([2.0, 1.7, 1.5, 1.4, 1.3, 1.2, 1.15, 1.1, 1.05, 1.02, 1.0])
    },
    "Half-pipe": {
        "L_over_Di": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "F_value": np.array([1.4, 1.3, 1.25, 1.2, 1.15, 1.1, 1.08, 1.05, 1.03, 1.02, 1.0])
    },
    "Vane-type": {
        "L_over_Di": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "F_value": np.array([1.3, 1.2, 1.15, 1.1, 1.08, 1.05, 1.03, 1.02, 1.01, 1.0, 1.0])
    },
    "Cyclonic": {
        "L_over_Di": np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        "F_value": np.array([1.2, 1.1, 1.05, 1.03, 1.02, 1.01, 1.0, 1.0, 1.0, 1.0, 1.0])
    }
}

def get_f_factor(inlet_device, L_over_Di, has_perforated_plate):
    """
    Calculates the F factor (Actual Velocity/Average (Plug Flow) Velocity)
    using linear interpolation based on the inlet device and L/Di from Figure 2.
    Applies perforated plate adjustment if selected.
    """
    if inlet_device not in F_FACTOR_DATA:
        st.warning(f"Unknown inlet device for F-factor: '{inlet_device}'. Defaulting F factor to 1.0.")
        return 1.0

    data = F_FACTOR_DATA[inlet_device]
    x_values = data["L_over_Di"]
    y_values = data["F_value"]

    # Clamp L_over_Di to the defined range
    clamped_L_over_Di = max(min(L_over_Di, x_values.max()), x_values.min())

    f_value = np.interp(clamped_L_over_Di, x_values, y_values)

    if has_perforated_plate:
        # Apply perforated plate adjustment: F_effective = F - 0.5 * (F - 1)
        f_value_adjusted = f_value - 0.5 * (f_value - 1)
        # Ensure F_value_adjusted is not less than 1.0 (perfect plug flow)
        return float(max(1.0, f_value_adjusted))

    return float(f_value)


# Table 3: Mesh Pad K Deration Factors as a Function of Pressure
K_DERATION_DATA = {
    "pressure_psig": np.array([0, 100, 200, 400, 600, 800, 1000, 1200]),
    "k_factor_percent": np.array([100, 93, 88, 83, 80, 78, 76, 75])
}

def get_k_deration_factor(pressure_psig):
    """
    Calculates the K deration factor based on pressure using linear interpolation from Table 3.
    """
    if pressure_psig < K_DERATION_DATA["pressure_psig"].min():
        # Clamp to min pressure, use max K factor
        return K_DERATION_DATA["k_factor_percent"].max() / 100.0
    elif pressure_psig > K_DERATION_DATA["pressure_psig"].max():
        # Clamp to max pressure, use min K factor
        return K_DERATION_DATA["k_factor_percent"].min() / 100.0

    k_factor_percent = np.interp(pressure_psig, K_DERATION_DATA["pressure_psig"], K_DERATION_DATA["k_factor_percent"])
    return float(k_factor_percent / 100.0)

# Table 2: Mesh Pad Design and Construction Parameters (FPS units for internal use)
MESH_PAD_PARAMETERS = {
    "Standard mesh pad": {
        "density_lb_ft3": 9,
        "voidage_percent": 98.5,
        "wire_diameter_in": 0.011,
        "specific_surface_area_ft2_ft3": 85,
        "Ks_ft_sec": 0.35,
        "liquid_load_gal_min_ft2": 0.75,
        "thickness_in": 6 # Typical thickness as per Fig 9 example
    },
    "High-capacity mesh pad": {
        "density_lb_ft3": 5,
        "voidage_percent": 99.0,
        "wire_diameter_in": 0.011,
        "specific_surface_area_ft2_ft3": 45,
        "Ks_ft_sec": 0.4,
        "liquid_load_gal_min_ft2": 1.5,
        "thickness_in": 6
    },
    "High-efficiency co-knit mesh pad": {
        "density_lb_ft3": 12,
        "voidage_percent": 96.2,
        "wire_diameter_in": 0.011, # Assuming 0.011 x 0.0008 means effective wire diameter is 0.011
        "specific_surface_area_ft2_ft3": 83, # Using 83, not 1100, as 1100 seems like a typo for a different unit
        "Ks_ft_sec": 0.25,
        "liquid_load_gal_min_ft2": 0.5,
        "thickness_in": 6
    }
}

# Table 4: Vane-Pack Design and Construction Parameters (FPS units for internal use)
VANE_PACK_PARAMETERS = {
    "Simple vane": { # Assuming upflow as default, horizontal is a variant
        "flow_direction": "Upflow", # This is a choice, not a fixed parameter
        "number_of_bends": 5, # Using 5-8, pick 5 as a representative
        "vane_spacing_in": 0.75, # Using 0.5-1, pick 0.75 as a representative
        "bend_angle_degree": 45, # Using 30-60, pick 45 as common
        "Ks_ft_sec_upflow": 0.5, # From table
        "Ks_ft_sec_horizontal": 0.65, # From table
        "liquid_load_gal_min_ft2": 2
    },
    "High-capacity pocketed vane": {
        "flow_direction": "Upflow", # This is a choice, not a fixed parameter
        "number_of_bends": 5,
        "vane_spacing_in": 0.75,
        "bend_angle_degree": 45,
        "Ks_ft_sec_upflow": 0.82, # Using 0.82-1.15, pick 0.82
        "Ks_ft_sec_horizontal": 0.82, # Using 0.82-1.15, pick 0.82
        "liquid_load_gal_min_ft2": 5
    }
}

# Table 5: Typical Demisting Axial-Flow Cyclone Design and Construction Parameters (FPS units for internal use)
CYCLONE_PARAMETERS = {
    "2.0 in. cyclones": { # Only one type given, so this is the default
        "cyclone_inside_diameter_in": 2.0,
        "cyclone_length_in": 10,
        "inlet_swirl_angle_degree": 45,
        "cyclone_to_cyclone_spacing_diameters": 1.75,
        "Ks_ft_sec_bundle_face_area": 0.8, # Using ~0.8-1, pick 0.8
        "liquid_load_gal_min_ft2_bundle_face_area": 10
    }
}


# Figure 8: Single-wire droplet capture efficiency (Ew) vs. Stokes' number (Stk)
# Curve fit given by Eq. 13: Ew = (-0.105 + 0.995 * Stk^0.0493) / (0.6261 + Stk^1.00493)
def calculate_single_wire_efficiency(Stk):
    """Calculates single-wire impaction efficiency using Eq. 13."""
    if Stk <= 0: # Handle Stk=0 or negative to avoid math domain errors
        return 0.0
    numerator = -0.105 + 0.995 * (Stk**1.00493)
    denominator = 0.6261 + (Stk**1.00493)
    if denominator == 0: # Avoid division by zero
        return 0.0
    Ew = numerator / denominator
    return max(0.0, min(1.0, Ew)) # Ensure efficiency is between 0 and 1

# --- Mist Extractor Efficiency Functions ---

def mesh_pad_efficiency_func(dp_fps, V_g_eff_sep_fps, rho_l_fps, rho_g_fps, mu_g_fps, mesh_pad_type_params_fps):
    """
    Calculates the droplet removal efficiency for a mesh pad using Equations 12, 13, and 14.
    All inputs in FPS units.
    """
    if V_g_eff_sep_fps <= 0 or mu_g_fps <= 0 or dp_fps <= 0:
        return 0.0, 0.0, 0.0 # No impaction if no gas flow or zero droplet/gas viscosity

    Dw_fps = mesh_pad_type_params_fps["wire_diameter_in"] * IN_TO_FT
    pad_thickness_fps = mesh_pad_type_params_fps["thickness_in"] * IN_TO_FT
    specific_surface_area_fps = mesh_pad_type_params_fps["specific_surface_area_ft2_ft3"]

    # Eq. 12: Stokes' number
    # Note: Article states some literature uses 9 in denominator instead of 18. Using 18 as per Eq. 12.
    if Dw_fps == 0: return 0.0, 0.0, 0.0  # Avoid division by zero if wire diameter is zero
    Stk = ((rho_l_fps - rho_g_fps) * (dp_fps**2) * V_g_eff_sep_fps) / (18 * mu_g_fps * Dw_fps)

    # Eq. 13: Single-wire capture efficiency
    Ew = calculate_single_wire_efficiency(Stk)

    # Eq. 14: Mesh-pad removal efficiency
    exponent = -0.238 * specific_surface_area_fps * pad_thickness_fps * Ew
    E_pad = 1 - np.exp(exponent)
    E_pad = max(0.0, min(1.0, E_pad))
    return E_pad, Stk, Ew
    #return max(0.0, min(1.0, E_pad)) # Ensure efficiency is between 0 and 1

def vane_type_efficiency_func(dp_fps, V_g_eff_sep_fps, rho_l_fps, rho_g_fps, mu_g_fps, vane_type_params_fps):
    """
    Calculates the droplet separation efficiency for a vane-type mist extractor using Eq. 15.
    All inputs in FPS units.
    """
    if V_g_eff_sep_fps <= 0 or mu_g_fps <= 0 or dp_fps <= 0:
        return 0.0

    num_bends = vane_type_params_fps["number_of_bends"]
    vane_spacing_fps = vane_type_params_fps["vane_spacing_in"] * IN_TO_FT
    bend_angle_rad = np.deg2rad(vane_type_params_fps["bend_angle_degree"])
    bend_angle_deg = (vane_type_params_fps["bend_angle_degree"])

    numerator = num_bends * (dp_fps**2) * (rho_l_fps - rho_g_fps) * V_g_eff_sep_fps * bend_angle_deg
    # Eq. 15: Evane = 1 - exp[ - (n * dp^3 * (rho_l - rho_g) * Vg_eff_sep) / (515.7 * mu_g * b * cos^2(theta)) ]
    #numerator = num_bends * (dp_fps**2) * (rho_l_fps - rho_g_fps) * V_g_eff_sep_fps * np.cos(bend_angle_rad) theta should be retained in degrees , no cos function
    denominator = 515.7 * mu_g_fps * vane_spacing_fps * (np.cos(bend_angle_rad)**2)

    if denominator == 0:
        return 0.0

    exponent = - (numerator / denominator)
    E_vane = 1 - np.exp(exponent)
    E_vane = max(0.0, min(1.0, E_vane))

    return E_vane

    #return max(0.0, min(1.0, E_vane)) # Ensure efficiency is between 0 and 1

def demisting_cyclone_efficiency_func(dp_fps, V_g_eff_sep_fps, rho_l_fps, rho_g_fps, mu_g_fps, cyclone_type_params_fps):
    """
    Calculates the droplet removal efficiency for an individual axial-flow cyclone tube
    using Eq. 16 and the associated Stokes' number definition.
    All inputs in FPS units.
    """
    if V_g_eff_sep_fps <= 0 or mu_g_fps <= 0 or dp_fps <= 0:
        return 0.0, 0.0

    Dcycl_fps = cyclone_type_params_fps["cyclone_inside_diameter_in"] * IN_TO_FT
    Lcycl_fps = cyclone_type_params_fps["cyclone_length_in"] * IN_TO_FT
    in_swirl_angle_rad = np.deg2rad(cyclone_type_params_fps["inlet_swirl_angle_degree"])

    # Eq. 16: Stk_cycl = ( (rho_l - rho_g) * dp^2 * Vg_cycl ) / (18 * mu_g * Dcycl)
    # Vg_cycl is superficial gas velocity through a single cyclone tube.
    # Assuming V_g_eff_sep_fps is the superficial velocity through the cyclone bundle face area,
    # and this can be used as Vg_cycl for a single cyclone for efficiency calculation.
    Vg_cycl = V_g_eff_sep_fps # Approximation for simplicity as bundle area is not easily translated to single tube area without more info.

    if Dcycl_fps == 0: return 0.0, 0.0 # Avoid division by zero
    Stk_cycl = ((rho_l_fps - rho_g_fps) * (dp_fps**2) * Vg_cycl) / (18 * mu_g_fps * Dcycl_fps)

    # Eq. 16: E_cycl = 1 - exp[ -8 * Stk_cycl * (Lcycl / (Dcycl * tan(alpha))) ]
    # Ensure tan(alpha) is not zero or near zero for 90 degree swirl angle etc.
    if np.tan(in_swirl_angle_rad) == 0:
        return 0.0, 0.0 # No swirl, no separation

    exponent = -8 * Stk_cycl * (Lcycl_fps / (Dcycl_fps * np.tan(in_swirl_angle_rad)**2))
    E_cycl = 1 - np.exp(exponent)
    E_cycl = max(0.0, min(1.0, E_cycl))

    return E_cycl, Stk_cycl

# --- PDF Report Generation Function ---
class PDF(FPDF):
    def header(self):
        # Draw a border rectangle (x, y, width, height)
        self.set_line_width(0.5)
        self.rect(5.0, 5.0, self.w - 10.0, self.h - 10.0)  # 5 units margin on all sides
        if self.page_no() > 1:
            self.set_font('Arial', 'B', 12)
            self.cell(0, 10, 'Oil and Gas Separation: Liquid in Gas Carry Over Report', 0, 1, 'C')
            self.ln(5)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', 0, 0, 'C')

    def chapter_title(self, title):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(2)

    def chapter_body(self, body):
        self.set_font('Arial', '', 10)
        indent_amount = 5 # Adjust this value as needed
        self.cell(indent_amount, 5, "", 0, 0)

        self.multi_cell(self.w - self.l_margin - self.r_margin - indent_amount, 4, body)
        self.ln()

    def chapter_body1(self, body):
        self.set_font('Arial', 'B', 10)
        indent_amount = 5 # Adjust this value as needed
        self.cell(indent_amount, 5, "", 0, 0)

        self.multi_cell(self.w - self.l_margin - self.r_margin - indent_amount, 4, body)
        self.ln()

    def add_table(self, headers, data, col_widths, title=None):
        # Calculate the total width of the table
        total_table_width = sum(col_widths)
        # Calculate the left margin needed to center the table on the page
        left_margin = (self.w - total_table_width) / 2
        
        if title:
            self.set_font('Arial', 'B', 10)
            # Set the x-coordinate for the centered title
            self.set_x(left_margin)
            self.cell(total_table_width, 7, title, 0, 1, 'C') # Centered within the table width
            self.ln(5)

        # Set font for table headers
        self.set_font('Arial', 'B', 9)
        self.set_x(left_margin) # Position the header row
        for i, header in enumerate(headers):
            self.cell(col_widths[i], 7, header, 1, 0, 'C')
        self.ln()

        # Set font for table data
        self.set_font('Arial', '', 8)
        for row in data:
            self.set_x(left_margin) # Position each data row
            for i, item in enumerate(row):
                self.cell(col_widths[i], 6, str(item), 1, 0, 'C')
            self.ln()
        self.ln(5)

#@st.cache_data
#def generate_pdf_report(inputs, results, plot_image_buffer_original, plot_image_buffer_adjusted, plot_data_original, plot_data_adjusted, plot_data_after_gravity, plot_data_after_mist_extractor):
def generate_pdf_report(
    inputs, results,
    plot_image_buffer_original, plot_image_buffer_adjusted,
    plot_image_buffer_after_gravity, plot_image_buffer_after_me,
    plot_data_original, plot_data_adjusted,
    plot_data_after_gravity, plot_data_after_mist_extractor,
    client_name, contractor_name, equipment_tag, prepared_by, rev_no, project_name
):
    pdf = PDF()
    pdf.alias_nb_pages()
    pdf.add_page()

    # Title Page
    pdf.ln(20)
    pdf.set_font('Arial', 'B', 18)
    pdf.cell(0, 20, 'Oil and Gas Separation: Liquid in Gas Carry Over Report', 0, 1, 'C')
    pdf.ln(50)

    try:
        pdf.image("Sep.png", x=pdf.w/2 - 25, y=80, w=50)  # Centered, adjust w as needed
    except Exception as e:
        print("Could not add Sep.png to PDF:", e)

    pdf.ln(20)
    pdf.set_font('Arial', 'B', 18)
    pdf.cell(0, 10, f"EQUIPMENT TAG", 0, 1, 'C')
    pdf.cell(0, 10, f"{equipment_tag}", 0, 1, 'C')
    pdf.ln(50)
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, f"Project Name: {project_name}", 0, 1, 'L')
    pdf.cell(0, 10, f"Client Name: {client_name}", 0, 1, 'L')
    pdf.cell(0, 10, f"Contractor Name: {contractor_name}", 0, 1, 'L')
    
    pdf.cell(0, 10, f"Prepared By: {prepared_by}", 0, 1, 'L')
    pdf.cell(0, 10, f"Revision No.: {rev_no}", 0, 1, 'L')
    pdf.cell(0, 10, f"Date: {st.session_state.report_date}", 0, 1, 'L')
    pdf.ln(10)

    # --- Summary of Results ---
    pdf.add_page()
    pdf.chapter_title('Results')                   

    entrainment_after_inlet_device_kg_hr = st.session_state.plot_data_adjusted.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
    inlet_device_efficiency_percent = st.session_state.calculation_results.get('inlet_device_efficiency_percent', 0.0)
    entrainment_after_me_kg_hr = st.session_state.plot_data_after_mist_extractor.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
    entrainment_after_gravity_kg_hr = st.session_state.plot_data_after_gravity.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
    mist_extractor_efficiency = st.session_state.calculation_results.get('mist_extractor_efficiency', 0.0)
    gasgravity_extractor_efficiency = st.session_state.calculation_results.get('gasgravity_extractor_efficiency', 0.0)
    final_carryover_gal_mmscf = st.session_state.results.get('final_carryover_gal_mmscf', 0.0)
    flow_regime_result = st.session_state.results.get('flow_regime_td_result', 'N/A')
    
    # Define the data for the summary table
    summary_headers = ["Parameter", "Value"]
    summary_data = [
        ["Mean particle diameter dv50 at the inlet (estimated)", f"{results['dv50_original_fps'] * FT_TO_MICRON:.2f} um"],
        ["Inlet flow regime", flow_regime_result], 
        ["Entrainment at inlet pipe", f"{results['E_fraction']*100:.2f} %"],
        ["Entrainment after Inlet Device", f"{entrainment_after_inlet_device_kg_hr:.2f} kg/hr"],
        ["Inlet device efficiency", f"{inlet_device_efficiency_percent:.2f} %"],
        ["Entrainment after Gravity Settling Section", f"{entrainment_after_gravity_kg_hr:.2f} kg/hr"],
        ["Gravity Settling Section efficiency", f"{gasgravity_extractor_efficiency:.2f} %"],
        ["Entrainment after Mist Extractor Section", f"{entrainment_after_me_kg_hr:.2f} kg/hr"],
        ["Mist extractor efficiency", f"{mist_extractor_efficiency:.2f} %"],
        ["Separator Performance - Liquid Carryover ", f"{final_carryover_gal_mmscf:.4f} US Gal / MMSCF"]
    ]

    # Define column widths for the table
    summary_col_widths = [100, 70] # Adjust as needed

    # Add the table to the PDF
    pdf.add_table(summary_headers, summary_data, summary_col_widths, title="Separator Performance Summary")                      

    # --- Input Parameters ---
    pdf.add_page()
    pdf.chapter_title('1. Input Parameters (SI Units)')
    # Replace problematic characters for PDF output
    pdf.chapter_body(f"Pipe Inside Diameter (D): {inputs['D_input']:.4f} m")
    pdf.chapter_body(f"Liquid Density (rho_l): {inputs['rho_l_input']:.2f} kg/m³")
    pdf.chapter_body(f"Liquid Viscosity (mu_l): {inputs['mu_l_input']:.8f} Pa.s")
   #pdf.chapter_body(f"Gas Velocity (Vg): {inputs['V_g_input']:.2f} m/s")
    pdf.chapter_body(f"Gas Density (rho_g): {inputs['rho_g_input']:.5f} kg/m³")
    pdf.chapter_body(f"Gas Viscosity (mu_g): {inputs['mu_g_input']:.9f} Pa.s")

    sigma_display_val = inputs['sigma_custom'] # Directly use sigma_custom for display
    pdf.chapter_body(f"Liquid Surface Tension (sigma): {sigma_display_val:.3f} N/m")
    pdf.chapter_body(f"Selected Inlet Device: {inputs['inlet_device']}")
    pdf.chapter_body(f"Total Liquid Mass Flow Rate: {inputs['Q_liquid_mass_flow_rate_input']:.2f} kg/s") # New input
    pdf.chapter_body(f"Number of Points for Distribution: {inputs['num_points_distribution']}") # New input
    pdf.ln(5)

    pdf.chapter_body1(f"-----Separator Type: {inputs['separator_type']}-----")
    if inputs['separator_type'] == "Horizontal":
        pdf.chapter_body(f"Gas Space Height (hg): {inputs['h_g_input']:.3f} m")
        pdf.chapter_body(f"Effective Separation Length (Le): {inputs['L_e_input']:.3f} m")
    else: # Vertical
        pdf.chapter_body(f"Separator Diameter: {inputs['D_separator_input']:.3f} m")

    pdf.chapter_body(f"Length from Inlet Device to Mist Extractor (L_to_ME): {inputs['L_to_ME_input']:.3f} m")
    pdf.chapter_body(f"Perforated Plate Used: {'Yes' if inputs['perforated_plate_option'] else 'No'}")
    pdf.chapter_body(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} psig")
    pdf.ln(5)

    pdf.chapter_body1(f"-----Mist Extractor Type: {inputs['mist_extractor_type']}-----")
    if inputs['mist_extractor_type'] == "Mesh Pad":
        pdf.chapter_body(f"Mesh Pad Type: {inputs['mesh_pad_type']}")
        pdf.chapter_body(f"Mesh Pad Thickness: {inputs['mesh_pad_thickness_in']:.2f} in")
    elif inputs['mist_extractor_type'] == "Vane-Type":
        pdf.chapter_body(f"Vane Type: {inputs['vane_type']}")
        pdf.chapter_body(f"Flow Direction: {inputs['vane_flow_direction']}")
        pdf.chapter_body(f"Number of Bends: {inputs['vane_num_bends']}")
        pdf.chapter_body(f"Vane Spacing: {inputs['vane_spacing_in']:.2f} in")
        pdf.chapter_body(f"Bend Angle: {inputs['vane_bend_angle_deg']:.1f} deg")
    elif inputs['mist_extractor_type'] == "Cyclonic":
        pdf.chapter_body(f"Cyclone Type: {inputs['cyclone_type']}")
        pdf.chapter_body(f"Cyclone Diameter: {inputs['cyclone_diameter_in']:.2f} in")
        pdf.chapter_body(f"Cyclone Length: {inputs['cyclone_length_in']:.2f} in")
        pdf.chapter_body(f"Inlet Swirl Angle: {inputs['cyclone_swirl_angle_deg']:.1f} deg")
    pdf.ln(5)

# --- Start of Flow Regime Section ---
    pdf.add_page()
    pdf.chapter_title('2. Flow Regime Analysis')

    rho_L = st.session_state.inputs['rho_l_input']
    mu_L = st.session_state.inputs['mu_l_input']
    rho_G = st.session_state.inputs['rho_g_input']
    mu_G = st.session_state.inputs['mu_g_input']
    Q_liquid_mass_flow_rate = st.session_state.inputs['Q_liquid_mass_flow_rate_input']
    Q_gas_mass_flow_rate = st.session_state.inputs['Q_gas_mass_flow_rate_input']
    D = st.session_state.inputs['D_input'] 
    
    # Calculate superficial velocities from mass flow rates
    u_Ls = (Q_liquid_mass_flow_rate / rho_L) / (0.25 * math.pi * D**2) if rho_L > 0 and D > 0 else 0
    u_Gs = (Q_gas_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0

       
    # Perform calculations by calling functions from the flow_regime module
    nu_L = mu_L / rho_L
    X = fr.compute_X(rho_L, rho_G, mu_L, mu_G, u_Ls, u_Gs, D)
    alpha = 0 # Assuming alpha is a constant, can be adjusted based on user input or conditions later versions
    F = fr.compute_F(rho_L, rho_G, u_Gs, D, alpha=alpha)
    K = fr.compute_K(F, D, u_Ls, nu_L)
            
    x_d = [1.62, 2.23086, 3.565971, 5.566819, 9.155307, 15.47836, 23.27883, 37.2084, 61.40991, 100.0441, 172.6473, 275.9255, 413.666, 740.6663, 1226.437, 1896.173, 2938.08]
    y_d = [1.25, 1.160695, 1.085948, 1.032387, 0.97697, 0.85925, 0.798867, 0.727388, 0.679639, 0.579151, 0.497267, 0.430495, 0.377634, 0.307546, 0.253466, 0.221436, 0.178456]
    line_d = (x_d, y_d)

    x_b = [1.58, 1.58]
    y_b = [0.15, 9.71]
    line_b = (x_b, y_b)

    x_a = [0.0001, 0.0042, 0.00675, 0.01079, 0.01653, 0.02878, 0.04599, 0.0735, 0.1175, 0.1877, 0.2999, 0.4792, 0.7656, 1.2228, 1.8322, 2.6304, 3.8576, 5.5376, 8.4734, 11.907, 16.034, 21.592, 29.071, 36.726, 46.386]
    y_a = [2.5, 1.642, 1.555, 1.432, 1.327, 1.209, 1.075, 0.932, 0.805, 0.663, 0.529, 0.397, 0.286, 0.180, 0.121, 0.0789, 0.0503, 0.0311, 0.0181, 0.0114, 0.0072, 0.00458, 0.00272, 0.00181, 0.00110]
    line_a = (x_a, y_a)

    x_c = [0.0189, 0.0306, 0.0498, 0.0823, 0.1350, 0.2383, 0.4006, 0.6174, 1.0734, 1.6345, 2.8599, 4.7267, 7.6305, 12.698, 22.888, 41.846]
    k_c = [2.0759, 2.6344, 3.3149, 4.1470, 4.9008, 5.6003, 5.9443, 6.0132, 5.8264, 5.3906, 4.6319, 3.8222, 3.1353, 2.4390, 1.6844, 1.1973]
   
    flow_regime = fr.get_flow_regime(X, F, line_a, line_b, line_d)

    pdf.set_font('Arial', '', 10)
    pdf.chapter_body(f"Martinelli Parameter (X) = {X:.2f}")
    pdf.chapter_body(f"Froude Number-based Parameter (F) = {F:.2f}")
    pdf.chapter_body(f"K Parameter (K) = {K:.2f}")
    pdf.ln(5)
    pdf.set_font('Arial', 'B', 10)
    pdf.chapter_body(f"The calculated flow regime is: {flow_regime}")
    pdf.ln(5)

    # Prepare data for plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    ax2 = ax.twinx()

    plot_data = {
        'X': X,
        'K': K,
        'F': F,
        'line_a': line_a,
        'line_b': line_b,
        'line_c_x': x_c,
        'line_c_k': k_c,
        'line_d': line_d
    }

    ax.plot(line_a[0], line_a[1], 'k-', label='Stratified Transition (Line A-A)')
    ax.plot(line_b[0], line_b[1], 'k--', label='Dispersed Bubble Transition (Line B-B)')
    ax2.plot(x_c, k_c, 'b--', label='K-based Transition (Line C-C)')
    ax.plot(line_d[0], line_d[1], 'k-.', label='Annular Transition (Line D-D)')

    # Then, use an if/else statement to perform the plotting
    if plot_data['K'] < 6:
        ax2.plot(plot_data['X'], plot_data['K'], 'go', markersize=10, label=f"Computed Point (K={plot_data['K']:.2f})")
    else:
        ax.plot(plot_data['X'], plot_data['F'], 'ro', markersize=10, label=f"Computed Point (X={plot_data['X']:.2f}, F={plot_data['F']:.2f})")
    
    # Create the plot buffer
    flow_regime_buf = create_flow_regime_plot_buffer(plot_data)

    # Add the plot to the PDF
    temp_image_path = "flow_regime_temp.png"
    with open(temp_image_path, "wb") as f:
        f.write(flow_regime_buf.getvalue())
        
    # Add the plot to the PDF using the file path
    image_width_mm = 150
    image_x = (pdf.w - image_width_mm) / 2
    pdf.image(temp_image_path, x=image_x, w=image_width_mm)
    pdf.ln(10)

    # Clean up: delete the temporary file
    import os
    os.remove(temp_image_path)
    # --- End of Flow Regime Section ---

    # --- Calculation Steps ---
    pdf.add_page()
    pdf.chapter_title('3. Calculation Results')

    # Define unit labels for SI system for report (using ASCII-safe versions)
    len_unit_pdf = "m"
    dens_unit_pdf = "kg/m³"
    vel_unit_pdf = "m/s"
    visc_unit_pdf = "Pa.s"
    momentum_unit_pdf = "Pa"
    micron_unit_label_pdf = "um"
    mass_flow_unit_pdf = "kg/s"
    vol_flow_unit_pdf = "m³/s" # New unit for PDF
    pressure_unit_pdf = "psig"

    D_pipe_si = inputs['D_input']
    Q_gas_mass_flow_rate_input_si = inputs['Q_gas_mass_flow_rate_input']
    rho_g_input_si = inputs['rho_g_input']
    Vol_flow_si = Q_gas_mass_flow_rate_input_si / rho_g_input_si if rho_g_input_si > 0 else 0.0
    Ug_si = Vol_flow_si / (0.785 * (D_pipe_si)**2) if D_pipe_si > 0 else 0.0
    Q_liquid_mass_flow_rate_input_si = inputs['Q_liquid_mass_flow_rate_input']
    rho_l_input_si = inputs['rho_l_input']
    V_g_input_fps = Ug_si * M_TO_FT


    def create_latex_image(latex_string, filename):
        """Generates a PNG image from a LaTeX string."""
        fig = plt.figure()
        fig.text(0.5, 0.5, f"${latex_string}$", fontsize=30, ha='center', va='center')
        plt.axis('off')
        plt.savefig(filename, bbox_inches='tight', transparent=True)
        plt.close(fig)
    # Example usage for your equation
    dv50_eq = "d_{v50} = 0.01 \\cdot \\frac{\\sigma}{\\rho_g V_g^2} \\cdot Re_g^{2/3} \\cdot \\left(\\frac{\\rho_g}{\\rho_l}\\right)^{-1/3} \\cdot \\left(\\frac{\\mu_g}{\\mu_l}\\right)^{2/3}"
    create_latex_image(dv50_eq, "dv50_eq.png")
    Re_g_eq = "Re_g = \\frac{D \\cdot V_g \\cdot \\rho_g}{\\mu_g}"
    create_latex_image(Re_g_eq, "Re_g_eq.png")
    ent_eq1 = r"""\frac{\left( \frac{E}{E_M} \right)}{1 - \left( \frac{E}{E_M} \right)} = A_2 \left( \frac{D U_G^3 \rho_L^{0.5} \rho_G^{0.5}}{\sigma} \right)\left( \frac{\rho_G^{1-m} \mu_G^{m}}{d^{1+m} g \rho_L} \right)^{\frac{1}{2-m}}"""
    create_latex_image(ent_eq1, "ent_eq1.png")
    # Second Equation
    ent_eq2 = "\\left( \\frac{\\rho_G U_G^2 d_{32}}{\\sigma} \\right)\\left( \\frac{d_{32}}{D} \\right) = 0.0091"
    create_latex_image(ent_eq2, "ent_eq2.png")
    

    pdf.set_font('Arial', 'B', 10)
    pdf.chapter_body1("Inputs Used for Calculation (Converted to FPS for internal calculation):")
    pdf.set_font('Arial', '', 10)
    pdf.chapter_body(f"Pipe Inside Diameter (D): {to_fps(inputs['D_input'], 'length'):.2f} ft")
    pdf.chapter_body(f"Liquid Density (rho_l): {to_fps(inputs['rho_l_input'], 'density'):.2f} lb/ft^3")
    pdf.chapter_body(f"Liquid Viscosity (mu_l): {to_fps(inputs['mu_l_input'], 'viscosity'):.7f} lb/ft-sec")
    pdf.chapter_body(f"Superficial Gas Velocity (Vg): {V_g_input_fps:.2f} ft/sec")
    pdf.chapter_body(f"Gas Density (rho_g): {to_fps(inputs['rho_g_input'], 'density'):.4f} lb/ft^3")
    pdf.chapter_body(f"Gas Viscosity (mu_g): {to_fps(inputs['mu_g_input'], 'viscosity'):.8f} lb/ft-sec")
    pdf.chapter_body(f"Liquid Surface Tension (sigma): {inputs['sigma_fps']:.4f} poundal/ft")
    pdf.chapter_body(f"Liquid Mass Flow Rate: {inputs['Q_liquid_mass_flow_rate_input']:.2f} {mass_flow_unit_pdf}") # New input
    pdf.chapter_body(f"Gas Mass Flow Rate: {inputs['Q_gas_mass_flow_rate_input']:.2f} {mass_flow_unit_pdf}") # New input
    pdf.chapter_body(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} {pressure_unit_pdf}")
    pdf.ln(5)

    # Step 1
    pdf.set_font('Arial', 'B', 10)
    pdf.chapter_body1("Step 1: Calculate Superficial Gas Reynolds Number (Re_g)")
    pdf.set_font('Arial', '', 10)
    image_y = pdf.get_y() - 15
    image_height = 50
    pdf.image("Re_g_eq.png", x=pdf.w/2 - 105, y=image_y, w=65, h=image_height)
    pdf.set_y(image_y + image_height-15)

    #pdf.chapter_body(f"Equation: Re_g = (D * V_g * rho_g) / mu_g")
    pdf.chapter_body(f"Calculation (FPS): Re_g = ({to_fps(inputs['D_input'], 'length'):.2f} ft * {V_g_input_fps:.2f} ft/sec * {to_fps(inputs['rho_g_input'], 'density'):.4f} lb/ft^3) / {to_fps(inputs['mu_g_input'], 'viscosity'):.8f} lb/ft-sec = {results['Re_g']:.2f}")
    pdf.chapter_body(f"Result: Superficial Gas Reynolds Number (Re_g) = {results['Re_g']:.2f} (dimensionless)")
    pdf.ln(5)

    # Step 2
    pdf.set_font('Arial', 'B', 10)
    pdf.chapter_body1("Step 2: Calculate Initial Volume Median Diameter (d_v50) (Kataoka et al., 1983)")
    pdf.set_font('Arial', '', 10)
    #pdf.chapter_body(f"Equation: d_v50 = 0.01 * (sigma / (rho_g V_g^2)) * Re_g^(2/3) * (rho_g / rho_l)^(-1/3) * (mu_g / mu_l)^(2/3)")
    # Place the image 5mm below the current cursor position
    image_y = pdf.get_y() - 15
    # The image height will be determined by its original aspect ratio,
    # or you can set a specific height (e.g., h=20).
    image_height = 50
    # Use the current cursor position for the y-coordinate
    pdf.image("dv50_eq.png", x=pdf.w/2 - 90, y=image_y, w=100, h=image_height)
    # Manually advance the cursor to the bottom of the image
    pdf.set_y(image_y + image_height-15)
    #pdf.image("equation.png", x=pdf.w/2 - 25, y=60, w=50)  # Centered, adjust w as needed
    
    pdf.chapter_body(f"Calculation (FPS): d_v50 = 0.01 * ({inputs['sigma_fps']:.4f} / ({to_fps(inputs['rho_g_input'], 'density'):.4f} * {V_g_input_fps:.2f}^2)) * ({results['Re_g']:.2f})^(2/3) * ({to_fps(inputs['rho_g_input'], 'density'):.4f} / {to_fps(inputs['rho_l_input'], 'density'):.2f})^(-1/3) * ({to_fps(inputs['mu_g_input'], 'viscosity'):.8f} / {to_fps(inputs['mu_l_input'], 'viscosity'):.7f})^(2/3) = {results['dv50_original_fps']:.6f} ft")
    pdf.chapter_body(f"Result: Initial Volume Median Diameter (d_v50) = {results['dv50_original_fps'] * FT_TO_MICRON:.2f} {micron_unit_label_pdf} ({from_fps(results['dv50_original_fps'], 'length'):.6f} {len_unit_pdf})")
    pdf.ln(5)

    rho_L = st.session_state.inputs['rho_l_input']
    mu_L = st.session_state.inputs['mu_l_input']
    rho_G = st.session_state.inputs['rho_g_input']
    mu_G = st.session_state.inputs['mu_g_input']
    Q_liquid_mass_flow_rate = st.session_state.inputs['Q_liquid_mass_flow_rate_input']
    Q_gas_mass_flow_rate = st.session_state.inputs['Q_gas_mass_flow_rate_input']
    D = st.session_state.inputs['D_input']  
    u_Gs = (Q_gas_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0
    u_ls = (Q_liquid_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0

    Q_gas_mass_lb_s = Q_gas_mass_flow_rate * 2.20462
    Q_liquid_mass_lb_s = Q_liquid_mass_flow_rate * 2.20462
    Q_gas_vol_ft3_s = Q_gas_mass_lb_s / (rho_G * 0.062428)  # Convert gas density from kg/m³ to lb/ft³
    Q_liquid_vol_ft3_s = Q_liquid_mass_lb_s / (rho_L * 0.062428)  # Convert liquid density from kg/m³ to lb/ft³
        
    total_vol_flow = Q_gas_vol_ft3_s + Q_liquid_vol_ft3_s
    # Volume fractions
    alpha_g = Q_gas_vol_ft3_s / total_vol_flow
    alpha_l = 1 - alpha_g

    rho_mix_fps = alpha_g * (rho_G *0.062428) + alpha_l * (rho_L *0.062428)  # Convert to lb/ft³
    # Convert velocities to ft/s    
    area_ft2 = math.pi * ((D * M_TO_FT) ** 2) / 4.0
    V_mix_fps = total_vol_flow / area_ft2 if area_ft2 > 0 else 0.0

    rho_v_squared_fps = rho_mix_fps * V_mix_fps**2


    # Step 3
    pdf.set_font('Arial', 'B', 10)
    pdf.chapter_body1("Step 3: Calculate Inlet Momentum (rho_m V_m^2)")
    pdf.set_font('Arial', '', 10)
    pdf.chapter_body(f"Equation: rho_m * V_m^2")
    #pdf.chapter_body(f"Calculation (FPS): rho_g V_g^2 = {to_fps(inputs['rho_g_input'], 'density'):.4f} lb/ft^3 * ({to_fps(inputs['V_g_input'], 'velocity'):.2f} ft/sec)^2 = {results['rho_v_squared_fps']:.2f} lb/ft-sec^2")
    pdf.chapter_body(f"Calculation (FPS): rho_m V_m^2 = {rho_mix_fps:.4f} lb/ft³ * ({V_mix_fps:.2f}2^2) = {results['rho_v_squared_fps']:.2f} lb/ft-sec^2")
    pdf.chapter_body(f"Result: Inlet Momentum (rho_m V_m^2) = {from_fps(results['rho_v_squared_fps'], 'momentum'):.2f} {momentum_unit_pdf}")
    pdf.ln(5)

    # Step 4
    pdf.add_page()
    pdf.chapter_body1("Step 4: Liquid Separation Efficiency / Droplet Size Distribution Shift Factor")
    
    pdf.chapter_body(f"Selected Inlet Device: {inputs['inlet_device']}")
    pdf.chapter_body(f"Estimated Shift Factor : {results['shift_factor']:.3f}")
    pdf.chapter_body(f"Equation: d_v50,adjusted = d_v50,original * Shift Factor")
    pdf.chapter_body(f"Calculation (FPS): d_v50,adjusted = {results['dv50_original_fps']:.6f} ft * {results['shift_factor']:.3f} = {results['dv50_adjusted_fps']:.6f} ft")
    pdf.chapter_body(f"Result: Adjusted Volume Median Diameter (d_v50) = {results['dv50_adjusted_fps'] * FT_TO_MICRON:.2f} {micron_unit_label_pdf} ({from_fps(results['dv50_adjusted_fps'], 'length'):.6f} {len_unit_pdf})")
    pdf.ln(5)

    # Step 5
    
    pdf.chapter_body1("Step 5: Calculate Parameters for Upper-Limit Log Normal Distribution")
    
    pdf.chapter_body(f"Using typical values, a = {A_DISTRIBUTION} and delta = {DELTA_DISTRIBUTION}. (Pan and Hanratty)")
    pdf.chapter_body(f"For Original d_v50:")
    pdf.chapter_body(f"Equation: d_max, original = a * d_v50, original")
    pdf.chapter_body(f"Calculation (FPS): d_max,original = {A_DISTRIBUTION} * {results['dv50_original_fps']:.6f} ft = {results['d_max_original_fps']:.6f} ft")
   #pdf.chapter_body(f"Calculation (FPS): d_max,original = {A_DISTRIBUTION} * {results['dv50_original_fps']:.6f} = {results['d_max_original_fps']:.6f} ft")
    pdf.chapter_body(f"Result: Maximum Droplet Size (Original d_max) = {results['d_max_original_fps'] * FT_TO_MICRON:.2f} {micron_unit_label_pdf} ({from_fps(results['d_max_original_fps'], 'length'):.6f} {len_unit_pdf})")
    pdf.ln(2) # Small line break for readability

    pdf.chapter_body(f"For Adjusted d_v50:")
    pdf.chapter_body(f"Equation: d_max, adjusted = a * d_v50,adjusted")
    pdf.chapter_body(f"Calculation (FPS): d_max,adjusted = {A_DISTRIBUTION} * {results['dv50_adjusted_fps']:.6f} ft = {results['d_max_adjusted_fps']:.6f} ft")
    pdf.chapter_body(f"Result: Maximum Droplet Size (Adjusted d_max) = {results['d_max_adjusted_fps'] * FT_TO_MICRON:.2f} {micron_unit_label_pdf} ({from_fps(results['d_max_adjusted_fps'], 'length'):.6f} {len_unit_pdf})")
    pdf.ln(5)

    # Step 6: Entrainment Fraction (E) Calculation
    
    pdf.chapter_body1("Step 6: Calculate Entrainment Fraction (E)")
    pdf.chapter_body(f"The entrainment fraction is calculated using the following correlation (Pan and Hanratty (2002)) :")
    image_y = pdf.get_y() - 15
    image_height = 50
    pdf.image("ent_eq1.png", x=pdf.w/2 - 50, y=image_y, w=100, h=image_height)
    pdf.set_y(image_y + image_height-15)

    image_y = pdf.get_y() - 15
    image_height = 50
    pdf.image("ent_eq2.png", x=pdf.w/2 - 40, y=image_y, w=65, h=image_height)
    pdf.set_y(image_y + image_height-15)


    pdf.chapter_body(f"Gas Velocity (UG): {u_Gs:.2f} m/s")
    pdf.chapter_body(f"Liquid Loading (WL): {inputs['Q_liquid_mass_flow_rate_input']:.2f} {mass_flow_unit_pdf}")
    pdf.chapter_body(f"Result: Entrainment Fraction (E) = {results['E_fraction']:.4f} (dimensionless)")
    pdf.chapter_body(f"Result: Total Entrained Liquid Mass Flow Rate = {results['Q_entrained_total_mass_flow_rate_si']:.4f} {mass_flow_unit_pdf}")
    pdf.chapter_body(f"Result: Total Entrained Liquid Volume Flow Rate = {results['Q_entrained_total_volume_flow_rate_si']:.6f} {vol_flow_unit_pdf}") # New total volume flow
    pdf.ln(20)

    # Step 7: Calculate F-factor and Effective Gas Velocity
    
    pdf.chapter_body1("Step 7: Calculate F-factor and Effective Gas Velocity in Separator")
    
    pdf.chapter_body(f"L/Di Ratio (L_to_ME / D_pipe): {results['L_over_Di']:.2f}")
    pdf.chapter_body(f"Inlet Device: {inputs['inlet_device']}")
    pdf.chapter_body(f"Perforated Plate Used: {'Yes' if inputs['perforated_plate_option'] else 'No'}")
    pdf.chapter_body(f"Calculated F-factor: {results['F_factor']:.3f}")
    pdf.chapter_body(f"Effective Gas Velocity = Gas Velocity * F-Factor")
    pdf.chapter_body(f"Effective Gas Velocity in Separator (V_g_effective_separator): {from_fps(results['V_g_effective_separator_fps'], 'velocity'):.2f} {vel_unit_pdf}")
    pdf.ln(5)

    # Step 8: Gas Gravity Separation Section Efficiency
    pdf.add_page()
    pdf.chapter_body1("Step 8: Gas Gravity Separation Section Efficiency")
   
    if inputs['separator_type'] == "Horizontal":
        pdf.chapter_body1(f"Separator Type: Horizontal")
        pdf.chapter_body(f"Gas Space Height (hg): {inputs['h_g_input']:.3f} {len_unit_pdf}")
        pdf.chapter_body(f"Effective Separation Length (Le): {inputs['L_e_input']:.3f} {len_unit_pdf}")
    else: # Vertical
        pdf.chapter_body1(f"Separator Type: Vertical")
        pdf.chapter_body(f"Separator Diameter: {inputs['D_separator_input']:.3f} {len_unit_pdf}")
    pdf.chapter_body(f"Overall Separation Efficiency of Gravity Section: {results['gravity_separation_efficiency']:.2%}")
    hrly = plot_data_after_gravity['total_entrained_mass_flow_rate_si']
    pdf.chapter_body(f"Total Entrained Liquid Mass Flow Rate After Gravity Settling: {plot_data_after_gravity['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit_pdf}, ({hrly*3600:.4f} kg/hr)")
    hrly = plot_data_after_gravity['total_entrained_volume_flow_rate_si']
    pdf.chapter_body(f"Total Entrained Liquid Volume Flow Rate After Gravity Settling: {plot_data_after_gravity['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit_pdf}, ({hrly*3600:.4f} m3/hr)")
    pdf.ln(5)

    # Debug print statement for PDF generation context
    print(f"PDF Gen: Length of gravity_details_table_data: {len(plot_data_after_gravity['gravity_details_table_data']) if plot_data_after_gravity and 'gravity_details_table_data' in plot_data_after_gravity else 'N/A'}")

    # Display detailed table for gravity separation
    if plot_data_after_gravity and plot_data_after_gravity['gravity_details_table_data']:
        pdf.set_font('Arial', 'B', 10)
        # Add a new page if the table might overflow
        if pdf.get_y() + 10 + (len(plot_data_after_gravity['gravity_details_table_data']) + 1) * 6 > pdf.page_break_trigger:
            pdf.add_page()
            pdf.chapter_title('3. Calculation Results (Continued)') # Add a continued title
            pdf.ln(5) # Some space after continued title
        

        pdf.add_table(
            headers=["Droplet Size (um)", "Vt (ft/s)", "Cd", "Re_p", "Flow Regime", "Time Settle (s)", "h_max_settle (ft)", "Efficiency (Edp)"],
            data=[
                [
                    f"{row_dict['dp_microns']:.0f}",
                    f"{row_dict['Vt_ftps']:.2f}",
                    f"{row_dict['Cd']:.4f}",
                    f"{row_dict['Re_p']:.1f}",
                    row_dict['Flow Regime'],
                    f"{row_dict['Time Settle (s)']:.3f}",
                    f"{row_dict['h_max_settle (ft)']:.3f}",
                    f"{row_dict['Edp']:.2%}"
                ] for row_dict in plot_data_after_gravity['gravity_details_table_data']
            ],
            col_widths=[30, 20, 15, 18, 25, 25, 30, 28], # Adjust these widths as needed
            title='Detailed Separation Performance in Gas Gravity Section'
        )
    else:
        pdf.chapter_body("Detailed droplet separation data for gravity section not available.")
    pdf.ln(5) # Add spacing after the table or message

    # Step 9: Mist Extractor Performance
    pdf.add_page()
    pdf.chapter_body1("Step 9: Mist Extractor Performance")
    
    pdf.chapter_body(f"Mist Extractor Type: {inputs['mist_extractor_type']}")
   #pdf.chapter_body(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} {pressure_unit_pdf}")
   #pdf.chapter_body(f"K-Deration Factor (from Table 3): {results['k_deration_factor']:.3f}")

    if inputs['mist_extractor_type'] == "Mesh Pad":
        pdf.chapter_body(f"Mesh Pad Type: {inputs['mesh_pad_type']}")
        pdf.chapter_body(f"Mesh Pad Thickness: {inputs['mesh_pad_thickness_in']:.2f} in")
        pdf.chapter_body(f"Wire Diameter: {results['mesh_pad_params']['wire_diameter_in']:.3f} in")
        pdf.chapter_body(f"Specific Surface Area: {results['mesh_pad_params']['specific_surface_area_ft2_ft3']:.1f} ft^2/ft^3")
      # pdf.chapter_body(f"Base K_s: {results['mesh_pad_params']['Ks_ft_sec']:.2f} ft/sec")
      # pdf.chapter_body(f"Liquid Load Capacity: {results['mesh_pad_params']['liquid_load_gal_min_ft2']:.2f} gal/min/ft^2")
    elif inputs['mist_extractor_type'] == "Vane-Type":
        pdf.chapter_body(f"Vane Type: {inputs['vane_type']}")
        pdf.chapter_body(f"Flow Direction: {inputs['vane_flow_direction']}")
        pdf.chapter_body(f"Number of Bends: {inputs['vane_num_bends']}")
        pdf.chapter_body(f"Vane Spacing: {inputs['vane_spacing_in']:.2f} in")
        pdf.chapter_body(f"Bend Angle: {inputs['vane_bend_angle_deg']:.1f} deg")
       #pdf.chapter_body(f"Base K_s (Upflow): {results['vane_type_params']['Ks_ft_sec_upflow']:.2f} ft/sec")
       #pdf.chapter_body(f"Base K_s (Horizontal): {results['vane_type_params']['Ks_ft_sec_horizontal']:.2f} ft/sec")
       #pdf.chapter_body(f"Liquid Load Capacity: {results['vane_type_params']['liquid_load_gal_min_ft2']:.2f} gal/min/ft^2")
    elif inputs['mist_extractor_type'] == "Cyclonic":
        pdf.chapter_body(f"Cyclone Type: {inputs['cyclone_type']}")
        pdf.chapter_body(f"Cyclone Diameter: {inputs['cyclone_diameter_in']:.2f} in")
        pdf.chapter_body(f"Cyclone Length: {inputs['cyclone_length_in']:.2f} in")
        pdf.chapter_body(f"Inlet Swirl Angle: {inputs['cyclone_swirl_angle_deg']:.1f} deg")
       #pdf.chapter_body(f"Base K_s: {results['cyclone_type_params']['Ks_ft_sec_bundle_face_area']:.2f} ft/sec")
       #pdf.chapter_body(f"Liquid Load Capacity: {results['cyclone_type_params']['liquid_load_gal_min_ft2_bundle_face_area']:.2f} gal/min/ft^2")

    pdf.chapter_body(f"Overall Separation Efficiency of Mist Extractor: {results['mist_extractor_separation_efficiency']:.2%}")
    pdf.chapter_body(f"Total Entrained Liquid Mass Flow Rate After Mist Extractor: {plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit_pdf}")
    pdf.chapter_body(f"Total Entrained Liquid Volume Flow Rate After Mist Extractor: {plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit_pdf}")
    pdf.ln(5)

    # Display detailed table for gravity separation
    if plot_data_after_mist_extractor and plot_data_after_mist_extractor['mist_extractor_details_table_data']:
        pdf.set_font('Arial', 'B', 10)
        # Add a new page if the table might overflow
        if pdf.get_y() + 10 + (len(plot_data_after_mist_extractor['mist_extractor_details_table_data']) + 1) * 6 > pdf.page_break_trigger:
            pdf.add_page()
            pdf.chapter_title('3. Calculation Results (Continued)') # Add a continued title
            pdf.ln(5) # Some space after continued title

        if inputs['mist_extractor_type'] == "Mesh Pad":
            pdf.add_table(
                headers=["Droplet Size (um)", "Ew", "Epad"],
                data=[
                    [
                        f"{row_dict['dp_microns']:.0f}",
                        f"{row_dict['Ew']:.2%}",
                        f"{row_dict['Epad']:.2%}"
                    ] for row_dict in plot_data_after_mist_extractor['mist_extractor_details_table_data']
                ],
                col_widths=[30, 20, 15], # Adjust these widths as needed
                title='Detailed Separation Performance in Mist Extractor Section'
            )
        elif inputs['mist_extractor_type'] == "Vane-Type":
            pdf.add_table(
                headers=["Droplet Size (um)", "E_vane"],
                data=[
                    [
                        f"{row_dict['dp_microns']:.0f}",
                        f"{row_dict['Evane']:.2%}"
                    ] for row_dict in plot_data_after_mist_extractor['mist_extractor_details_table_data']
                ],
                col_widths=[30, 20], # Adjust these widths as needed
                title='Detailed Separation Performance in Mist Extractor Section'
            )
        elif inputs['mist_extractor_type'] == "Cyclonic":
            pdf.add_table(
                headers=["Droplet Size (um)", "Stk", "E_cyclone"],
                data=[
                    [
                        f"{row_dict['dp_microns']:.0f}",
                        f"{row_dict['Stk']:.1f}",
                        f"{row_dict['E_cycl']:.2%}"
                    ] for row_dict in plot_data_after_mist_extractor['mist_extractor_details_table_data']
                ],
                col_widths=[30, 30, 30], # Adjust these widths as needed
                title='Detailed Separation Performance in Mist Extractor Section'
            )    
    else:
        pdf.chapter_body("Detailed droplet separation data for Mist Extractor section not available.")
    pdf.ln(5) # Add spacing after the table or message
    
    pdf.chapter_body1("Final Carry-Over from Separator Outlet:")
    
    pdf.chapter_body(f"Total Carry-Over Mass Flow Rate: {plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit_pdf}")
    pdf.chapter_body(f"Total Carry-Over Volume Flow Rate: {plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit_pdf}")
    pdf.ln(5)


    # --- Droplet Distribution Plots ---
    pdf.add_page() # Start a new page for the plots
    pdf.chapter_title('4. Droplet Distribution Results')

    pdf.chapter_body("The following graphs show the calculated entrainment droplet size distribution:")

    # Original Distribution Plot
    pdf.set_font('Arial', 'B', 10)
    pdf.cell(0, 7, '4.1. Distribution at Inlet', 0, 1, 'L')
    pdf.ln(2)
    if plot_image_buffer_original:
        #pdf.image(plot_image_buffer_original, x=10, y=pdf.get_y(), w=pdf.w - 20)
        import tempfile
        # ...inside generate_pdf_report...
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
            tmpfile.write(plot_image_buffer_original.getbuffer())
            tmpfile.flush()
            pdf.image(tmpfile.name, x=10, y=pdf.get_y(), w=pdf.w - 20)
    pdf.ln(5)

    # Adjusted Distribution Plot (after inlet device)
    pdf.add_page() # Ensure the second plot is on a new page
    pdf.set_font('Arial', 'B', 10)
    pdf.cell(0, 7, '4.2. Distribution after Inlet Device', 0, 1, 'L')
    pdf.ln(2)
    if plot_image_buffer_adjusted:
        #pdf.image(plot_image_buffer_adjusted, x=10, y=pdf.get_y(), w=pdf.w - 20)
        import tempfile
        # ...inside generate_pdf_report...
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
            tmpfile.write(plot_image_buffer_adjusted.getbuffer())
            tmpfile.flush()
            pdf.image(tmpfile.name, x=10, y=pdf.get_y(), w=pdf.w - 20)
    pdf.ln(5)

    # Distribution After Gravity Settling Plot
    pdf.add_page() # Ensure this plot is on a new page
    pdf.set_font('Arial', 'B', 10)
    pdf.cell(0, 7, '4.3. Distribution after Gravity Settling Section', 0, 1, 'L')
    pdf.ln(2)
    # Generate and add plot for after gravity settling
    if plot_data_after_gravity and plot_data_after_gravity['dp_values_ft'].size > 0:
        fig_after_gravity, ax_after_gravity = plt.subplots(figsize=(10, 6))
        dp_values_microns_after_gravity = plot_data_after_gravity['dp_values_ft'] * FT_TO_MICRON

        ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
       #ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_after_gravity.set_xlabel(f'Droplet Size ({micron_unit_label_pdf})', fontsize=12)
        ax_after_gravity.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_after_gravity.tick_params(axis='y', labelcolor='black')
        ax_after_gravity.set_ylim(0, 1.2)
        ax_after_gravity.set_xlim(0, max(dp_values_microns_after_gravity) * 1.0 if dp_values_microns_after_gravity.size > 0 else 1000)

        ax2_after_gravity = ax_after_gravity.twinx()
        ax2_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
        ax2_after_gravity.set_ylabel('Volume Fraction', color='black', fontsize=12)
        ax2_after_gravity.tick_params(axis='y', labelcolor='black')
        max_norm_fv_after_gravity = max(plot_data_after_gravity['volume_fraction']) if plot_data_after_gravity['volume_fraction'].size > 0 else 0.1
        ax2_after_gravity.set_ylim(0, max_norm_fv_after_gravity * 1.2)

        lines_after_gravity, labels_after_gravity = ax_after_gravity.get_legend_handles_labels()
        lines2_after_gravity, labels2_after_gravity = ax2_after_gravity.get_legend_handles_labels()
        ax2_after_gravity.legend(lines_after_gravity + lines2_after_gravity, labels_after_gravity + labels2_after_gravity, loc='upper left', fontsize=10)
        ax_after_gravity.axhline(y=1, color='b', linestyle='--')

        plt.title('Entrainment Droplet Size Distribution (after Gravity Settling)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)

        buf_after_gravity = io.BytesIO()
        fig_after_gravity.savefig(buf_after_gravity, format="png", dpi=300)
        buf_after_gravity.seek(0)
        #pdf.image(buf_after_gravity, x=10, y=pdf.get_y(), w=pdf.w - 20)
        import tempfile
        # ...inside generate_pdf_report...
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
            tmpfile.write(buf_after_gravity.getbuffer())
            tmpfile.flush()
            pdf.image(tmpfile.name, x=10, y=pdf.get_y(), w=pdf.w - 20)
        pdf.ln(5)
        plt.close(fig_after_gravity) # Close the plot to free memory
    else:
        pdf.chapter_body("No data available for distribution after gravity settling. Please check your input parameters.")

    # Distribution After Mist Extractor Plot
    pdf.add_page() # Ensure this plot is on a new page
    pdf.set_font('Arial', 'B', 10)
    pdf.cell(0, 7, '4.4. Distribution after Mist Extractor', 0, 1, 'L')
    pdf.ln(2)
    if plot_data_after_mist_extractor and plot_data_after_mist_extractor['dp_values_ft'].size > 0:
        fig_after_me, ax_after_me = plt.subplots(figsize=(10, 6))
        dp_values_microns_after_me = plot_data_after_mist_extractor['dp_values_ft'] * FT_TO_MICRON

        ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
       #ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_after_me.set_xlabel(f'Droplet Size ({micron_unit_label_pdf})', fontsize=12)
        ax_after_me.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_after_me.tick_params(axis='y', labelcolor='black')
        ax_after_me.set_ylim(0, 1.2)
        ax_after_me.set_xlim(0, max(dp_values_microns_after_me) * 0.5 if dp_values_microns_after_me.size > 0 else 1000)

        ax2_after_me = ax_after_me.twinx()
        ax2_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['volume_fraction'], 'o-', label='Volume Fraction', markersize=2, color='#2ca02c')
        ax2_after_me.set_ylabel('Volume Fraction', color='black', fontsize=12)
        ax2_after_me.tick_params(axis='y', labelcolor='black')
        max_norm_fv_after_me = max(plot_data_after_mist_extractor['volume_fraction']) if plot_data_after_mist_extractor['volume_fraction'].size > 0 else 0.1
        ax2_after_me.set_ylim(0, max_norm_fv_after_me * 1.2)

        lines_after_me, labels_after_me = ax_after_me.get_legend_handles_labels()
        lines2_after_me, labels2_after_me = ax2_after_me.get_legend_handles_labels()
        ax2_after_me.legend(lines_after_me + lines2_after_me, labels_after_me + labels2_after_me, loc='upper left', fontsize=10)
        ax_after_me.axhline(y=1, color='b', linestyle='--')

        plt.title('Entrainment Droplet Size Distribution (after Mist Extractor)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)

        buf_after_me = io.BytesIO()
        fig_after_me.savefig(buf_after_me, format="png", dpi=300)
        buf_after_me.seek(0)
        #pdf.image(buf_after_me, x=10, y=pdf.get_y(), w=pdf.w - 20)
        import tempfile
        # ...inside generate_pdf_report...
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
            tmpfile.write(buf_after_me.getbuffer())
            tmpfile.flush()
            pdf.image(tmpfile.name, x=10, y=pdf.get_y(), w=pdf.w - 20)
        pdf.ln(5)
        plt.close(fig_after_me) # Close the plot to free memory
    else:
        pdf.chapter_body("No data available for distribution after mist extractor. Please check your input parameters.")


    # --- Volume Fraction Data Tables ---
    pdf.add_page() # Start a new page for the tables
    pdf.chapter_title('5. Volume Fraction Data Tables')

    # Original Data Table
    if plot_data_original and 'dp_values_ft' in plot_data_original and len(plot_data_original['dp_values_ft']) > 0:
        headers = ["Droplet Size (um)", "Volume Fraction", "Cumulative Undersize", "Entrained Mass Flow (kg/s)", "Entrained Volume Flow (m³/s)"]

        # Original Data Table
        full_data_original = []
        for i in range(len(plot_data_original['dp_values_ft'])): # Iterate using 'dp_values_ft'
            full_data_original.append([
                f"{plot_data_original['dp_values_ft'][i] * FT_TO_MICRON:.2f}", # Convert to microns here
                f"{plot_data_original['volume_fraction'][i]:.4f}",
                f"{plot_data_original['cumulative_volume_undersize'][i]:.4f}",
                f"{plot_data_original['entrained_mass_flow_rate_per_dp'][i]:.6f}",
                f"{plot_data_original['entrained_volume_flow_rate_per_dp'][i]:.9f}"
            ])
        col_widths = [28, 28, 35, 42, 48]
        pdf.add_table(headers, full_data_original, col_widths, title='5.1. Distribution at Inlet')
    else:
        pdf.chapter_body("No data available for original distribution table. Please check your input parameters.")

    # Adjusted Data Table
    pdf.add_page()
    if plot_data_adjusted and 'dp_values_ft' in plot_data_adjusted and len(plot_data_adjusted['dp_values_ft']) > 0:
        # Adjusted Data Table
        full_data_adjusted = []
        for i in range(len(plot_data_adjusted['dp_values_ft'])): # Iterate using 'dp_values_ft'
            full_data_adjusted.append([
                f"{plot_data_adjusted['dp_values_ft'][i] * FT_TO_MICRON:.2f}", # Convert to microns here
                f"{plot_data_adjusted['volume_fraction'][i]:.4f}",
                f"{plot_data_adjusted['cumulative_volume_undersize'][i]:.4f}",
                f"{plot_data_adjusted['entrained_mass_flow_rate_per_dp'][i]:.6f}",
                f"{plot_data_adjusted['entrained_volume_flow_rate_per_dp'][i]:.9f}"
            ])
        col_widths = [28, 28, 35, 42, 48]
        pdf.add_table(headers, full_data_adjusted, col_widths, title='5.2. Distribution after Inlet Device')
    else:
        pdf.chapter_body("No data available for adjusted distribution table. Please check your input parameters.")

    # Data Table After Gravity Settling
    pdf.add_page()
    if plot_data_after_gravity and 'dp_values_ft' in plot_data_after_gravity and len(plot_data_after_gravity['dp_values_ft']) > 0:
        full_data_after_gravity = []
        for i in range(len(plot_data_after_gravity['dp_values_ft'])):
            full_data_after_gravity.append([
                f"{plot_data_after_gravity['dp_values_ft'][i] * FT_TO_MICRON:.2f}",
                f"{plot_data_after_gravity['volume_fraction'][i]:.4f}",
                f"{plot_data_after_gravity['cumulative_volume_undersize'][i]:.4f}",
                f"{plot_data_after_gravity['entrained_mass_flow_rate_per_dp'][i]:.6f}",
                f"{plot_data_after_gravity['entrained_volume_flow_rate_per_dp'][i]:.9f}"
            ])
        col_widths = [28, 28, 35, 42, 48]
        pdf.add_table(headers, full_data_after_gravity, col_widths, title='5.3. Distribution after Gravity Settling Section')
    else:
        pdf.chapter_body("No data available for distribution table after gravity settling. Please check your input parameters.")

    # Data Table After Mist Extractor
    pdf.add_page()
    if plot_data_after_mist_extractor and 'dp_values_ft' in plot_data_after_mist_extractor and len(plot_data_after_mist_extractor['dp_values_ft']) > 0:
        full_data_after_me = []
        for i in range(len(plot_data_after_mist_extractor['dp_values_ft'])):
            full_data_after_me.append([
                f"{plot_data_after_mist_extractor['dp_values_ft'][i] * FT_TO_MICRON:.2f}",
                f"{plot_data_after_mist_extractor['volume_fraction'][i]:.4f}",
                f"{plot_data_after_mist_extractor['cumulative_volume_undersize'][i]:.4f}",
                f"{plot_data_after_mist_extractor['entrained_mass_flow_rate_per_dp'][i]:.6f}",
                f"{plot_data_after_mist_extractor['entrained_volume_flow_rate_per_dp'][i]:.9f}"
            ])
        col_widths = [28, 28, 35, 42, 48]
        pdf.add_table(headers, full_data_after_me, col_widths, title='5.4. Distribution After Mist Extractor')
    else:
        pdf.chapter_body("No data available for distribution table after mist extractor. Please check your input parameters.")

    # Meshpad Efficiency Plot
    pdf.add_page() # Ensure this plot is on a new page
    pdf.chapter_title('6. Mesh Pad Efficiency & Pressure Drop')
    if inputs['mist_extractor_type'] == "Mesh Pad":
        Q_gas_mass_lb_s = inputs['Q_gas_mass_flow_rate_input'] * 2.20462
        rho_g_fps = to_fps(inputs['rho_g_input'], "density")
        Q_gas_vol_ft3_s = Q_gas_mass_lb_s / rho_g_fps
        A_installed_ft2 = inputs['mesh_pad_A_installed_ft2']
        mesh_pad_params = MESH_PAD_PARAMETERS[inputs['mesh_pad_type']]
        mesh_pad_params = mesh_pad_params.copy()
        mesh_pad_params["thickness_in"] = inputs['mesh_pad_thickness_in']
        K_factor = mesh_pad_params["Ks_ft_sec"]
        K_dp = inputs['mesh_pad_K_dp']
        rho_l_fps = to_fps(inputs['rho_l_input'], "density")
        mu_g_fps = to_fps(inputs['mu_g_input'], "viscosity")
        fig = plot_mesh_pad_efficiency_with_pressure(
            Q_gas_vol_ft3_s, A_installed_ft2, rho_l_fps, rho_g_fps, mu_g_fps,
            mesh_pad_params, K_factor, K_dp, results
        )
             
        pdf.chapter_body(f"Actual face velocity before Mesh Pad: {results['V_face_actual_ft_s']:.2f} ft/s; ({results['V_face_actual_ft_s']*FT_TO_M:.2f} m/s)")
        pdf.chapter_body(f"Allowable velocity before Mesh Pad: {results['V_allow_ft_s']:.2f} ft/s; ({results['V_allow_ft_s']*FT_TO_M:.2f} m/s)")
        pdf.chapter_body(f"Mesh Pad Estimated Pressure Drop: {results['dp_mesh_pad_actual']:.2f} Pa; ({results['dp_mesh_pad_actual']/1000 :.2f} kPa)")
        if results['V_face_actual_ft_s'] > results['V_allow_ft_s']:
                carryover = results['carryover_percent'] / 100.0 * plot_data_after_gravity['total_entrained_volume_flow_rate_si'] * 3600.0 # Convert to hourly
                pdf.chapter_body(f"Estimated Carryover percent : {results['carryover_percent']:.1f} % ")
                pdf.chapter_body(f"Estimated Carryover : {results['carryover_percent']:.1f} % of Inlet flow rate {plot_data_after_gravity['total_entrained_volume_flow_rate_si']*3600:.3f} m3/hr = {carryover:.3f} m3/hr")
                pdf.chapter_body1("!!FLOODING RISK: OPERATING VELOCITY EXCEEDS ALLOWABLE VELOCITY FOR SELECTED MESHPAD GAS LOAD FACTOR (Ks).")     
        else:
                pdf.ln(0)
            #pdf.chapter_body(f"Operating Region within Mesh Pad Estimated Pressure Drop: {results['dp_mesh_pad_actual']:.2f} Pa; ({results['dp_mesh_pad_actual']/1000 :.2f} kPa)")

        pdf.chapter_body("The following graph shows the estimated efficiency and pressure drop across the Mesh Pad:")
        import tempfile
        buf_meshpad = io.BytesIO()
        fig.savefig(buf_meshpad, format="png", dpi=300)
        buf_meshpad.seek(0)
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
            tmpfile.write(buf_meshpad.getbuffer())
            tmpfile.flush()
            pdf.image(tmpfile.name, x=10, y=pdf.get_y(), w=pdf.w - 20)
        plt.close(fig)
        pdf.ln(5)


    return bytes(pdf.output(dest='S')) # Return PDF as bytes directly - Use this for GitHub
    #return pdf.output(dest='S').encode('latin1') # Use this for localhost testing
    #return pdf.output(dest='S')

# --- Streamlit App Layout ---

st.set_page_config(layout="wide", page_title="Oil & Gas Separation App")

tit1, tit2 = st.columns(2)
#with tit1:
    #st.header("🛢️ Oil and Gas Separator")
    #st.header("Oil and Gas Separator")
with tit2:
    st.markdown("<p style='text-align: right;'><b>🛢️ Liquid in Gas CarryOver Prediction (LICA Version 1.0)</b></p>", unsafe_allow_html=True)

#st.subheader("Liquid in Gas CarryOver Prediction")
st.markdown("---")
# --- Helper function to generate distribution data for a given dv50 and d_max ---

#@st.cache_data
def _generate_initial_distribution_data(dv50_value_fps, d_max_value_fps, num_points, E_fraction, Q_liquid_mass_flow_rate_si, rho_l_input_si):
    """
    Generates initial particle size distribution data (volume fraction, cumulative, entrained flow)
    for a given dv50, d_max, and other flow parameters.
    """
    plot_data = {
        'dp_values_ft': np.array([]),
        'volume_fraction': np.array([]),
        'cumulative_volume_undersize': np.array([]),
        'cumulative_volume_oversize': np.array([]),
        'entrained_mass_flow_rate_per_dp': np.array([]),
        'entrained_volume_flow_rate_per_dp': np.array([]),
        'total_entrained_mass_flow_rate_si': 0.0,
        'total_entrained_volume_flow_rate_si': 0.0
    }

    # Add a check for valid dv50 and d_max values
    if dv50_value_fps <= 0 or d_max_value_fps <= 0 or d_max_value_fps < dv50_value_fps:
        st.warning(f"Warning: Invalid droplet size range calculated (dv50: {dv50_value_fps:.6f} ft, d_max: {d_max_value_fps:.6f} ft). Check input parameters.")
        return plot_data # Return empty data

    dp_min_calc_fps = dv50_value_fps * 0.01
    dp_max_calc_fps = d_max_value_fps * 0.999

    if dp_min_calc_fps >= dp_max_calc_fps:
        st.warning("Calculated min droplet size is greater than or equal to max droplet size. Distribution cannot be generated.")
        return plot_data # Return empty data
        
    dp_values_ft = np.geomspace(dp_min_calc_fps, dp_max_calc_fps, num_points)

    volume_fraction_pdf_values = []

    for dp in dp_values_ft:
        if dp >= d_max_value_fps:
            z_val = np.inf
        else:
            # Ensure the argument to log is positive and finite
            log_arg = (A_DISTRIBUTION * dp) / (d_max_value_fps - dp)
            if log_arg <= 0: # Handle cases where dp is too small or d_max_value_fps - dp is too small
                z_val = -np.inf
            else:
                z_val = np.log(log_arg)

        # Handle potential division by zero for fv_dp if dp or (d_max_value_fps - dp) is zero
        if dp == 0 or (d_max_value_fps - dp) == 0:
            fv_dp = 0
        else:
           #denominator = np.sqrt(np.pi * dp * (d_max_value_fps - dp))
            denominator = np.sqrt(np.pi) * dp * (d_max_value_fps - dp) # reviewed 12.08.25
            
            if denominator == 0:
                fv_dp = 0
            else:
                fv_dp = ((DELTA_DISTRIBUTION * d_max_value_fps) / denominator) * np.exp(-DELTA_DISTRIBUTION**2 * z_val**2)

        volume_fraction_pdf_values.append(fv_dp)

    volume_fraction_pdf_values_array = np.array(volume_fraction_pdf_values)
    sum_of_pdf_values = np.sum(volume_fraction_pdf_values_array)

    normalized_volume_fraction = np.zeros_like(volume_fraction_pdf_values_array)
    if sum_of_pdf_values > 1e-9:
        normalized_volume_fraction = volume_fraction_pdf_values_array / sum_of_pdf_values

    cumulative_volume_undersize = np.cumsum(normalized_volume_fraction)
    cumulative_volume_oversize = 1 - cumulative_volume_undersize

    plot_data['dp_values_ft'] = dp_values_ft
    plot_data['volume_fraction'] = normalized_volume_fraction
    plot_data['cumulative_volume_undersize'] = cumulative_volume_undersize
    plot_data['cumulative_volume_oversize'] = cumulative_volume_oversize

    # Calculate total entrained liquid mass and volume flow rates
    Q_entrained_total_volume_flow_rate_si = E_fraction * Q_liquid_mass_flow_rate_si / rho_l_input_si
    Q_entrained_total_mass_flow_rate_si = 0.0
    if rho_l_input_si > 0:
        Q_entrained_total_mass_flow_rate_si = Q_entrained_total_volume_flow_rate_si * rho_l_input_si

    # Calculate entrained mass flow rate per droplet size interval using the normalized volume fraction
    plot_data['entrained_mass_flow_rate_per_dp'] = [
        fv_norm * Q_entrained_total_mass_flow_rate_si for fv_norm in normalized_volume_fraction
    ]

    # Calculate entrained volume flow rate per droplet size interval
    plot_data['entrained_volume_flow_rate_per_dp'] = [
        fv_norm * Q_entrained_total_volume_flow_rate_si for fv_norm in normalized_volume_fraction
    ]
    plot_data['total_entrained_mass_flow_rate_si'] = Q_entrained_total_mass_flow_rate_si
    plot_data['total_entrained_volume_flow_rate_si'] = Q_entrained_total_volume_flow_rate_si

    return plot_data

#@st.cache_data
def _calculate_and_apply_separation(
    initial_plot_data,
    _separation_stage_efficiency_func=None, # Function to apply for separation
    is_gravity_stage=False, # Flag to indicate if this is the gravity separation stage
    V_g_eff_sep_fps=0.0, # Required for gravity and mist extractor
    rho_l_fps=0.0, rho_g_fps=0.0, mu_g_fps=0.0, # Required for terminal velocity calc
    h_g_sep_fps=0.0, L_e_sep_fps=0.0, # Required for horizontal gravity
    separator_type="Horizontal", # Required for gravity
    mist_extractor_type=None, # Added to differentiate mist extractor types
    **kwargs_for_efficiency_func # Arguments for the efficiency function
):
    """
    Applies a separation efficiency function to an existing droplet distribution
    and calculates the new entrained mass/volume flow rates.
    If is_gravity_stage is True, it also collects detailed per-droplet data.
    """
    if not initial_plot_data or not initial_plot_data['dp_values_ft'].size > 0:
        return {
            'dp_values_ft': np.array([]),
            'volume_fraction': np.array([]),
            'cumulative_volume_undersize': np.array([]),
            'cumulative_volume_oversize': np.array([]),
            'entrained_mass_flow_rate_per_dp': np.array([]),
            'entrained_volume_flow_rate_per_dp': np.array([]),
            'total_entrained_mass_flow_rate_si': 0.0,
            'total_entrained_volume_flow_rate_si': 0.0,
            'overall_separation_efficiency': 0.0,
            'gravity_details_table_data': [], # Added for detailed gravity data
            'mist_extractor_details_table_data': [] # Added for detailed mist extractor data
        }

    dp_values_ft = initial_plot_data['dp_values_ft']
    initial_volume_fraction = initial_plot_data['volume_fraction']
    initial_entrained_mass_flow_rate_per_dp = initial_plot_data['entrained_mass_flow_rate_per_dp']
    initial_entrained_volume_flow_rate_per_dp = initial_plot_data['entrained_volume_flow_rate_per_dp']

    separated_entrained_mass_flow_rate_per_dp = np.zeros_like(initial_entrained_mass_flow_rate_per_dp)
    separated_entrained_volume_flow_rate_per_dp = np.zeros_like(initial_entrained_volume_flow_rate_per_dp)

    # Calculate initial total entrained flow rates from the provided initial_plot_data
    initial_total_entrained_mass_flow_rate_si = np.sum(initial_entrained_mass_flow_rate_per_dp)
    initial_total_entrained_volume_flow_rate_si = np.sum(initial_entrained_volume_flow_rate_per_dp)

    gravity_details_table_data = [] # To store details for Step 8 table
    mist_extractor_details_table_data = [] # To store details for mist extractor table
    

    # Apply separation efficiency for each droplet size
    for i, dp in enumerate(dp_values_ft):
        efficiency = 0.0
        Vt = 0.0
        Cd = 0.0
        Re_p = 0.0
        flow_regime = "N/A"
        time_settle = 0.0 # Time for droplet to fall h_g or L_e
        h_max_settle = 0.0 # Max height droplet can fall in gas residence time (horizontal) or effective separation height (vertical)

        if _separation_stage_efficiency_func:
            # For gravity stage, we need more detailed returns from the efficiency function
            if is_gravity_stage:
                if separator_type == "Horizontal":
                    efficiency, Vt, Cd, Re_p, h_max_settle_calc = gravity_efficiency_func_horizontal(
                        dp_fps=dp, V_g_eff_sep_fps=V_g_eff_sep_fps, h_g_sep_fps=h_g_sep_fps,
                        L_e_sep_fps=L_e_sep_fps, rho_l_fps=rho_l_fps, rho_g_fps=rho_g_fps, mu_g_fps=mu_g_fps
                    )
                    # Time for droplet to fall h_g
                    if Vt > 1e-9: # Avoid division by zero
                        time_settle = h_g_sep_fps / Vt
                    else:
                        time_settle = float('inf')
                    h_max_settle = h_max_settle_calc # This is the h_max_settle calculated within the function
                else: # Vertical
                    efficiency, Vt, Cd, Re_p = gravity_efficiency_func_vertical(
                        dp_fps=dp, V_g_eff_sep_fps=V_g_eff_sep_fps, rho_l_fps=rho_l_fps,
                        rho_g_fps=rho_g_fps, mu_g_fps=mu_g_fps
                    )
                    # For vertical, h_max_settle is effectively the height of the gas gravity section if separated
                    # Time for droplet to fall L_e (gas gravity section height)
                    if Vt > 1e-9:
                        time_settle = L_e_sep_fps / Vt
                    else:
                        time_settle = float('inf')
                    h_max_settle = L_e_sep_fps if efficiency > 0 else 0.0 # If separated, it effectively settles through L_e

                # Determine flow regime
                if Re_p < 2:
                    flow_regime = "Stokes'"
                elif 2 <= Re_p <= 500:
                    flow_regime = "Intermediate"
                else:
                    flow_regime = "Newton's"

                gravity_details_table_data.append({
                    "dp_microns": dp * FT_TO_MICRON,
                    "Vt_ftps": Vt,
                    "Cd": Cd,
                    "Re_p": Re_p,
                    "Flow Regime": flow_regime,
                    "Time Settle (s)": time_settle,
                    "h_max_settle (ft)": h_max_settle,
                    "Edp": efficiency # Individual droplet efficiency
                })

            else:  
                # For mist extractor stage (collect details based on type)
                if _separation_stage_efficiency_func == mesh_pad_efficiency_func:
                   efficiency, Stk, Ew = _separation_stage_efficiency_func(
                        dp_fps=dp,
                        V_g_eff_sep_fps=V_g_eff_sep_fps,
                        rho_l_fps=rho_l_fps,
                        rho_g_fps=rho_g_fps,
                        mu_g_fps=mu_g_fps,
                        **kwargs_for_efficiency_func
                   )
                   E_pad = efficiency # For mesh pad, efficiency is E_pad
                   mist_extractor_details_table_data.append({
                        "dp_microns": dp * FT_TO_MICRON,
                        "Stokes Number": Stk,
                        "Ew": Ew,
                        "Epad": E_pad
                   })
                elif _separation_stage_efficiency_func == vane_type_efficiency_func:
                   efficiency = _separation_stage_efficiency_func(
                        dp_fps=dp,
                        V_g_eff_sep_fps=V_g_eff_sep_fps,
                        rho_l_fps=rho_l_fps,
                        rho_g_fps=rho_g_fps,
                        mu_g_fps=mu_g_fps,
                        **kwargs_for_efficiency_func
                   )
                   E_vane = efficiency # For vane type, efficiency is E_vane
                   mist_extractor_details_table_data.append({
                        "dp_microns": dp * FT_TO_MICRON,
                        "Evane": E_vane
                   })
                elif _separation_stage_efficiency_func == demisting_cyclone_efficiency_func:
                   efficiency, Stk_cycl = _separation_stage_efficiency_func(
                        dp_fps=dp,
                        V_g_eff_sep_fps=V_g_eff_sep_fps,
                        rho_l_fps=rho_l_fps,
                        rho_g_fps=rho_g_fps,
                        mu_g_fps=mu_g_fps,
                        **kwargs_for_efficiency_func
                   )
                   E_cycl = efficiency
                   mist_extractor_details_table_data.append({
                        "dp_microns": dp * FT_TO_MICRON,
                        "Stk": Stk_cycl,
                        "E_cycl": E_cycl
                   })
 
            # Ensure efficiency is between 0 and 1
            efficiency = max(0.0, min(1.0, efficiency))

            separated_entrained_mass_flow_rate_per_dp[i] = initial_entrained_mass_flow_rate_per_dp[i] * (1 - efficiency)
            separated_entrained_volume_flow_rate_per_dp[i] = initial_entrained_volume_flow_rate_per_dp[i] * (1 - efficiency)
        else:
            # If no separation function, just carry over the initial values
            separated_entrained_mass_flow_rate_per_dp[i] = initial_entrained_mass_flow_rate_per_dp[i]
            separated_entrained_volume_flow_rate_per_dp[i] = initial_entrained_volume_flow_rate_per_dp[i]

    # Calculate new total entrained flow rates after this separation stage
    final_total_entrained_mass_flow_rate_si = np.sum(separated_entrained_mass_flow_rate_per_dp)
    final_total_entrained_volume_flow_rate_si = np.sum(separated_entrained_volume_flow_rate_per_dp)

    # Calculate overall separation efficiency for this stage
    overall_separation_efficiency = 0.0
    if initial_total_entrained_mass_flow_rate_si > 1e-9: # Avoid division by near-zero
        overall_separation_efficiency = 1.0 - (final_total_entrained_mass_flow_rate_si / initial_total_entrained_mass_flow_rate_si)

    # Recalculate normalized volume fraction based on the remaining entrained mass flow
    # This represents the *new* distribution of the *remaining* droplets
    new_volume_fraction = np.zeros_like(separated_entrained_volume_flow_rate_per_dp)
    if final_total_entrained_mass_flow_rate_si > 1e-9:
        new_volume_fraction = separated_entrained_volume_flow_rate_per_dp / final_total_entrained_volume_flow_rate_si

    new_cumulative_volume_undersize = np.cumsum(new_volume_fraction)
    new_cumulative_volume_oversize = 1 - new_cumulative_volume_undersize

    return {
        'dp_values_ft': dp_values_ft,
        'volume_fraction': new_volume_fraction, # This is the new normalized distribution of remaining droplets
        'cumulative_volume_undersize': new_cumulative_volume_undersize,
        'cumulative_volume_oversize': new_cumulative_volume_oversize,
        'entrained_mass_flow_rate_per_dp': separated_entrained_mass_flow_rate_per_dp,
        'entrained_volume_flow_rate_per_dp': separated_entrained_volume_flow_rate_per_dp,
        'total_entrained_mass_flow_rate_si': final_total_entrained_mass_flow_rate_si,
        'total_entrained_volume_flow_rate_si': final_total_entrained_volume_flow_rate_si,
        'overall_separation_efficiency': overall_separation_efficiency,
        'gravity_details_table_data': gravity_details_table_data, # Include detailed data for gravity stage
        'mist_extractor_details_table_data': mist_extractor_details_table_data # Include detailed data for mist extractor stage
    }


# --- Gravity Settling Efficiency Functions ---
def gravity_efficiency_func_horizontal(dp_fps, V_g_eff_sep_fps, h_g_sep_fps, L_e_sep_fps, rho_l_fps, rho_g_fps, mu_g_fps):
    """
    Calculates separation efficiency for a horizontal separator's gas gravity section.
    Assumes uniform droplet release over h_g.
    Returns: efficiency, Vt, Cd, Re_p, h_max_settle
    """
    if V_g_eff_sep_fps <= 0 or h_g_sep_fps <= 0 or L_e_sep_fps <= 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0 # No separation if no gas flow or no settling height/length

    V_t, Cd, Re_p = calculate_terminal_velocity(dp_fps, rho_l_fps, rho_g_fps, mu_g_fps)

    # Calculate maximum height from which a droplet can settle
    if V_g_eff_sep_fps > 1e-9: # Avoid division by near zero
        h_max_settle = (V_t * L_e_sep_fps) / V_g_eff_sep_fps
    else:
        h_max_settle = float('inf') # If gas velocity is zero, droplet can settle from infinite height

    # Efficiency is the fraction of h_g from which droplets of this size will settle
    efficiency = min(1.0, h_max_settle / h_g_sep_fps)
    return efficiency, V_t, Cd, Re_p, h_max_settle

def gravity_efficiency_func_vertical(dp_fps, V_g_eff_sep_fps, rho_l_fps, rho_g_fps, mu_g_fps):
    """
    Calculates separation efficiency for a vertical separator's gas gravity section.
    Sharp cutoff: 100% if Vt > V_g_eff_sep, 0% otherwise.
    Returns: efficiency, Vt, Cd, Re_p
    """
    V_t, Cd, Re_p = calculate_terminal_velocity(dp_fps, rho_l_fps, rho_g_fps, mu_g_fps)

    if V_t > V_g_eff_sep_fps:
        efficiency = 1.0 # Droplet settles
    else:
        efficiency = 0.0 # Droplet is carried over

    # For vertical, h_max_settle is conceptually the entire gravity section height if it settles
    # We'll handle this detail in the _calculate_and_apply_separation function for table consistency
    return efficiency, V_t, Cd, Re_p


# --- Function to perform all main calculations ---
def _perform_main_calculations(inputs):
    """Performs all scalar calculations and returns results dictionary."""
    results = {}

    D_pipe_si = inputs['D_input']
    Q_gas_mass_flow_rate_input_si = inputs['Q_gas_mass_flow_rate_input']
    rho_g_input_si = inputs['rho_g_input']
    Vol_flow_si = Q_gas_mass_flow_rate_input_si / rho_g_input_si if rho_g_input_si > 0 else 0.0
    Ug_si = Vol_flow_si / (0.785 * (D_pipe_si)**2) if D_pipe_si > 0 else 0.0
    Q_liquid_mass_flow_rate_input_si = inputs['Q_liquid_mass_flow_rate_input']
    rho_l_input_si = inputs['rho_l_input']
    #V_g_input_fps = Ug_si * M_TO_FT

    # Convert all SI inputs to FPS for consistent calculation
    D_pipe_fps = to_fps(inputs['D_input'], "length")
    rho_l_fps = to_fps(inputs['rho_l_input'], "density")
    mu_l_fps = to_fps(inputs['mu_l_input'], "viscosity")
    V_g_input_fps = Ug_si * M_TO_FT #to_fps(inputs['V_g_input'], "velocity") # Superficial gas velocity in feed pipe
    rho_g_fps = to_fps(inputs['rho_g_input'], "density")
    mu_g_fps = to_fps(inputs['mu_g_input'], "viscosity")
    sigma_fps = inputs['sigma_fps'] # This is already in poundal/ft

    # Separator specific inputs
    L_to_ME_fps = to_fps(inputs['L_to_ME_input'], 'length')
    D_separator_fps = to_fps(inputs['D_separator_input'], 'length')
    h_g_input_fps = to_fps(inputs['h_g_input'], 'length')
    L_e_input_fps = to_fps(inputs['L_e_input'], 'length')

    # Step 1: Calculate Superficial Gas Reynolds Number (Re_g) in feed pipe
    if mu_g_fps == 0: raise ValueError("Gas viscosity (μg) cannot be zero for Reynolds number calculation.")
    Re_g = (D_pipe_fps * V_g_input_fps * rho_g_fps) / mu_g_fps
    results['Re_g'] = Re_g

    # Step 2: Calculate Volume Median Diameter ($d_{v50}$) without inlet device effect
    if V_g_input_fps == 0 or rho_g_fps == 0 or rho_l_fps == 0 or mu_l_fps == 0:
        raise ValueError("Gas velocity, gas density, liquid density, and liquid viscosity must be non-zero for $d_{v50}$ calculation.")

    dv50_original_fps = 0.01 * (sigma_fps / (rho_g_fps * V_g_input_fps**2)) * (Re_g**(2/3)) * ((rho_g_fps / rho_l_fps)**(-1/3)) * ((mu_g_fps / mu_l_fps)**(2/3))
    results['dv50_original_fps'] = dv50_original_fps

    # Calculate d_max for the original distribution
    d_max_original_fps = A_DISTRIBUTION * dv50_original_fps
    results['d_max_original_fps'] = d_max_original_fps


    # Step 3: Determine Inlet Momentum (rho_m V_m^2)
    Q_gas_mass_lb_s = Q_gas_mass_flow_rate_input_si * 2.20462
    Q_liquid_mass_lb_s = Q_liquid_mass_flow_rate_input_si * 2.20462
    Q_gas_vol_ft3_s = Q_gas_mass_lb_s / rho_g_fps
    Q_liquid_vol_ft3_s = Q_liquid_mass_lb_s / rho_l_fps
        
    total_vol_flow = Q_gas_vol_ft3_s + Q_liquid_vol_ft3_s
    # Volume fractions
    alpha_g = Q_gas_vol_ft3_s / total_vol_flow
    alpha_l = 1 - alpha_g

    rho_mix_fps = alpha_g * rho_g_fps + alpha_l * rho_l_fps
    area_ft2 = math.pi * (D_pipe_fps ** 2) / 4.0
    V_mix_fps = total_vol_flow / area_ft2 if area_ft2 > 0 else 0.0

    rho_v_squared_fps = rho_mix_fps * V_mix_fps**2
    results['rho_v_squared_fps'] = rho_v_squared_fps

    # Step 4: Apply Inlet Device "Droplet Size Distribution Shift Factor"
    shift_factor = get_shift_factor(inputs['inlet_device'], rho_v_squared_fps)
    dv50_adjusted_fps = dv50_original_fps * shift_factor
    results['shift_factor'] = shift_factor
    results['dv50_adjusted_fps'] = dv50_adjusted_fps

    # Step 5: Calculate parameters for Upper-Limit Log Normal Distribution for adjusted dv50
    d_max_adjusted_fps = A_DISTRIBUTION * dv50_adjusted_fps
    results['d_max_adjusted_fps'] = d_max_adjusted_fps

    # Step 6: Entrainment Fraction (E) Calculation
   
    E_fraction = calculate_e_value(rho_l_fps, rho_g_fps, mu_l_fps, mu_g_fps, sigma_fps, D_pipe_fps,Ug_si, Q_liquid_mass_flow_rate_input_si)
        
    Q_entrained_total_mass_flow_rate_si = E_fraction * Q_liquid_mass_flow_rate_input_si

    Q_entrained_total_volume_flow_rate_si = 0.0
    if rho_l_input_si > 0:
        Q_entrained_total_volume_flow_rate_si = Q_entrained_total_mass_flow_rate_si / rho_l_input_si

    results['Wl_for_e_calc'] = Q_liquid_mass_flow_rate_input_si
    results['E_fraction'] = E_fraction
    results['Q_entrained_total_mass_flow_rate_si'] = Q_entrained_total_mass_flow_rate_si
    results['Q_entrained_total_volume_flow_rate_si'] = Q_entrained_total_volume_flow_rate_si

    # Step 7: Calculate F-factor and Effective Gas Velocity in Separator
    L_over_Di = L_to_ME_fps / D_separator_fps
    F_factor = get_f_factor(inputs['inlet_device'], L_over_Di, inputs['perforated_plate_option'])
    results['L_over_Di'] = L_over_Di
    results['F_factor'] = F_factor

    V_g_effective_separator_fps = 0.0
    if inputs['separator_type'] == "Vertical":
        # For vertical, gas velocity in separator is (Q_g_feed / A_separator_gas)
        # Q_g_feed = V_g_input_fps (feed pipe velocity) * A_pipe_fps
        #A_pipe_fps = np.pi * (D_pipe_fps / 2)**2
        A_separator_gas_vertical_fps = np.pi * (D_separator_fps / 2)**2
        if A_separator_gas_vertical_fps > 0:
            V_g_superficial_separator_fps = (Q_gas_vol_ft3_s) / A_separator_gas_vertical_fps
            V_g_effective_separator_fps = V_g_superficial_separator_fps * F_factor
        else:
            raise ValueError("Separator diameter cannot be zero for vertical separator gas velocity calculation.")
    else: # Horizontal
        # For horizontal, assume V_g_input is superficial velocity in separator gas section
        # and F_factor adjusts it.
        #############Start Added 11.08.2025
        R = D_separator_fps/2.0
        h_g = h_g_input_fps * D_separator_fps
        h_liq = D_separator_fps - h_g
        # clamp h_liq
        h = max(0.0, min(D_separator_fps, h_liq))
        if h == 0.0:
            A_liq = 0.0
        elif h == D:
            A_liq = math.pi * R**2
        else:
            A_liq = R**2 * math.acos((R - h)/R) - (R - h) * math.sqrt(2*R*h - h**2)
        A_g = math.pi * R**2 - A_liq
        V_g_superficial_separator_fps = (Q_gas_vol_ft3_s / A_g)
        #############End - Added 11.08.2025
        #V_g_superficial_separator_fps = V_g_input_fps # Assuming V_g_input is now superficial in separator for horizontal
        
                              
        V_g_effective_separator_fps = V_g_superficial_separator_fps * F_factor

    results['V_g_superficial_separator_fps'] = V_g_superficial_separator_fps
    results['V_g_effective_separator_fps'] = V_g_effective_separator_fps

    # Step 8: Gas Gravity Separation Section Efficiency
    # This will be used to calculate plot_data_after_gravity
    gravity_separation_efficiency = 0.0
    if inputs['separator_type'] == "Horizontal":
        if V_g_effective_separator_fps <= 0 or h_g_input_fps <= 0 or L_e_input_fps <= 0:
            st.warning("Invalid horizontal separator dimensions for gravity settling calculation. Efficiency set to 0.")
            gravity_separation_efficiency = 0.0
        else:
            # For reporting, calculate an average efficiency or a representative one
            # This will be the overall efficiency from the _calculate_and_apply_separation call
            pass # Calculated later in _calculate_and_apply_separation
    else: # Vertical
        if V_g_effective_separator_fps <= 0:
            st.warning("Invalid vertical separator gas velocity for gravity settling calculation. Efficiency set to 0.")
            gravity_separation_efficiency = 0.0
        else:
            # For reporting, calculate an average efficiency or a representative one
            pass # Calculated later in _calculate_and_apply_separation

    results['gravity_separation_efficiency'] = gravity_separation_efficiency # This will be updated after calling _calculate_and_apply_separation

    # Step 9: Mist Extractor Performance
    # Get K-deration factor based on pressure
    k_deration_factor = get_k_deration_factor(inputs['pressure_psig_input'])
    results['k_deration_factor'] = k_deration_factor

    mist_extractor_separation_efficiency = 0.0

    if inputs['mist_extractor_type'] == "Mesh Pad":
        mesh_pad_params = MESH_PAD_PARAMETERS[inputs['mesh_pad_type']]
        # Override thickness with user input
        mesh_pad_params_with_user_thickness = mesh_pad_params.copy()
        mesh_pad_params_with_user_thickness["thickness_in"] = inputs['mesh_pad_thickness_in']
        results['mesh_pad_params'] = mesh_pad_params_with_user_thickness # Store for reporting

        # The efficiency function will be called in _calculate_and_apply_separation
        # It needs rho_l_fps, rho_g_fps, mu_g_fps, V_g_effective_separator_fps, and mesh_pad_params_with_user_thickness
        pass # Efficiency calculated later

    elif inputs['mist_extractor_type'] == "Vane-Type":
        vane_type_params = VANE_PACK_PARAMETERS[inputs['vane_type']]
        # Override with user inputs for flow_direction, num_bends, spacing, angle
        vane_type_params_with_user_inputs = vane_type_params.copy()
        vane_type_params_with_user_inputs["flow_direction"] = inputs['vane_flow_direction']
        vane_type_params_with_user_inputs["number_of_bends"] = inputs['vane_num_bends']
        vane_type_params_with_user_inputs["vane_spacing_in"] = inputs['vane_spacing_in']
        vane_type_params_with_user_inputs["bend_angle_degree"] = inputs['vane_bend_angle_deg']
        results['vane_type_params'] = vane_type_params_with_user_inputs # Store for reporting

        pass # Efficiency calculated later

    elif inputs['mist_extractor_type'] == "Cyclonic":
        cyclone_type_params = CYCLONE_PARAMETERS[inputs['cyclone_type']]
        # Override with user inputs for diameter, length, swirl angle
        cyclone_type_params_with_user_inputs = cyclone_type_params.copy()
        cyclone_type_params_with_user_inputs["cyclone_inside_diameter_in"] = inputs['cyclone_diameter_in']
        cyclone_type_params_with_user_inputs["cyclone_length_in"] = inputs['cyclone_length_in']
        cyclone_type_params_with_user_inputs["inlet_swirl_angle_degree"] = inputs['cyclone_swirl_angle_deg']
        results['cyclone_type_params'] = cyclone_type_params_with_user_inputs # Store for reporting

        pass # Efficiency calculated later

    results['mist_extractor_separation_efficiency'] = mist_extractor_separation_efficiency # Updated after separation call

    return results


# Initialize session state for inputs and results if not already present
if 'inputs' not in st.session_state:
    st.session_state.inputs = {
        'D_input': 0.7, # m (1 ft) - Feed pipe diameter
        'rho_l_input': 730.7, # kg/m3 (40 lb/ft3)
        'mu_l_input': 0.000743, # Pa.s (0.0005 lb/ft-sec)
        #'V_g_input': 15.24, # m/s (50 ft/sec) - Superficial gas velocity in feed pipe
        'rho_g_input': 35.0, # kg/m3 (0.1 lb/ft3)
        'mu_g_input': 0.00001488, # Pa.s (0.00001 lb/ft-sec)
        'sigma_custom': 0.022, # Default N/m (from user request)
        'inlet_device': "No inlet device",
        'Q_liquid_mass_flow_rate_input': 180.0, # New input: kg/s (example value)
        'Q_gas_mass_flow_rate_input': 50.0, # New input: kg/s (example value)
        'gas_molecular_weight_input': 16.0, # g/mol 
        'num_points_distribution': 20, # Default number of points
        'separator_type': "Horizontal", # New input
        'h_g_input': 1.6, # m (for horizontal)
        'L_e_input': 12.0, # m (for horizontal)
        'D_separator_input': 2.0, # m (for vertical, or vessel diameter for horizontal)
        'L_to_ME_input': 1.0, # m (Length from Inlet Device to Mist Extractor)
        'perforated_plate_option': False, # New input
        'pressure_psig_input': 500.0, # psig (example value)
        'mist_extractor_type': "Mesh Pad", # New input
        'mesh_pad_type': "Standard mesh pad", # New input
        'mesh_pad_thickness_in': 6.0, # New input (default 6 inches)
        'vane_type': "Simple vane", # New input
        'vane_flow_direction': "Upflow", # New input
        'vane_num_bends': 5, # New input
        'vane_spacing_in': 0.75, # New input
        'vane_bend_angle_deg': 45.0, # New input
        'cyclone_type': "2.0 in. cyclones", # New input
        'cyclone_diameter_in': 2.0, # New input
        'cyclone_length_in': 10.0, # New input
        'cyclone_swirl_angle_deg': 45.0, # New input
    }
    # Initialize sigma_fps based on the new default sigma_custom
    st.session_state.inputs['sigma_fps'] = to_fps(st.session_state.inputs['sigma_custom'], "surface_tension")

if 'calculation_results' not in st.session_state:
    st.session_state.calculation_results = None
if 'plot_data_original' not in st.session_state:
    st.session_state.plot_data_original = None
if 'plot_data_adjusted' not in st.session_state:
    st.session_state.plot_data_adjusted = None
if 'plot_data_after_gravity' not in st.session_state:
    st.session_state.plot_data_after_gravity = None
if 'plot_data_after_mist_extractor' not in st.session_state:
    st.session_state.plot_data_after_mist_extractor = None
if 'report_date' not in st.session_state:
    st.session_state.report_date = ""


# Sidebar for navigation
#st.sidebar.markdown("## **Oil and Gas Separation Performance**")
#tit1, tit2 = st.columns(2)
#with tit1:
    #st.header("🛢️ Oil and Gas Separator")

st.sidebar.markdown("<h3 style='font-size: 20px;'><b>Oil and Gas Separator</b></h3>", unsafe_allow_html=True)


with st.sidebar:
    col1, col2, col3 = st.columns([3, 8, 3])  # Adjust the ratios as needed

    with col2:
        st.image("Sep.png", caption="", width=100)
        
st.sidebar.markdown("<h3 style='font-size: 16px;'><b>Liquid in Gas CarryOver Prediction</b></h3>", unsafe_allow_html=True)
#st.sidebar.markdown("---")

#with tit2:
#st.sidebar.markdown("<p style='font-size: 12px;'><b>🛢️ LICA Version 1.0</b></p>", unsafe_allow_html=True)

#st.subheader("Liquid in Gas CarryOver Prediction")


#st.sidebar.image("Sep.png", caption="", width=150)
st.sidebar.markdown("<h3 style='font-size: 16px;'><b>NAVIGATION BAR</b></h3>", unsafe_allow_html=True)
page = st.sidebar.radio("",["App Overview", "Input Parameters", "Flow Regime", "Calculation Results", "Carry Over Plots", "Summary of Results", "Generate Report"]) #, "References"])

import streamlit as st

st.sidebar.markdown("---")
st.sidebar.caption("Connect with me")
st.sidebar.caption("Feel free to connect with me on LinkedIn for further details or collaboration.")
st.sidebar.markdown("[![My LinkedIn Profile](https://img.shields.io/badge/LinkedIn-Profile-blue?style=for-the-badge&logo=linkedin)](https://www.linkedin.com/in/navaneethakrishnan-kannan/)")

# --- Page: App Overview ---
if page == "App Overview":
    import streamlit as st

    st.subheader("App Overview: Liquid in Gas Carry Over Prediction Calculator")

    st.markdown("""
    This application provides a **rigorous** and **transparent** framework for quantifying the performance of gas-liquid separators.
    It is built on a foundation of established first principles and validated empirical correlations, with methodologies detailed in the References section.
    """)

    st.subheader("Methodology")
    st.markdown("""
    While the **K-factor** remains a widely used approach for separator sizing, this tool employs a more granular **droplet settling velocity** analysis.
    This method provides a more detailed assessment of a separator's efficiency by modeling the behavior of entrained liquid droplets throughout the vessel.

    By shifting the focus from a single, aggregate factor to the **fundamental physics of droplet dynamics**, the application offers a more precise and defensible basis for design and performance evaluation.
    """)

    st.subheader("Key Sections and Functionality")
    st.markdown("""
    The app breaks down the separator's performance into a **multi-stage analysis**, allowing for a comprehensive evaluation of liquid carryover at each critical point:
    """)
    st.markdown("""
    - **Inlet Section**: Evaluates the initial droplet size distribution and the effectiveness of the inlet device in handling bulk liquid and large droplet separation.
    - **Gas Gravity Section**: Assesses the primary separation zone, where larger droplets are removed by gravity as the gas slows down.
    - **Mist Extractor Section**: Models the final stage of separation, analyzing the performance of the mist extractor (e.g., vanes, wire mesh) to capture sub-micron-sized droplets and minimize liquid carryover.
    """)

    st.subheader("Benefits and Applications")
    st.markdown("""
    This tool is designed for **process engineers**, **design specialists**, and **operations personnel** who require a more detailed and accurate method for:
    - Sizing separators
    - Troubleshooting performance issues
    - Validating vendor specifications

    For comprehensive details on the underlying methodology and data sources, please refer to the referenced articles.
    """)

    with st.expander("References", expanded=True):
        st.markdown("""
        - Bothamley, Mark. (2015). *Gas/Liquid Separators: Quantifying Separation Performance - Part 1*. Oil and Gas Facilities. 2. 21-29. [DOI: 10.2118/0813-0021-OGF](https://doi.org/10.2118/0813-0021-OGF)  
        - Bothamley, Mark. (2015). *Gas/Liquids Separators: Quantifying Separation Performance - Part 2*. Oil and Gas Facilities. 2. 35-47. [DOI: 10.2118/1013-0035-OGF](https://doi.org/10.2118/1013-0035-OGF)  
        - Pan, Lei & Hanratty, Thomas. (2002). *Correlation of entrainment for annular flow in horizontal pipes*. International Journal of Multiphase Flow. 28. 385-408. [DOI: 10.1016/S0301-9322(01)00074-X](https://doi.org/10.1016/S0301-9322%2801%2900074-X)  
        """)

    st.subheader("Disclaimer")
    st.markdown("""
    The purpose of this application is to serve as a **preliminary engineering tool** for evaluating gas-liquid separators.
    The results generated by this tool are based on the empirical correlations and methodologies outlined in the references.

    While this app is a powerful preliminary tool, **Computational Fluid Dynamics (CFD)** offers a more advanced method for analyzing complex liquid carryover and entrainment.

    ---

    **Important Considerations**:
    - **Correlation Accuracy**: The accuracy of the results depends on the range and conditions for which the original correlations were developed. Extrapolating beyond these conditions may lead to inaccuracies.
    - **Physical Properties**: The results are highly sensitive to the accuracy of the fluid physical properties (e.g., density, viscosity, surface tension) provided by the user.
    - **Not a Substitute for Professional Engineering**: Results should be reviewed and validated by an experienced professional before being used for final design, construction, or operational decisions.
    - **Limitation of Use**: This application should be used for educational and informational purposes only. The developers and contributors are not liable for any damages or losses resulting from the use of this software.

    For comprehensive details on the underlying methodology and data sources, please refer to the referenced articles.
    """)


# --- Page: Input Parameters ---
if page == "Input Parameters":
    st.subheader("1. Input Parameters (SI Units)")

    # Define unit labels for SI system
    len_unit = "m"
    dens_unit = "kg/m³"
    vel_unit = "m/s"
    visc_unit = "Pa·s"
    surf_tens_input_unit = "N/m"
    mass_flow_unit = "kg/s"
    pressure_unit = "psig"
    in_unit = "in" # For mist extractor dimensions

    st.markdown("#### **Separator Inlet Pipe Parameters**")
    #st.subheader("Separator Inlet Pipe Conditions")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.session_state.inputs['Q_liquid_mass_flow_rate_input'] = st.number_input(f"Liquid Mass Flow Rate ({mass_flow_unit})", min_value=0.0, value=st.session_state.inputs['Q_liquid_mass_flow_rate_input'], format="%.2f", key='Q_liquid_mass_flow_rate_input_widget',
                                help="The total mass flow rate of the liquid phase entering the system.")
        st.session_state.inputs['rho_l_input'] = st.number_input(f"Liquid Density ({dens_unit})", min_value=0.1, value=st.session_state.inputs['rho_l_input'], format="%.2f", key='rho_l_input_widget',
                                help="Density of the liquid phase.")
        st.session_state.inputs['mu_l_input'] = st.number_input(f"Liquid Viscosity ({visc_unit})", min_value=1e-8, value=st.session_state.inputs['mu_l_input'], format="%.8f", key='mu_l_input_widget',
                                help=f"Viscosity of the liquid phase.")
         # Construct the tooltip string from SURFACE_TENSION_TABLE_DYNE_CM
        tooltip_text = "Typical Liquid Surface Tension Values (N/m):\n"
        for fluid, value_dyne_cm in SURFACE_TENSION_TABLE_DYNE_CM.items():
            value_nm = value_dyne_cm * DYNE_CM_TO_NM
            tooltip_text += f"- {fluid}: {value_nm:.3f} N/m\n"

        st.session_state.inputs['sigma_custom'] = st.number_input(
            f"Liquid Surface Tension ({surf_tens_input_unit})",
            min_value=0.0001,
            value=st.session_state.inputs['sigma_custom'],
            format="%.4f",
            key='sigma_custom_input',
            help=tooltip_text # Use the constructed tooltip text here
        )
        st.session_state.inputs['sigma_fps'] = to_fps(st.session_state.inputs['sigma_custom'], "surface_tension")
    with col2:
        st.session_state.inputs['Q_gas_mass_flow_rate_input'] = st.number_input(f"Gas Mass Flow Rate ({mass_flow_unit})", min_value=0.0, value=st.session_state.inputs['Q_gas_mass_flow_rate_input'], format="%.2f", key='Q_gas_mass_flow_rate_input_widget',
                                help="The total mass flow rate of the gas phase entering the system.")
        st.session_state.inputs['rho_g_input'] = st.number_input(f"Gas Density ({dens_unit})", min_value=1e-5, value=st.session_state.inputs['rho_g_input'], format="%.5f", key='rho_g_input_widget',
                                help="Density of the gas phase.")
        st.session_state.inputs['mu_g_input'] = st.number_input(f"Gas Viscosity ({visc_unit})", min_value=1e-9, value=st.session_state.inputs['mu_g_input'], format="%.9f", key='mu_g_input_widget',
                                help=f"Viscosity of the gas phase.")
        st.session_state.inputs['gas_molecular_weight_input'] = st.number_input(f"Gas Molecular Weight", min_value=1.0, value=st.session_state.inputs['gas_molecular_weight_input'], format="%.2f", key='gas_molecular_weight_input_widget',
                                help=f"Gas Molecular Weight.")
                
    with col3:
        st.session_state.inputs['D_input'] = st.number_input(f"Pipe Inside Diameter ({len_unit})", min_value=0.001, value=st.session_state.inputs['D_input'], format="%.4f", key='D_input_widget',
                               help="Diameter of the feed pipe to the separator.")
       # st.session_state.inputs['V_g_input'] = st.number_input(f"Gas Velocity in Feed Pipe ({vel_unit})", min_value=0.01, value=st.session_state.inputs['V_g_input'], format="%.2f", key='V_g_input_widget',
        #                      help="Superficial gas velocity in the feed pipe.")
    
    st.markdown("---")
    st.markdown("#### **Separator Details**")
    col_id, col_st = st.columns(2)

    #with col_st:
        # Update sigma_fps based on the new default sigma_custom
        #st.info(f"**Current Liquid Surface Tension:** {st.session_state.inputs['sigma_custom']:.3f} {surf_tens_input_unit}")
    
    with col_id:
        #st.subheader("Separator Inlet Device")
        st.session_state.inputs['separator_type'] = st.radio(
            "Select Separator Type",
            options=["Horizontal", "Vertical"],
            index=0 if st.session_state.inputs['separator_type'] == "Horizontal" else 1,
            key='separator_type_radio'
        )
        st.markdown("---") 
        current_inlet_device_index = ["No inlet device", "Diverter plate", "Half-pipe", "Vane-type", "Cyclonic"].index(st.session_state.inputs['inlet_device'])
        st.session_state.inputs['inlet_device'] = st.radio(
            "Choose Inlet Device Type",
            options=["No inlet device", "Diverter plate", "Half-pipe", "Vane-type", "Cyclonic"],
            #index=0,
            index=current_inlet_device_index,
            key='inlet_device_select',
            help="The inlet device influences the droplet size distribution downstream."
        )
        st.markdown("---") 
        st.session_state.inputs['perforated_plate_option'] = st.checkbox("Use Perforated Plate (Flow Straightening)", value=st.session_state.inputs['perforated_plate_option'], key='perforated_plate_checkbox',
                                help="Check if a perforated plate is used to improve gas velocity profile.")
        
    with col_st:
        
        if st.session_state.inputs['separator_type'] == "Horizontal":
            st.session_state.inputs['h_g_input'] = st.number_input(f"Gas Space Height (h_g) ({len_unit})", min_value=0.01, value=st.session_state.inputs['h_g_input'], format="%.3f", key='h_g_input_widget',
                                    help="Vertical height of the gas phase in the horizontal separator.")
            st.session_state.inputs['L_e_input'] = st.number_input(f"Effective Separation Length (L_e) ({len_unit})", min_value=0.01, value=st.session_state.inputs['L_e_input'], format="%.3f", key='L_e_input_widget',
                                    help="Horizontal length available for gas-liquid separation in the horizontal separator.")
            st.session_state.inputs['D_separator_input'] = st.number_input(f"Horizontal Separator Vessel Diameter ({len_unit})", min_value=0.1, value=st.session_state.inputs['D_separator_input'], format="%.3f", key='D_separator_input_widget',
                                    help="Diameter of the horizontal separator vessel. Used for context, not directly in gravity settling calculations if h_g is provided.")
        else: # Vertical
            st.session_state.inputs['D_separator_input'] = st.number_input(f"Vertical Separator Diameter ({len_unit})", min_value=0.1, value=st.session_state.inputs['D_separator_input'], format="%.3f", key='D_separator_input_widget',
                                    help="Diameter of the vertical separator vessel.")
            st.session_state.inputs['L_e_input'] = st.number_input(f"Gas Gravity Section Height (L_e) ({len_unit})", min_value=0.01, value=st.session_state.inputs['L_e_input'], format="%.3f", key='L_e_input_widget',
                                    help="Vertical height of the gas gravity separation section in the vertical separator. This is used as 'h_g' for vertical settling.")
            # For vertical, h_g_input is effectively L_e_input for gravity calculations
            st.session_state.inputs['h_g_input'] = st.session_state.inputs['L_e_input']

        st.session_state.inputs['L_to_ME_input'] = st.number_input(f"Length from Inlet Device to Mist Extractor (L_to_ME) ({len_unit})", min_value=0.0, value=st.session_state.inputs['L_to_ME_input'], format="%.3f", key='L_to_ME_input_widget',
                                help="The distance from the inlet device outlet to the mist extractor. Used for F-factor calculation (L/Di).")
        st.session_state.inputs['pressure_psig_input'] = st.number_input(f"Operating Pressure ({pressure_unit})", min_value=0.0, value=st.session_state.inputs['pressure_psig_input'], format="%.1f", key='pressure_psig_input_widget',
                                help="Operating pressure of the separator. Used for K-deration of mist extractors.")

    st.markdown("---")
    st.markdown("###### **Mist Extractor Details**")
    col_id1, col_st1 = st.columns(2)
    
    with col_id1:
        st.session_state.inputs['mist_extractor_type'] = st.radio(
        "Select Mist Extractor Type",
        options=["Mesh Pad", "Vane-Type", "Cyclonic"],
        index=["Mesh Pad", "Vane-Type", "Cyclonic"].index(st.session_state.inputs['mist_extractor_type']),
        key='mist_extractor_type_select'
        )
    with col_st1:
        if st.session_state.inputs['mist_extractor_type'] == "Mesh Pad":
            current_mesh_pad_index = list(MESH_PAD_PARAMETERS.keys()).index(st.session_state.inputs['mesh_pad_type'])
            st.session_state.inputs['mesh_pad_type'] = st.selectbox(
            "Mesh Pad Type",
            options=list(MESH_PAD_PARAMETERS.keys()),
            index=current_mesh_pad_index,
            key='mesh_pad_type_select'
            )
            st.session_state.inputs['mesh_pad_thickness_in'] = st.number_input(
            f"Mesh Pad Thickness ({in_unit})",
            min_value=1.0, value=st.session_state.inputs['mesh_pad_thickness_in'], format="%.2f", key='mesh_pad_thickness_in_input',
            help="Typical thickness for mesh pads is around 6 inches."
            )
            # NEW: Installed mesh pad area (ft²)
            st.session_state.inputs['mesh_pad_A_installed_ft2'] = st.number_input(
                "Installed Mesh Pad Area (ft²)",
                min_value=0.1, value=st.session_state.inputs.get('mesh_pad_A_installed_ft2', 25.0), format="%.2f",
                key='mesh_pad_A_installed_ft2_input',
                help="Frontal area of the mesh pad installed (ft²)."
            )
            # NEW: K_dp (pressure drop factor)
            st.session_state.inputs['mesh_pad_K_dp'] = st.number_input(
                "Mesh Pad Pressure Drop Factor (K_dp)",
                min_value=0.01, max_value= 0.1, value=st.session_state.inputs.get('mesh_pad_K_dp', 0.1), format="%.2f",
                key='mesh_pad_K_dp_input',
                help="Pressure drop factor for mesh pad (dimensionless, typical 0.01–0.1)."
            )

        elif st.session_state.inputs['mist_extractor_type'] == "Vane-Type":
            current_vane_type_index = list(VANE_PACK_PARAMETERS.keys()).index(st.session_state.inputs['vane_type'])
            st.session_state.inputs['vane_type'] = st.selectbox(
            "Vane Type",
            options=list(VANE_PACK_PARAMETERS.keys()),
            index=current_vane_type_index,
            key='vane_type_select'
            )
            current_vane_flow_direction_index = ["Upflow", "Horizontal"].index(st.session_state.inputs['vane_flow_direction'])
            st.session_state.inputs['vane_flow_direction'] = st.selectbox(
            "Flow Direction",
            options=["Upflow", "Horizontal"],
            index=current_vane_flow_direction_index,
            key='vane_flow_direction_select'
            )
            st.session_state.inputs['vane_num_bends'] = st.number_input(
            "Number of Bends",
            min_value=1, max_value=10, value=st.session_state.inputs['vane_num_bends'], step=1, key='vane_num_bends_input',
            help="Typical range is 5-8 bends."
            )
            st.session_state.inputs['vane_spacing_in'] = st.number_input(
            f"Vane Spacing ({in_unit})",
            min_value=0.1, value=st.session_state.inputs['vane_spacing_in'], format="%.2f", key='vane_spacing_in_input',
            help="Typical range is 0.5-1 inch."
            )
            st.session_state.inputs['vane_bend_angle_deg'] = st.number_input(
            f"Bend Angle (degrees)",
            min_value=1.0, max_value=90.0, value=st.session_state.inputs['vane_bend_angle_deg'], format="%.1f", key='vane_bend_angle_deg_input',
            help="Typical range is 30-60 degrees, 45 degrees is most common."
            )
        elif st.session_state.inputs['mist_extractor_type'] == "Cyclonic":
            current_cyclone_type_index = list(CYCLONE_PARAMETERS.keys()).index(st.session_state.inputs['cyclone_type'])
            st.session_state.inputs['cyclone_type'] = st.selectbox(
            "Cyclone Type",
            options=list(CYCLONE_PARAMETERS.keys()),
            index=current_cyclone_type_index,
            key='cyclone_type_select'
            )
            st.session_state.inputs['cyclone_diameter_in'] = st.number_input(
            f"Cyclone Diameter ({in_unit})",
            min_value=0.1, value=st.session_state.inputs['cyclone_diameter_in'], format="%.2f", key='cyclone_diameter_in_input',
            help="Inside diameter of individual cyclone tube."
            )
            st.session_state.inputs['cyclone_length_in'] = st.number_input(
            f"Cyclone Length ({in_unit})",
            min_value=1.0, value=st.session_state.inputs['cyclone_length_in'], format="%.2f", key='cyclone_length_in_input',
            help="Length of individual cyclone tube."
            )
            st.session_state.inputs['cyclone_swirl_angle_deg'] = st.number_input(
            f"Inlet Swirl Angle (degrees)",
            min_value=1.0, max_value=90.0, value=st.session_state.inputs['cyclone_swirl_angle_deg'], format="%.1f", key='cyclone_swirl_angle_deg_input',
            help="Inlet swirl angle of the cyclone."
            )
    
    # When inputs on this page change, trigger recalculation for initial state
    import datetime
    st.session_state.report_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    try:
        # Only perform main scalar calculations here
        st.session_state.calculation_results = _perform_main_calculations(st.session_state.inputs)

        # Generate initial distribution data (after inlet device, before gravity settling)
        st.session_state.plot_data_original = _generate_initial_distribution_data(
            st.session_state.calculation_results['dv50_original_fps'],
            st.session_state.calculation_results['d_max_original_fps'],
            st.session_state.inputs['num_points_distribution'],
            st.session_state.calculation_results['E_fraction'],
            st.session_state.inputs['Q_liquid_mass_flow_rate_input'],
            st.session_state.inputs['rho_l_input']
        )

        st.session_state.plot_data_adjusted = _generate_initial_distribution_data(
            st.session_state.calculation_results['dv50_adjusted_fps'],
            st.session_state.calculation_results['d_max_adjusted_fps'],
            st.session_state.inputs['num_points_distribution'],
            st.session_state.calculation_results['E_fraction'],
            st.session_state.inputs['Q_liquid_mass_flow_rate_input'],
            st.session_state.inputs['rho_l_input']
        )

        # Calculate and apply gravity settling
        if st.session_state.inputs['separator_type'] == "Horizontal":
            st.session_state.plot_data_after_gravity = _calculate_and_apply_separation(
                st.session_state.plot_data_adjusted, # Input is the adjusted distribution
                _separation_stage_efficiency_func=gravity_efficiency_func_horizontal,
                is_gravity_stage=True,
                V_g_eff_sep_fps=st.session_state.calculation_results['V_g_effective_separator_fps'],
                h_g_sep_fps=to_fps(st.session_state.inputs['h_g_input'], 'length'),
                L_e_sep_fps=to_fps(st.session_state.inputs['L_e_input'], 'length'),
                rho_l_fps=to_fps(st.session_state.inputs['rho_l_input'], 'density'),
                rho_g_fps=to_fps(st.session_state.inputs['rho_g_input'], 'density'),
                mu_g_fps=to_fps(st.session_state.inputs['mu_g_input'], 'viscosity'),
                separator_type=st.session_state.inputs['separator_type']
            )
        else: # Vertical
            st.session_state.plot_data_after_gravity = _calculate_and_apply_separation(
                st.session_state.plot_data_adjusted, # Input is the adjusted distribution
                _separation_stage_efficiency_func=gravity_efficiency_func_vertical,
                is_gravity_stage=True,
                V_g_eff_sep_fps=st.session_state.calculation_results['V_g_effective_separator_fps'],
                # For vertical, h_g_input is effectively L_e_input for gravity calculations
                h_g_sep_fps=to_fps(st.session_state.inputs['L_e_input'], 'length'), 
                L_e_sep_fps=to_fps(st.session_state.inputs['L_e_input'], 'length'), # Pass L_e_input for vertical time_settle calc
                rho_l_fps=to_fps(st.session_state.inputs['rho_l_input'], 'density'),
                rho_g_fps=to_fps(st.session_state.inputs['rho_g_input'], 'density'),
                mu_g_fps=to_fps(st.session_state.inputs['mu_g_input'], 'viscosity'),
                separator_type=st.session_state.inputs['separator_type']
            )

        # Update the overall gravity separation efficiency in results for reporting
        if st.session_state.plot_data_after_gravity:
            st.session_state.calculation_results['gravity_separation_efficiency'] = st.session_state.plot_data_after_gravity['overall_separation_efficiency']
        else:
            st.session_state.calculation_results['gravity_separation_efficiency'] = 0.0

        # Calculate and apply mist extractor efficiency
        if st.session_state.plot_data_after_gravity and st.session_state.plot_data_after_gravity['dp_values_ft'].size > 0:
            if st.session_state.inputs['mist_extractor_type'] == "Mesh Pad":
                mesh_pad_params = MESH_PAD_PARAMETERS[st.session_state.inputs['mesh_pad_type']]
                mesh_pad_params_with_user_thickness = mesh_pad_params.copy()
                mesh_pad_params_with_user_thickness["thickness_in"] = st.session_state.inputs['mesh_pad_thickness_in']

                st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                    st.session_state.plot_data_after_gravity,
                    _separation_stage_efficiency_func=mesh_pad_efficiency_func,
                    V_g_eff_sep_fps=st.session_state.calculation_results['V_g_effective_separator_fps'],
                    rho_l_fps=to_fps(st.session_state.inputs['rho_l_input'], 'density'),
                    rho_g_fps=to_fps(st.session_state.inputs['rho_g_input'], 'density'),
                    mu_g_fps=to_fps(st.session_state.inputs['mu_g_input'], 'viscosity'),
                    #mist_extractor_type_str="Mesh Pad", # Pass the type string
                    mesh_pad_type_params_fps=mesh_pad_params_with_user_thickness
                )
            elif st.session_state.inputs['mist_extractor_type'] == "Vane-Type":
                vane_type_params = VANE_PACK_PARAMETERS[st.session_state.inputs['vane_type']]
                vane_type_params_with_user_inputs = vane_type_params.copy()
                vane_type_params_with_user_inputs["flow_direction"] = st.session_state.inputs['vane_flow_direction']
                vane_type_params_with_user_inputs["number_of_bends"] = st.session_state.inputs['vane_num_bends']
                vane_type_params_with_user_inputs["vane_spacing_in"] = st.session_state.inputs['vane_spacing_in']
                vane_type_params_with_user_inputs["bend_angle_degree"] = st.session_state.inputs['vane_bend_angle_deg']

                st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                    st.session_state.plot_data_after_gravity,
                    _separation_stage_efficiency_func=vane_type_efficiency_func,
                    V_g_eff_sep_fps=st.session_state.calculation_results['V_g_effective_separator_fps'],
                    rho_l_fps=to_fps(st.session_state.inputs['rho_l_input'], 'density'),
                    rho_g_fps=to_fps(st.session_state.inputs['rho_g_input'], 'density'),
                    mu_g_fps=to_fps(st.session_state.inputs['mu_g_input'], 'viscosity'),
                   # mist_extractor_type_str="Vane-Type", # Pass the type string
                    vane_type_params_fps=vane_type_params_with_user_inputs
                )
            elif st.session_state.inputs['mist_extractor_type'] == "Cyclonic":
                cyclone_type_params = CYCLONE_PARAMETERS[st.session_state.inputs['cyclone_type']]
                cyclone_type_params_with_user_inputs = cyclone_type_params.copy()
                cyclone_type_params_with_user_inputs["cyclone_inside_diameter_in"] = st.session_state.inputs['cyclone_diameter_in']
                cyclone_type_params_with_user_inputs["cyclone_length_in"] = st.session_state.inputs['cyclone_length_in']
                cyclone_type_params_with_user_inputs["inlet_swirl_angle_degree"] = st.session_state.inputs['cyclone_swirl_angle_deg']

                st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                    st.session_state.plot_data_after_gravity,
                    _separation_stage_efficiency_func=demisting_cyclone_efficiency_func,
                    V_g_eff_sep_fps=st.session_state.calculation_results['V_g_effective_separator_fps'],
                    rho_l_fps=to_fps(st.session_state.inputs['rho_l_input'], 'density'),
                    rho_g_fps=to_fps(st.session_state.inputs['rho_g_input'], 'density'),
                    mu_g_fps=to_fps(st.session_state.inputs['mu_g_input'], 'viscosity'),
                   # mist_extractor_type_str="Cyclonic", # Pass the type string
                    cyclone_type_params_fps=cyclone_type_params_with_user_inputs
                )
            else:
                st.session_state.plot_data_after_mist_extractor = st.session_state.plot_data_after_gravity # No mist extractor selected, so no change

            # Update the overall mist extractor separation efficiency in results for reporting
            if st.session_state.plot_data_after_mist_extractor:
                st.session_state.calculation_results['mist_extractor_separation_efficiency'] = st.session_state.plot_data_after_mist_extractor['overall_separation_efficiency']
                 # Store mist extractor details if available
                st.session_state.calculation_results['mist_extractor_details_table_data'] = st.session_state.plot_data_after_mist_extractor['mist_extractor_details_table_data']
            else:
                st.session_state.calculation_results['mist_extractor_separation_efficiency'] = 0.0

        else:
            st.session_state.plot_data_after_mist_extractor = None


    except Exception as e:
        st.error(f"An error occurred during calculation: {e}")
        st.session_state.calculation_results = None
        st.session_state.plot_data_original = None
        st.session_state.plot_data_adjusted = None
        st.session_state.plot_data_after_gravity = None
        st.session_state.plot_data_after_mist_extractor = None

# --- Page: Flow Regime ---
elif page == "Flow Regime":
        st.subheader("2. Flow Regime Map (Inlet)")  
        st.markdown("<h3 style='font-size: 24px;'><b>Taitel & Dukler Flow Regime Map</b></h3>", unsafe_allow_html=True)
        #st.markdown("**Taitel & Dukler Flow Regime Map**")

        st.markdown("---")
        rho_L = st.session_state.inputs['rho_l_input']
        mu_L = st.session_state.inputs['mu_l_input']
        rho_G = st.session_state.inputs['rho_g_input']
        mu_G = st.session_state.inputs['mu_g_input']
        Q_liquid_mass_flow_rate = st.session_state.inputs['Q_liquid_mass_flow_rate_input']
        Q_gas_mass_flow_rate = st.session_state.inputs['Q_gas_mass_flow_rate_input']
        D = st.session_state.inputs['D_input']  
        # Use st.columns to organize input widgets
        col1, col2 = st.columns(2)

        with col1:
            #rho_L = 800 #st.number_input("Liquid Density (ρ_L, kg/m³)", value=rho_l_input, min_value=0.01)
            #rho_G = 50 #st.number_input("Gas Density (ρ_G, kg/m³)", value=rho_g_fps, min_value=0.01)
            #mu_L = st.number_input("Liquid Viscosity (µ_L, Pa·s)", value=mu_L1, min_value=0.00001, format="%.6f")
            #mu_G = 0.00001 #st.number_input("Gas Viscosity (µ_G, Pa·s)", value=mu_g_fps, min_value=0.000001, format="%.8f")
            st.write(f"Liquid Density (ρ_L, kg/m³): {rho_L:.2f}")
            st.write(f"Liquid Viscosity (µ_L, Pa·s): {mu_L:.6f}")
            st.write(f"Gas Density (ρ_G, kg/m³): {rho_G:.2f}")
            st.write(f"Gas Viscosity (µ_G, Pa·s): {mu_G:.6f}")
            st.write(f"Liquid Mass Flow Rate (kg/s): {Q_liquid_mass_flow_rate:.2f}")
            st.write(f"Gas Mass Flow Rate (kg/s): {Q_gas_mass_flow_rate:.2f}")

        #with col2:
            #Q_liquid_mass_flow_rate = 0.5 #st.number_input("Liquid Mass Flow Rate (kg/s)", value= Q_liquid_mass_flow_rate, min_value=0.0)
            #Q_gas_mass_flow_rate = 1 #st.number_input("Gas Mass Flow Rate (kg/s)", value=Q_gas_mass_flow_rate, min_value=0.0)
            #D = 0.6 #st.number_input("Pipe Diameter (D, m)", value=D_pipe_fps, min_value=0.001)
            alpha = 0 #st.number_input("Pipe Inclination Angle (α, degrees)", value=0.0, min_value=-90.0, max_value=90.0)

        # Calculate superficial velocities from mass flow rates
        u_Ls = (Q_liquid_mass_flow_rate / rho_L) / (0.25 * math.pi * D**2) if rho_L > 0 and D > 0 else 0
        u_Gs = (Q_gas_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0

       
        # Perform calculations by calling functions from the flow_regime module
        nu_L = mu_L / rho_L
        X = fr.compute_X(rho_L, rho_G, mu_L, mu_G, u_Ls, u_Gs, D)
        F = fr.compute_F(rho_L, rho_G, u_Gs, D, alpha=alpha)
        K = fr.compute_K(F, D, u_Ls, nu_L)

        with col2:
            st.markdown("##### ---------Computed Parameters---------")
            st.write(f"**Superficial Liquid Velocity:** $U_{{Ls}} = {u_Ls:.4f}$ m/s")
            st.write(f"**Superficial Gas Velocity:** $U_{{Gs}} = {u_Gs:.4f}$ m/s")
            st.write(f"**Martinelli Parameter (X):** {X:.4f}")
            st.write(f"**Froude Number-based Parameter (F):** {F:.4f}")
            st.write(f"**K Parameter:** {K:.4f}")
        st.markdown("---")

        fr1, fr2, fr3 = st.columns([1, 5, 1])
        with fr2:
        # Plot the map and get the flow regime
        # The streamlit object 'st' is passed to the plot_map function
            regime = fr.plot_map(X, F, K, st)
        
        

        st.markdown("#### Flow Regime Classification")
        st.info(f"**The calculated point falls within the '{regime}' regime.**")
        
        if 'results' not in st.session_state:
            st.session_state.results = {}
        st.session_state.results['flow_regime_td_result'] = regime

# --- Page: Calculation Steps ---
elif page == "Calculation Results":
    st.subheader("3. Calculation Results")

    if st.session_state.calculation_results:
        results = st.session_state.calculation_results
        inputs = st.session_state.inputs

        # Define unit labels for SI system
        len_unit = "m"
        dens_unit = "kg/m³"
        vel_unit = "m/s"
        visc_unit = "Pa·s"
        momentum_unit = "Pa"
        micron_unit_label = "µm"
        mass_flow_unit = "kg/s"
        vol_flow_unit = "m³/s" # New unit for Streamlit display
        pressure_unit = "psig"
        in_unit = "in"

        # Display inputs used for calculation (original SI values)
        st.markdown("<h3 style='font-size: 24px;'><b>Inputs Used for Calculation (SI Units)</b></h3>", unsafe_allow_html=True)
        cl1,cl2,cl3 = st.columns(3)
        with cl1:
            #st.write(f"Pipe Inside Diameter (D): {inputs['D_input']:.4f} {len_unit}")
            st.write(f"Liquid Density (ρl): {inputs['rho_l_input']:.2f} {dens_unit}")
            st.write(f"Liquid Viscosity (μl): {inputs['mu_l_input']:.8f} {visc_unit}")
            sigma_display_val = inputs['sigma_custom'] # Use sigma_custom for display
            st.write(f"Liquid Surface Tension (σ): {sigma_display_val:.3f} N/m")
            st.write(f"Gas Density (ρg): {inputs['rho_g_input']:.5f} {dens_unit}")
            st.write(f"Gas Viscosity (μg): {inputs['mu_g_input']:.9f} {visc_unit}")
            
        rho_L = st.session_state.inputs['rho_l_input']
        mu_L = st.session_state.inputs['mu_l_input']
        rho_G = st.session_state.inputs['rho_g_input']
        mu_G = st.session_state.inputs['mu_g_input']
        Q_liquid_mass_flow_rate = st.session_state.inputs['Q_liquid_mass_flow_rate_input']
        Q_gas_mass_flow_rate = st.session_state.inputs['Q_gas_mass_flow_rate_input']
        D = st.session_state.inputs['D_input']  
        u_Gs = (Q_gas_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0
        u_ls = (Q_liquid_mass_flow_rate / rho_G) / (0.25 * math.pi * D**2) if rho_G > 0 and D > 0 else 0

        Q_gas_mass_lb_s = Q_gas_mass_flow_rate * 2.20462
        Q_liquid_mass_lb_s = Q_liquid_mass_flow_rate * 2.20462
        Q_gas_vol_ft3_s = Q_gas_mass_lb_s / (rho_G * 0.062428)  # Convert gas density from kg/m³ to lb/ft³
        Q_liquid_vol_ft3_s = Q_liquid_mass_lb_s / (rho_L * 0.062428)  # Convert liquid density from kg/m³ to lb/ft³
        
        total_vol_flow = Q_gas_vol_ft3_s + Q_liquid_vol_ft3_s
        # Volume fractions
        alpha_g = Q_gas_vol_ft3_s / total_vol_flow
        alpha_l = 1 - alpha_g

        rho_mix_fps = alpha_g * (rho_G *0.062428) + alpha_l * (rho_L *0.062428)  # Convert to lb/ft³
        # Convert velocities to ft/s    
        area_ft2 = math.pi * ((D * M_TO_FT) ** 2) / 4.0
        V_mix_fps = total_vol_flow / area_ft2 if area_ft2 > 0 else 0.0

        rho_v_squared_fps = rho_mix_fps * V_mix_fps**2

        # Display selected surface tension in SI units
        with cl2:
            st.write(f"Liquid Mass Flow Rate (Q_liquid): {inputs['Q_liquid_mass_flow_rate_input']:.2f} {mass_flow_unit}")
            st.write(f"Gas Mass Flow Rate (Q_gas): {inputs['Q_gas_mass_flow_rate_input']:.2f} {mass_flow_unit}")
            st.write(f"Pipe Inside Diameter (D): {inputs['D_input']:.4f} {len_unit}")

            st.write(f"**Superficial Gas Velocity(Calc) (Vg):** {u_Gs:.2f} {vel_unit}")
            st.write(f"**Superficial Liquid Velocity(Calc) (Vl):** {u_ls:.2f} {vel_unit}")
            #st.write(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} {pressure_unit}")
           #st.write(f"Gas Velocity in Feed Pipe (Vg): {u_Gs:.2f} {vel_unit}")
        with cl3:
            st.write(f"Separator Type: {inputs['separator_type']}")
            st.write(f"Selected Inlet Device: {inputs['inlet_device']}")
            if inputs['separator_type'] == "Horizontal":
                st.write(f"Gas Space Height (hg): {inputs['h_g_input']:.3f} {len_unit}")
                st.write(f"Effective Separation Length (Le): {inputs['L_e_input']:.3f} {len_unit}")
            else: # Vertical
                st.write(f"Separator Diameter: {inputs['D_separator_input']:.3f} {len_unit}")
                st.write(f"Gas Gravity Section Height (Le): {inputs['L_e_input']:.3f} {len_unit}")
            #st.write(f"Liquid Mass Flow Rate: {inputs['Q_liquid_mass_flow_rate_input']:.2f} {mass_flow_unit}") # New input
            #st.write(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} {pressure_unit}")
            st.write(f"Length from Inlet Device to Mist Extractor (L_to_ME): {inputs['L_to_ME_input']:.3f} {len_unit}")
            st.write(f"Perforated Plate Used: {'Yes' if inputs['perforated_plate_option'] else 'No'}")
            st.write(f"Mist Extractor Type: {inputs['mist_extractor_type']}")
            if inputs['mist_extractor_type'] == "Mesh Pad":
                st.write(f"  Mesh Pad Type: {inputs['mesh_pad_type']}")
                st.write(f"  Mesh Pad Thickness: {inputs['mesh_pad_thickness_in']:.2f} {in_unit}")
            elif inputs['mist_extractor_type'] == "Vane-Type":
                st.write(f"  Vane Type: {inputs['vane_type']}")
                st.write(f"  Flow Direction: {inputs['vane_flow_direction']}")
                st.write(f"  Number of Bends: {inputs['vane_num_bends']}")
                st.write(f"  Vane Spacing: {inputs['vane_spacing_in']:.2f} {in_unit}")
                st.write(f"  Bend Angle: {inputs['vane_bend_angle_deg']:.1f} deg")
            elif inputs['mist_extractor_type'] == "Cyclonic":
                st.write(f"  Cyclone Type: {inputs['cyclone_type']}")
                st.write(f"  Cyclone Diameter: {inputs['cyclone_diameter_in']:.2f} {in_unit}")
                st.write(f"  Cyclone Length: {inputs['cyclone_length_in']:.2f} {in_unit}")
                st.write(f"  Inlet Swirl Angle: {inputs['cyclone_swirl_angle_deg']:.1f} deg")
        
        st.markdown("---")

        # Step 1: Calculate Superficial Gas Reynolds Number (Re_g)
        st.markdown("#### Step 1: Calculate Superficial Gas Reynolds Number ($Re_g$)")
        D_pipe_fps = to_fps(inputs['D_input'], "length")
        V_g_input_fps = u_Gs * M_TO_FT #to_fps(inputs['V_g_input'], "velocity")
        rho_g_fps = to_fps(inputs['rho_g_input'], "density")
        mu_g_fps = to_fps(inputs['mu_g_input'], "viscosity")

        st.write(f"Equation: $Re_g = \\frac{{D \\cdot V_g \\cdot \\rho_g}}{{\\mu_g}}$")
        st.write(f"Calculation (FPS): $Re_g = \\frac{{{D_pipe_fps:.2f} \\text{{ ft}} \\cdot {V_g_input_fps:.2f} \\text{{ ft/sec}} \\cdot {rho_g_fps:.4f} \\text{{ lb/ft}}^3}}{{{mu_g_fps:.8f} \\text{{ lb/ft-sec}}}} = {results['Re_g']:.2f}$")
        st.success(f"**Result:** Superficial Gas Reynolds Number ($Re_g$) = **{results['Re_g']:.2f}** (dimensionless)")

        st.markdown("---")

        # Step 2: Calculate Volume Median Diameter ($d_{v50}$) without inlet device effect
        st.markdown("#### Step 2: Calculate Initial Volume Median Diameter ($d_{v50}$) (*Kataoka et al., 1983*)")
        rho_l_fps = to_fps(inputs['rho_l_input'], "density")
        mu_l_fps = to_fps(inputs['mu_l_input'], "viscosity")

        dv50_original_display = from_fps(results['dv50_original_fps'], "length")

        # Updated LaTeX formula for display
        st.write(f"Equation: $d_{{v50}} = 0.01 \\left(\\frac{{\\sigma}}{{\\rho_g V_g^2}}\\right) Re_g^{{2/3}} \\left(\\frac{{\\rho_g}}{{\\rho_l}}\\right)^{{-1/3}} \\left(\\frac{{\\mu_g}}{{\\mu_l}}\\right)^{{2/3}}$")
        st.write(f"Calculation (FPS): $d_{{v50}} = 0.01 \\left(\\frac{{{inputs['sigma_fps']:.4f}}}{{{rho_g_fps:.4f} \\cdot ({V_g_input_fps:.2f})^2}}\\right) ({results['Re_g']:.2f})^{{2/3}} \\left(\\frac{{{rho_g_fps:.4f}}}{{{rho_l_fps:.2f}}}\\right)^{{-0.333}} \\left(\\frac{{{mu_g_fps:.8f}}}{{{mu_l_fps:.7f}}}\\right)^{{0.667}}$")
        st.success(f"**Result:** Initial Volume Median Diameter ($d_{{v50}}$) = **{results['dv50_original_fps'] * FT_TO_MICRON:.2f} {micron_unit_label}** ({dv50_original_display:.6f} {len_unit})")

        st.markdown("---")


        # Step 3: Determine Inlet Momentum (rho_g V_g^2)
        st.markdown("#### Step 3: Calculate Inlet Momentum ($\\rho_m V_m^2$)")
        rho_v_squared_display = from_fps(results['rho_v_squared_fps'], "momentum")
        st.write(f"Equation: $\\rho_m V_m^2 = \\rho_m \\cdot V_m^2$")
        st.write(f"Calculation (FPS): $\\rho_m V_m^2 = {rho_mix_fps:.4f} \\text{{ lb/ft}}^3 \\cdot ({V_mix_fps:.2f} \\text{{ ft/sec}})^2 = {results['rho_v_squared_fps']:.2f} \\text{{ lb/ft-sec}}^2$")
        st.success(f"**Result:** Inlet Momentum ($\\rho_m V_m^2$) = **{rho_v_squared_display:.2f} {momentum_unit}**")

        st.markdown("---")

        # Step 4: Apply Inlet Device "Droplet Size Distribution Shift Factor"
        st.markdown("#### Step 4: Liquid Separation Efficiency / Droplet Size Distribution Shift Factor")
        st.write(f"Selected Inlet Device: **{inputs['inlet_device']}**")
        dv50_adjusted_display = from_fps(results['dv50_adjusted_fps'], "length")
        st.write(f"For an inlet momentum of {rho_v_squared_display:.2f} {momentum_unit} and a '{inputs['inlet_device']}' device, the estimated shift factor is **{results['shift_factor']:.3f}**.")
        st.write(f"Equation: $d_{{v50, adjusted}} = d_{{v50, original}} \\cdot \\text{{Shift Factor}}$")
        st.write(f"Calculation (FPS): $d_{{v50, adjusted}} = {results['dv50_original_fps']:.6f} \\text{{ ft}} \\cdot {results['shift_factor']:.3f} = {results['dv50_adjusted_fps']:.6f} \\text{{ ft}}$")
        st.success(f"**Result:** Adjusted Volume Median Diameter ($d_{{v50}}$) = **{results['dv50_adjusted_fps'] * FT_TO_MICRON:.2f} {micron_unit_label}** ({dv50_adjusted_display:.6f} {len_unit})")

        st.markdown("---")

        # Step 5: Calculate parameters for Upper-Limit Log Normal Distribution
        st.markdown("#### Step 5: Calculate Parameters for Upper-Limit Log Normal Distribution ")
        #st.caption("*Kataoka et al., 1983 and Simmons and Hanratty (2001)*")
        d_max_original_display = from_fps(results['d_max_original_fps'], "length")
        d_max_adjusted_display = from_fps(results['d_max_adjusted_fps'], "length")
        st.write(f"Using typical values: $a = {A_DISTRIBUTION}$ and $\\delta = {DELTA_DISTRIBUTION}$. *(Pan and Hanratty)*")
        st.write(f"For **Original** $d_{{v50}}$:")
        st.write(f"Equation: $d_{{max, original}} = a \\cdot d_{{v50, original}}$")
        st.write(f"Calculation (FPS): $d_{{max, original}} = {A_DISTRIBUTION} \\cdot {results['dv50_original_fps']:.6f} \\text{{ ft}} = {results['d_max_original_fps']:.6f} \\text{{ ft}}$")
        st.success(f"**Result:** Maximum Droplet Size (Original $d_{{max}}$) = **{results['d_max_original_fps'] * FT_TO_MICRON:.2f} {micron_unit_label}** ({d_max_original_display:.6f} {len_unit})")
        st.write(f"For **Adjusted** $d_{{v50}}$:")
        st.write(f"Equation: $d_{{max, adjusted}} = a \\cdot d_{{v50, adjusted}}$")
        st.write(f"Calculation (FPS): $d_{{max, adjusted}} = {A_DISTRIBUTION} \\cdot {results['dv50_adjusted_fps']:.6f} \\text{{ ft}} = {results['d_max_adjusted_fps']:.6f} \\text{{ ft}}$")
        st.success(f"**Result:** Maximum Droplet Size (Adjusted $d_{{max}}$) = **{results['d_max_adjusted_fps'] * FT_TO_MICRON:.2f} {micron_unit_label}** ({d_max_adjusted_display:.6f} {len_unit})")

        st.markdown("---")

        # Step 6: Entrainment Fraction (E) Calculation
        st.markdown("#### Step 6: Calculate Entrainment Fraction (E)")
        st.write(f"The entrainment fraction is calculated using the following correlation *(Pan and Hanratty (2002)) :*")
        st.latex(r"""\frac{\left( \frac{E}{E_M} \right)}{1 - \left( \frac{E}{E_M} \right)}
                = A_2 \left( \frac{D U_G^3 \rho_L^{0.5} \rho_G^{0.5}}{\sigma} \right)
                \left( \frac{\rho_G^{1-m} \mu_G^{m}}{d^{1+m} g \rho_L} \right)^{\frac{1}{2-m}}
                """)
        st.latex(r"""\left( \frac{\rho_G U_G^2 d_{32}}{\sigma} \right)\left( \frac{d_{32}}{D} \right) = 0.0091""")

        st.write(f"Gas Velocity ($U_{{G}}$): {u_Gs:.2f} {vel_unit}")
        hrly = Q_liquid_mass_flow_rate * 3600  # Convert kg/s to kg/hr
        st.write(f"Liquid Loading ($W_{{L}}$): {Q_liquid_mass_flow_rate:.2f} {mass_flow_unit}; ({hrly:.2f} kg/hr)")
        st.success(f"**Result:** Entrainment Fraction (E) = **{results['E_fraction']:.4f} ; ({results['E_fraction'] * 100:.2f}%)**")
        hrly = results['Q_entrained_total_mass_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
        st.success(f"**Result:** Total Entrained Liquid Mass Flow Rate = **{results['Q_entrained_total_mass_flow_rate_si']:.4f} {mass_flow_unit} ; ({hrly:.2f} kg/hr)**")
        hrly = results['Q_entrained_total_volume_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
        st.success(f"**Result:** Total Entrained Liquid Volume Flow Rate = **{results['Q_entrained_total_volume_flow_rate_si']:.6f} {vol_flow_unit} ; ({hrly:.6f} m³/hr)**") # New total volume flow
        st.markdown("---")

        # Step 7: Calculate F-factor and Effective Gas Velocity in Separator
        st.markdown("#### Step 7: Calculate F-factor and Effective Gas Velocity in Separator")
        D_separator_fps = to_fps(inputs['D_separator_input'], "length")
        L_to_ME_fps = to_fps(inputs['L_to_ME_input'], 'length')
        L_over_Di = L_to_ME_fps / D_separator_fps
        st.write(f"L/Di Ratio (Length from Inlet Device to Mist Extractor / Vessel Diameter): {L_to_ME_fps:.2f} ft / {D_separator_fps:.2f} ft = {L_over_Di:.2f}")
        st.write(f"Inlet Device: {inputs['inlet_device']}")
        st.write(f"Perforated Plate Used: {'Yes' if inputs['perforated_plate_option'] else 'No'}")
        st.write(f"Calculated F-factor : {results['F_factor']:.3f}")

        
        V_g_effective_separator_display = from_fps(results['V_g_effective_separator_fps'], 'velocity')
        if inputs['separator_type'] == "Vertical":
            D_separator_fps = to_fps(inputs['D_separator_input'], 'length')
            #A_pipe_fps = np.pi * (D_pipe_fps / 2)**2
            A_separator_gas_vertical_fps = np.pi * (D_separator_fps / 2)**2
            V_g_superficial_separator_fps = (Q_gas_vol_ft3_s) / A_separator_gas_vertical_fps
            st.write(f"Superficial Gas Velocity in Vertical Separator: {from_fps(V_g_superficial_separator_fps, 'velocity'):.2f} {vel_unit}")
            st.write(f"Equation: $V_{{g,effective}} = V_{{g,superficial}} * F$")
            st.write(f"Calculation (FPS): $V_{{g,effective}} = {V_g_superficial_separator_fps:.2f} \\text{{ ft/sec}} * {results['F_factor']:.3f} = {results['V_g_effective_separator_fps']:.2f} \\text{{ ft/sec}}$")
        else: # Horizontal
            st.write(f"Equation: $V_{{g,effective}} = V_{{g,input}} * F$ - Factor")
            st.write(f"Calculation (FPS): $V_{{g,effective}} = {results['V_g_superficial_separator_fps']:.2f} \\text{{ ft/sec}} * {results['F_factor']:.3f} = {results['V_g_effective_separator_fps']:.2f} \\text{{ ft/sec}}$")
        st.success(f"**Result:** Effective Gas Velocity in Separator ($V_{{g,effective}}$) = **{V_g_effective_separator_display:.2f} {vel_unit}**")
        st.markdown("---")

        # Step 8: Gas Gravity Separation Section Efficiency
        st.markdown("#### Step 8: Gas Gravity Separation Section Efficiency")
        st.write(f"Separator Type: **{inputs['separator_type']}**")
        if inputs['separator_type'] == "Horizontal":
            st.write(f"Gas Space Height (h_g): {inputs['h_g_input']:.3f} {len_unit}")
            st.write(f"Effective Separation Length (L_e): {inputs['L_e_input']:.3f} {len_unit}")
            st.write("For each droplet size, the separation efficiency is calculated based on its terminal velocity and the available settling time/distance.")
        else: # Vertical
            st.write(f"Separator Diameter: {inputs['D_separator_input']:.3f} {len_unit}")
            st.write(f"Gas Gravity Section Height (L_e): {inputs['L_e_input']:.3f} {len_unit}")
            st.write("For a vertical separator, a droplet is separated if its terminal settling velocity is greater than the effective upward gas velocity.")

        st.success(f"**Result:** Overall Separation Efficiency of Gas Gravity Section = **{results['gravity_separation_efficiency']:.2%}**")
        if st.session_state.plot_data_after_gravity:
            hrly = st.session_state.plot_data_after_gravity['total_entrained_mass_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Result:** Total Entrained Liquid Mass Flow Rate After Gravity Settling = **{st.session_state.plot_data_after_gravity['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit}; ({hrly:.4f} kg/hr)**")
            hrly = st.session_state.plot_data_after_gravity['total_entrained_volume_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Result:** Total Entrained Liquid Volume Flow Rate After Gravity Settling = **{st.session_state.plot_data_after_gravity['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit}; ({hrly:.4f} m³/hr)**")
        else:
            st.warning("Gravity settling results not available. Please check inputs and previous steps.")

        # Display detailed table for gravity separation
        if st.session_state.plot_data_after_gravity and st.session_state.plot_data_after_gravity['gravity_details_table_data']:
            st.markdown("##### Detailed Separation Performance in Gas Gravity Section")
            gravity_table_df = pd.DataFrame(st.session_state.plot_data_after_gravity['gravity_details_table_data'])

            # Format columns for display
            st.dataframe(gravity_table_df.style.format({
                "dp_microns": "{:.2f}",
                "Vt_ftps": "{:.2f}",
                "Cd": "{:.3f}",
                "Re_p": "{:.0f}", # Scientific notation for Reynolds number
                "Flow Regime": "{}",
                "Time Settle (s)": "{:.2f}",
                "h_max_settle (ft)": "{:.6f}",
                "Edp": "{:.2%}" # Percentage for efficiency
            }))
        else:
            st.info("Detailed separation performance for gravity section not available.")


        st.markdown("---")

        # Step 9: Mist Extractor Performance
        st.markdown("#### Step 9: Mist Extractor Performance")
        st.write(f"Mist Extractor Type: **{inputs['mist_extractor_type']}**")
        #st.write(f"Operating Pressure: {inputs['pressure_psig_input']:.1f} {pressure_unit}")
        #st.write(f"K-Deration Factor (from Table 3): {results['k_deration_factor']:.3f}")

        if inputs['mist_extractor_type'] == "Mesh Pad":
            mesh_pad_params_fps = results['mesh_pad_params']
            st.write(f"  Mesh Pad Type: {inputs['mesh_pad_type']}")
            st.write(f"  Mesh Pad Thickness: {inputs['mesh_pad_thickness_in']:.2f} {in_unit}")
            st.write(f"  Wire Diameter: {mesh_pad_params_fps['wire_diameter_in']:.3f} {in_unit}")
            st.write(f"  Specific Surface Area: {mesh_pad_params_fps['specific_surface_area_ft2_ft3']:.1f} ft²/ft³")
            #st.write(f"  Base K_s: {mesh_pad_params_fps['Ks_ft_sec']:.2f} ft/sec")
            #st.write(f"  Liquid Load Capacity: {mesh_pad_params_fps['liquid_load_gal_min_ft2']:.2f} gal/min/ft²")
            st.write("  Efficiency calculated using Stokes' number, single-wire efficiency, and mesh-pad removal efficiency.")
            # Display detailed table for mesh pad separation
            if st.session_state.calculation_results and st.session_state.calculation_results['mist_extractor_details_table_data']:
                st.markdown("##### Detailed Separation Performance in Mesh Pad")
                me_table_df = pd.DataFrame(st.session_state.calculation_results['mist_extractor_details_table_data'])
                st.dataframe(me_table_df.style.format({
                    "dp_microns": "{:.2f}",
                    "Stokes Number": "{:.2f}",
                    "Ew": "{:.2%}",
                    "Epad": "{:.2%}"
                }))
            else:
                st.info("Detailed droplet separation data for mesh pad not available.")
                
        elif inputs['mist_extractor_type'] == "Vane-Type":
            vane_type_params_fps = results['vane_type_params']
            st.write(f"  Vane Type: {inputs['vane_type']}")
            st.write(f"  Flow Direction: {inputs['vane_flow_direction']}")
            st.write(f"  Number of Bends: {inputs['vane_num_bends']}")
            st.write(f"  Vane Spacing: {inputs['vane_spacing_in']:.2f} {in_unit}")
            st.write(f"  Bend Angle: {inputs['vane_bend_angle_deg']:.1f} deg")
            st.write(f"  Base K_s (Upflow): {vane_type_params_fps['Ks_ft_sec_upflow']:.2f} ft/sec")
            st.write(f"  Base K_s (Horizontal): {vane_type_params_fps['Ks_ft_sec_horizontal']:.2f} ft/sec")
            st.write(f"  Liquid Load Capacity: {vane_type_params_fps['liquid_load_gal_min_ft2']:.2f} gal/min/ft²")
            st.write("  Efficiency calculated using Eq. 15.")
                        # Display detailed table for vane type separation
            if st.session_state.calculation_results and st.session_state.calculation_results['mist_extractor_details_table_data']:
                st.markdown("##### Detailed Separation Performance in Vane Type")
                me_table_df = pd.DataFrame(st.session_state.calculation_results['mist_extractor_details_table_data'])
                st.dataframe(me_table_df.style.format({
                    "dp_microns": "{:.2f}",
                    "Evane": "{:.2%}"
                }))
            else:
                st.info("Detailed droplet separation data for mesh pad not available.")

        elif inputs['mist_extractor_type'] == "Cyclonic":
            cyclone_type_params_fps = results['cyclone_type_params']
            st.write(f"  Cyclone Type: {inputs['cyclone_type']}")
            st.write(f"  Cyclone Diameter: {inputs['cyclone_diameter_in']:.2f} {in_unit}")
            st.write(f"  Cyclone Length: {inputs['cyclone_length_in']:.2f} {in_unit}")
            st.write(f"  Inlet Swirl Angle: {inputs['cyclone_swirl_angle_deg']:.1f} deg")
            st.write(f"  Base K_s: {cyclone_type_params_fps['Ks_ft_sec_bundle_face_area']:.2f} ft/sec")
            st.write(f"  Liquid Load Capacity: {cyclone_type_params_fps['liquid_load_gal_min_ft2_bundle_face_area']:.2f} gal/min/ft²")
            st.write("  Efficiency calculated using Eq. 16.")
                        # Display detailed table for cyclonic separation
            if st.session_state.calculation_results and st.session_state.calculation_results['mist_extractor_details_table_data']:
                st.markdown("##### Detailed Droplet Separation Performance in Cyclonic")
                me_table_df = pd.DataFrame(st.session_state.calculation_results['mist_extractor_details_table_data'])
                st.dataframe(me_table_df.style.format({
                    "dp_microns": "{:.2f}",
                    "Stk": "{:.1f}",
                    "E_cycl": "{:.2%}"
                }))
            else:
                st.info("Detailed droplet separation data for mesh pad not available.")

        st.success(f"**Result:** Overall Separation Efficiency of Mist Extractor = **{results['mist_extractor_separation_efficiency']:.2%}**")
        if st.session_state.plot_data_after_mist_extractor:
            hrly = st.session_state.plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Result:** Total Entrained Liquid Mass Flow Rate After Mist Extractor = **{st.session_state.plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit}; ({hrly:.4f} kg/hr)**")
            hrly = st.session_state.plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Result:** Total Entrained Liquid Volume Flow Rate After Mist Extractor = **{st.session_state.plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit}; ({hrly:.6f} m³/hr)**")
        else:
            st.warning("Mist extractor results not available. Please check inputs and previous steps.")

        st.markdown("---")
        st.subheader("Final Carry-Over from Separator Outlet")
        if st.session_state.plot_data_after_mist_extractor:
            hrly = st.session_state.plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Total Carry-Over Mass Flow Rate:** **{st.session_state.plot_data_after_mist_extractor['total_entrained_mass_flow_rate_si']:.4f} {mass_flow_unit}; ({hrly:.4f} kg/hr)**")
            hrly = st.session_state.plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si'] * 3600  # Convert kg/s to kg/hr
            st.success(f"**Total Carry-Over Volume Flow Rate:** **{st.session_state.plot_data_after_mist_extractor['total_entrained_volume_flow_rate_si']:.6f} {vol_flow_unit}; ({hrly:.6f} m³/hr)**")
        else:
            st.warning("Final carry-over results not available. Please ensure all previous steps are calculated.")

    else:
        st.warning("Please go to the 'Input Parameters' page and modify inputs to trigger calculations.")

# --- Page: Droplet Distribution Results ---
elif page == "Carry Over Plots":
    st.subheader("4. Droplet Size Distribution Plot")

    # Moved the input for num_points_distribution here
    st.markdown("<h3 style='font-size: 24px;'><b>Distribution Plot Settings</b></h3>", unsafe_allow_html=True)
    #st.subheader("Distribution Plot Settings")
    tbl1, tbl2, tbl3, tbl4 = st.columns(4)
    with tbl1:
        st.session_state.inputs['num_points_distribution'] = st.number_input(
            "Number of Points",
            min_value=10,
            max_value=100,
            value=st.session_state.inputs['num_points_distribution'],
            step=5,
            key='num_points_distribution_input',
            help="Adjust the number of data points used to generate the droplet size distribution curve and table (10-100)."
        )

    st.markdown("---") # Add a separator after the input

    # Recalculate plot_data specifically on this page after num_points_distribution is updated
    # Ensure calculation_results are available before proceeding
    if st.session_state.calculation_results:
        try:
            results = st.session_state.calculation_results
            inputs = st.session_state.inputs
            num_points = inputs['num_points_distribution']

            # Ensure required inputs for distribution generation are available
            if 'Q_liquid_mass_flow_rate_input' in inputs and \
               inputs['Q_liquid_mass_flow_rate_input'] is not None and \
               'rho_l_input' in inputs and \
               inputs['rho_l_input'] is not None:

                # Generate initial distribution (after inlet device, before gravity settling)
                st.session_state.plot_data_original = _generate_initial_distribution_data(
                    results['dv50_original_fps'],
                    results['d_max_original_fps'],
                    num_points,
                    results['E_fraction'],
                    inputs['Q_liquid_mass_flow_rate_input'],
                    inputs['rho_l_input']
                )

                st.session_state.plot_data_adjusted = _generate_initial_distribution_data(
                    results['dv50_adjusted_fps'],
                    results.get('d_max_adjusted_fps', results['d_max_original_fps']), # Use original if adjusted not present
                    num_points,
                    results['E_fraction'],
                    inputs['Q_liquid_mass_flow_rate_input'],
                    inputs['rho_l_input']
                )

                # Calculate and apply gravity settling
                if inputs['separator_type'] == "Horizontal":
                    st.session_state.plot_data_after_gravity = _calculate_and_apply_separation(
                        st.session_state.plot_data_adjusted, # Input is the adjusted distribution
                        _separation_stage_efficiency_func=gravity_efficiency_func_horizontal,
                        is_gravity_stage=True,
                        V_g_eff_sep_fps=results['V_g_effective_separator_fps'],
                        h_g_sep_fps=to_fps(inputs['h_g_input'], 'length'),
                        L_e_sep_fps=to_fps(inputs['L_e_input'], 'length'),
                        rho_l_fps=to_fps(inputs['rho_l_input'], 'density'),
                        rho_g_fps=to_fps(inputs['rho_g_input'], 'density'),
                        mu_g_fps=to_fps(inputs['mu_g_input'], 'viscosity'),
                        separator_type=inputs['separator_type']
                    )
                else: # Vertical
                    st.session_state.plot_data_after_gravity = _calculate_and_apply_separation(
                        st.session_state.plot_data_adjusted, # Input is the adjusted distribution
                        _separation_stage_efficiency_func=gravity_efficiency_func_vertical,
                        is_gravity_stage=True,
                        V_g_eff_sep_fps=results['V_g_effective_separator_fps'],
                        # For vertical, h_g_input is effectively L_e_input for gravity calculations
                        h_g_sep_fps=to_fps(inputs['L_e_input'], 'length'), 
                        L_e_sep_fps=to_fps(inputs['L_e_input'], 'length'), # Pass L_e_input for vertical time_settle calc
                        rho_l_fps=to_fps(inputs['rho_l_input'], 'density'),
                        rho_g_fps=to_fps(inputs['rho_g_input'], 'density'),
                        mu_g_fps=to_fps(inputs['mu_g_input'], 'viscosity'),
                        separator_type=inputs['separator_type']
                    )

                # Update the overall gravity separation efficiency in results for reporting
                if st.session_state.plot_data_after_gravity:
                    st.session_state.calculation_results['gravity_separation_efficiency'] = st.session_state.plot_data_after_gravity['overall_separation_efficiency']
                else:
                    st.session_state.calculation_results['gravity_separation_efficiency'] = 0.0

                # Calculate and apply mist extractor efficiency
                if st.session_state.plot_data_after_gravity and st.session_state.plot_data_after_gravity['dp_values_ft'].size > 0:
                    if st.session_state.inputs['mist_extractor_type'] == "Mesh Pad":
                        mesh_pad_params = MESH_PAD_PARAMETERS[st.session_state.inputs['mesh_pad_type']]
                        mesh_pad_params_with_user_thickness = mesh_pad_params.copy()
                        mesh_pad_params_with_user_thickness["thickness_in"] = st.session_state.inputs['mesh_pad_thickness_in']

                        st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                            st.session_state.plot_data_after_gravity,
                            _separation_stage_efficiency_func=mesh_pad_efficiency_func,
                            V_g_eff_sep_fps=results['V_g_effective_separator_fps'],
                            rho_l_fps=to_fps(inputs['rho_l_input'], 'density'),
                            rho_g_fps=to_fps(inputs['rho_g_input'], 'density'),
                            mu_g_fps=to_fps(inputs['mu_g_input'], 'viscosity'),
                            mesh_pad_type_params_fps=mesh_pad_params_with_user_thickness
                        )
                    elif st.session_state.inputs['mist_extractor_type'] == "Vane-Type":
                        vane_type_params = VANE_PACK_PARAMETERS[st.session_state.inputs['vane_type']]
                        vane_type_params_with_user_inputs = vane_type_params.copy()
                        vane_type_params_with_user_inputs["flow_direction"] = st.session_state.inputs['vane_flow_direction']
                        vane_type_params_with_user_inputs["number_of_bends"] = st.session_state.inputs['vane_num_bends']
                        vane_type_params_with_user_inputs["vane_spacing_in"] = st.session_state.inputs['vane_spacing_in']
                        vane_type_params_with_user_inputs["bend_angle_degree"] = st.session_state.inputs['vane_bend_angle_deg']

                        st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                            st.session_state.plot_data_after_gravity,
                            _separation_stage_efficiency_func=vane_type_efficiency_func,
                            V_g_eff_sep_fps=results['V_g_effective_separator_fps'],
                            rho_l_fps=to_fps(inputs['rho_l_input'], 'density'),
                            rho_g_fps=to_fps(inputs.get('rho_g_input', 0.0), 'density'), # Use .get with default for robustness
                            mu_g_fps=to_fps(inputs.get('mu_g_input', 0.0), 'viscosity'), # Use .get with default for robustness
                            vane_type_params_fps=vane_type_params_with_user_inputs
                        )
                    elif st.session_state.inputs['mist_extractor_type'] == "Cyclonic":
                        cyclone_type_params = CYCLONE_PARAMETERS[st.session_state.inputs['cyclone_type']]
                        cyclone_type_params_with_user_inputs = cyclone_type_params.copy()
                        cyclone_type_params_with_user_inputs["cyclone_inside_diameter_in"] = st.session_state.inputs['cyclone_diameter_in']
                        cyclone_type_params_with_user_inputs["cyclone_length_in"] = st.session_state.inputs['cyclone_length_in']
                        cyclone_type_params_with_user_inputs["inlet_swirl_angle_degree"] = st.session_state.inputs['cyclone_swirl_angle_deg']

                        st.session_state.plot_data_after_mist_extractor = _calculate_and_apply_separation(
                            st.session_state.plot_data_after_gravity,
                            _separation_stage_efficiency_func=demisting_cyclone_efficiency_func,
                            V_g_eff_sep_fps=results['V_g_effective_separator_fps'],
                            rho_l_fps=to_fps(inputs['rho_l_input'], 'density'),
                            rho_g_fps=to_fps(inputs.get('rho_g_input', 0.0), 'density'), # Use .get with default for robustness
                            mu_g_fps=to_fps(inputs.get('mu_g_input', 0.0), 'viscosity'), # Use .get with default for robustness
                            cyclone_type_params_fps=cyclone_type_params_with_user_inputs
                        )
                    else:
                        st.session_state.plot_data_after_mist_extractor = st.session_state.plot_data_after_gravity # No mist extractor selected, so no change

                    # Update the overall mist extractor separation efficiency in results for reporting
                    if st.session_state.plot_data_after_mist_extractor:
                        st.session_state.calculation_results['mist_extractor_separation_efficiency'] = st.session_state.plot_data_after_mist_extractor['overall_separation_efficiency']
                        st.session_state.calculation_results['mist_extractor_details_table_data'] = st.session_state.plot_data_after_mist_extractor['mist_extractor_details_table_data']
                    else:
                        st.session_state.calculation_results['mist_extractor_separation_efficiency'] = 0.0
                        st.session_state.calculation_results['mist_extractor_details_table_data'] = []
                else:
                    st.session_state.plot_data_after_mist_extractor = None # No gravity data, so no mist extractor data either

            else:
                st.warning("Required liquid flow rate or density inputs are missing in session state. Please check 'Input Parameters' page.")
                st.session_state.plot_data_original = None
                st.session_state.plot_data_adjusted = None
                st.session_state.plot_data_after_gravity = None
                st.session_state.plot_data_after_mist_extractor = None

        except Exception as e:
            st.error(f"An error occurred during plot data calculation: {e}")
            st.session_state.plot_data_original = None
            st.session_state.plot_data_adjusted = None
            st.session_state.plot_data_after_gravity = None
            st.session_state.plot_data_after_mist_extractor = None
    else:
        #st.warning("Please go to the 'Input Parameters' page and modify inputs to trigger calculations and generate the plot data.")
        st.write("")

    if st.session_state.plot_data_original and st.session_state.plot_data_adjusted and st.session_state.plot_data_after_gravity and st.session_state.plot_data_after_mist_extractor:
        plot_data_original = st.session_state.plot_data_original
        plot_data_adjusted = st.session_state.plot_data_adjusted
        plot_data_after_gravity = st.session_state.plot_data_after_gravity
        plot_data_after_mist_extractor = st.session_state.plot_data_after_mist_extractor

        # Define unit labels for plotting
        micron_unit_label = "µm" # Always SI for this version
        mass_flow_unit = "kg/s"
        vol_flow_unit = "m³/s" # New unit for Streamlit display

        # --- Plot for Original Distribution ---

        psd1, psd2, psd3 = st.columns([1, 5, 1])
        with psd2:
            #st.markdown("<h3 style='font-size: 24px;'><b>4. Particle Size Distribution</b></h3>", unsafe_allow_html=True)
            st.markdown("<h3 style='font-size: 18px;'><b>4.1. Distribution at Inlet</b></h3>", unsafe_allow_html=True)
            #st.subheader("3.1. Distribution Before Inlet Device")
            dp_values_microns_original = plot_data_original['dp_values_ft'] * FT_TO_MICRON
            fig_original, ax_original = plt.subplots(figsize=(10, 6))

            ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
            #ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
            ax_original.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
            ax_original.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
            ax_original.tick_params(axis='y', labelcolor='black')
            ax_original.axhline(y=1, color='b', linestyle='--')
            ax_original.set_ylim(0, 1.2)
            ax_original.set_xlim(0, max(dp_values_microns_original) * 1.0 if dp_values_microns_original.size > 0 else 1000)

            ax2_original = ax_original.twinx()
            ax2_original.plot(dp_values_microns_original, plot_data_original['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
            ax2_original.set_ylabel('Volume Fraction', color='black', fontsize=12)
            ax2_original.tick_params(axis='y', labelcolor='black')
            max_norm_fv_original = max(plot_data_original['volume_fraction']) if plot_data_original['volume_fraction'].size > 0 else 0.1
            ax2_original.set_ylim(0, max_norm_fv_original * 1.2)

            lines_original, labels_original = ax_original.get_legend_handles_labels()
            lines2_original, labels2_original = ax2_original.get_legend_handles_labels()
            ax2_original.legend(lines_original + lines2_original, labels_original + labels2_original, loc='upper left', fontsize=10)

            plt.title('Entrainment Droplet Size Distribution (Before Inlet Device)', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            st.pyplot(fig_original)
            plt.close(fig_original) # Close the plot to free memory

            # --- Plot for Adjusted Distribution ---
            st.markdown("<h3 style='font-size: 18px;'><b>4.2. Distribution after Inlet Device</b></h3>", unsafe_allow_html=True)
            #st.subheader("4.2. Distribution after Inlet Device")
            dp_values_microns_adjusted = plot_data_adjusted['dp_values_ft'] * FT_TO_MICRON
            fig_adjusted, ax_adjusted = plt.subplots(figsize=(8, 6))

            ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
            #ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
            ax_adjusted.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
            ax_adjusted.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
            ax_adjusted.tick_params(axis='y', labelcolor='black')
            ax_adjusted.axhline(y=1, color='b', linestyle='--')
            ax_adjusted.set_ylim(0, 1.2)
            ax_adjusted.set_xlim(0, max(dp_values_microns_adjusted) * 1.0 if dp_values_microns_adjusted.size > 0 else 1000)

            ax2_adjusted = ax_adjusted.twinx()
            ax2_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
            ax2_adjusted.set_ylabel('Volume Fraction', color='black', fontsize=12)
            ax2_adjusted.tick_params(axis='y', labelcolor='black')
            max_norm_fv_adjusted = max(plot_data_adjusted['volume_fraction']) if plot_data_adjusted['volume_fraction'].size > 0 else 0.1
            ax2_adjusted.set_ylim(0, max_norm_fv_adjusted * 1.2)

            lines_adjusted, labels_adjusted = ax_adjusted.get_legend_handles_labels()
            lines2_adjusted, labels2_adjusted = ax2_adjusted.get_legend_handles_labels()
            ax2_adjusted.legend(lines_adjusted + lines2_adjusted, labels_adjusted + labels2_adjusted, loc='upper left', fontsize=10)

            plt.title('Entrainment Droplet Size Distribution (after Inlet Device)', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            st.pyplot(fig_adjusted)
            plt.close(fig_adjusted) # Close the plot to free memory

            # --- Plot for After Gravity Settling ---
            st.markdown("<h3 style='font-size: 18px;'><b>4.3. Distribution after Gravity Settling Section</b></h3>", unsafe_allow_html=True)
            #st.subheader("3.3. Distribution After Gas Gravity Settling")
            dp_values_microns_after_gravity = plot_data_after_gravity['dp_values_ft'] * FT_TO_MICRON
            fig_after_gravity, ax_after_gravity = plt.subplots(figsize=(8, 6))

            ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
            #ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
            ax_after_gravity.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
            ax_after_gravity.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
            ax_after_gravity.tick_params(axis='y', labelcolor='black')
            ax_after_gravity.axhline(y=1, color='b', linestyle='--')
            ax_after_gravity.set_ylim(0, 1.2)
            ax_after_gravity.set_xlim(0, max(dp_values_microns_after_gravity) * 1.0 if dp_values_microns_after_gravity.size > 0 else 1000)

            ax2_after_gravity = ax_after_gravity.twinx()
            ax2_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
            ax2_after_gravity.set_ylabel('Volume Fraction', color='black', fontsize=12)
            ax2_after_gravity.tick_params(axis='y', labelcolor='black')
            max_norm_fv_after_gravity = max(plot_data_after_gravity['volume_fraction']) if plot_data_after_gravity['volume_fraction'].size > 0 else 0.1
            ax2_after_gravity.set_ylim(0, max_norm_fv_after_gravity * 1.2)

            lines_after_gravity, labels_after_gravity = ax_after_gravity.get_legend_handles_labels()
            lines2_after_gravity, labels2_after_gravity = ax2_after_gravity.get_legend_handles_labels()
            ax2_after_gravity.legend(lines_after_gravity + lines2_after_gravity, labels_after_gravity + labels2_after_gravity, loc='upper left', fontsize=10)

            plt.title('Entrainment Droplet Size Distribution (after Gravity Settling)', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            st.pyplot(fig_after_gravity)
            plt.close(fig_after_gravity) # Close the plot to free memory

            # --- Plot for After Mist Extractor ---
            st.markdown("<h3 style='font-size: 18px;'><b>4.4. Distribution after Mist Extractor</b></h3>", unsafe_allow_html=True)
            #st.subheader("3.4. Distribution After Mist Extractor")
            dp_values_microns_after_me = plot_data_after_mist_extractor['dp_values_ft'] * FT_TO_MICRON
            fig_after_me, ax_after_me = plt.subplots(figsize=(8, 6))

            ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
            #ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
            ax_after_me.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
            ax_after_me.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
            ax_after_me.tick_params(axis='y', labelcolor='black')
            ax_after_me.axhline(y=1, color='b', linestyle='--')
            ax_after_me.set_ylim(0, 1.2)
            ax_after_me.set_xlim(0, max(dp_values_microns_after_me) * 0.5 if dp_values_microns_after_me.size > 0 else 1000)

            ax2_after_me = ax_after_me.twinx()
            ax2_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
            ax2_after_me.set_ylabel('Volume Fraction', color='black', fontsize=12)
            ax2_after_me.tick_params(axis='y', labelcolor='black')
            max_norm_fv_after_me = max(plot_data_after_mist_extractor['volume_fraction']) if plot_data_after_mist_extractor['volume_fraction'].size > 0 else 0.1
            ax2_after_me.set_ylim(0, max_norm_fv_after_me * 1.2)

            lines_after_me, labels_after_me = ax_after_me.get_legend_handles_labels()
            lines2_after_me, labels2_after_me = ax2_after_me.get_legend_handles_labels()
            ax2_after_me.legend(lines_after_me + lines2_after_me, labels_after_me + labels2_after_me, loc='upper left', fontsize=10)

            plt.title('Entrainment Droplet Size Distribution (after Mist Extractor)', fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.7)
            st.pyplot(fig_after_me)
            plt.close(fig_after_me) # Close the plot to free memory

        st.markdown("---")
        # --- Volume Fraction Data Tables for Streamlit App ---
        st.subheader("5. Volume Fraction Data Tables")
        st.markdown("---")

        # Original Data Table
        st.markdown("<h3 style='font-size: 18px;'><b>5.1. Distribution at Inlet</b></h3>", unsafe_allow_html=True)
        #st.markdown("#### 4.1. Distribution Before Inlet Device")
            
        if plot_data_original['dp_values_ft'].size > 0:
            full_df_original = pd.DataFrame({
                "Droplet Size (µm)": plot_data_original['dp_values_ft'] * FT_TO_MICRON,
                "Volume Fraction": plot_data_original['volume_fraction'],
                "Cumulative Undersize": plot_data_original['cumulative_volume_undersize'],
                f"Entrained Mass Flow ({mass_flow_unit})": plot_data_original['entrained_mass_flow_rate_per_dp'],
                f"Entrained Volume Flow ({vol_flow_unit})": plot_data_original['entrained_volume_flow_rate_per_dp']
                })
                
            st.dataframe(full_df_original.style.format({
                "Droplet Size (µm)": "{:.2f}",
                "Volume Fraction": "{:.4f}",
                "Cumulative Undersize": "{:.4f}",
                f"Entrained Mass Flow ({mass_flow_unit})": "{:.6f}",
                f"Entrained Volume Flow ({vol_flow_unit})": "{:.9f}"
            })
            )
                
                
            #st.markdown(f"**Sum of Entrained Mass Flow in Table:** {np.sum(plot_data_original['entrained_mass_flow_rate_per_dp']):.6f} {mass_flow_unit}")
            hrly = st.session_state.calculation_results['Q_entrained_total_mass_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Mass Flow Rate:** {st.session_state.calculation_results['Q_entrained_total_mass_flow_rate_si']:.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            #st.markdown(f"**Sum of Entrained Volume Flow in Table:** {np.sum(plot_data_original['entrained_volume_flow_rate_per_dp']):.9f} {vol_flow_unit}")
            hrly = st.session_state.calculation_results['Q_entrained_total_volume_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Volume Flow Rate:** {st.session_state.calculation_results['Q_entrained_total_volume_flow_rate_si']:.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
            #st.info("Note: The sum of 'Entrained Flow' in the table should now precisely match the 'Total Entrained Liquid Flow Rate' from Step 6, as the volume frequency distribution is normalized and all calculated points are displayed.")
        else:
            st.info("No data available to display in the table for original distribution. Please check your input parameters.")
    
        # Adjusted Data Table
        st.markdown("<h3 style='font-size: 18px;'><b>5.2. Distribution after Inlet Device</b></h3>", unsafe_allow_html=True)
        #st.markdown("#### 4.2. Distribution After Inlet Device (Shift Factor Applied)")
        if plot_data_adjusted['dp_values_ft'].size > 0:
            full_df_adjusted = pd.DataFrame({
                "Droplet Size (µm)": plot_data_adjusted['dp_values_ft'] * FT_TO_MICRON,
                "Volume Fraction": plot_data_adjusted['volume_fraction'],
                "Cumulative Undersize": plot_data_adjusted['cumulative_volume_undersize'],
                f"Entrained Mass Flow ({mass_flow_unit})": plot_data_adjusted['entrained_mass_flow_rate_per_dp'],
                f"Entrained Volume Flow ({vol_flow_unit})": plot_data_adjusted['entrained_volume_flow_rate_per_dp']
            })
            st.dataframe(full_df_adjusted.style.format({
                "Droplet Size (µm)": "{:.2f}",
                "Volume Fraction": "{:.4f}",
                "Cumulative Undersize": "{:.4f}",
                f"Entrained Mass Flow ({mass_flow_unit})": "{:.6f}",
                f"Entrained Volume Flow ({vol_flow_unit})": "{:.9f}"
            }))
            hrly = np.sum(plot_data_adjusted['entrained_mass_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Mass Flow in Table:** {np.sum(plot_data_adjusted['entrained_mass_flow_rate_per_dp']):.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = plot_data_original['total_entrained_mass_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Mass Flow Rate (from previous stage):** {plot_data_original['total_entrained_mass_flow_rate_si']:.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = np.sum(plot_data_adjusted['entrained_volume_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Volume Flow in Table:** {np.sum(plot_data_adjusted['entrained_volume_flow_rate_per_dp']):.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
            hrly = plot_data_original['total_entrained_volume_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Volume Flow Rate (from previous stage):** {plot_data_original['total_entrained_volume_flow_rate_si']:.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
            #st.info("Note: The sum of 'Entrained Flow' in the table should now precisely match the 'Total Entrained Liquid Flow Rate' from the previous stage, as the volume frequency distribution is normalized and all calculated points are displayed.")
        else:
            st.info("No data available to display in the table for adjusted distribution. Please check your input parameters.")

        # Data Table After Gravity Settling
        st.markdown("<h3 style='font-size: 18px;'><b>5.3. Distribution after Gravity Settling Section</b></h3>", unsafe_allow_html=True)
        #st.markdown("#### 4.3. Distribution After Gas Gravity Settling")
        if plot_data_after_gravity['dp_values_ft'].size > 0:
            full_df_after_gravity = pd.DataFrame({
                "Droplet Size (µm)": plot_data_after_gravity['dp_values_ft'] * FT_TO_MICRON,
                "Volume Fraction": plot_data_after_gravity['volume_fraction'],
                "Cumulative Undersize": plot_data_after_gravity['cumulative_volume_undersize'],
                f"Entrained Mass Flow ({mass_flow_unit})": plot_data_after_gravity['entrained_mass_flow_rate_per_dp'],
                f"Entrained Volume Flow ({vol_flow_unit})": plot_data_after_gravity['entrained_volume_flow_rate_per_dp']
            })
            st.dataframe(full_df_after_gravity.style.format({
                "Droplet Size (µm)": "{:.2f}",
                "Volume Fraction": "{:.4f}",
                "Cumulative Undersize": "{:.4f}",
                f"Entrained Mass Flow ({mass_flow_unit})": "{:.6f}",
                f"Entrained Volume Flow ({vol_flow_unit})": "{:.9f}"
            }))
            hrly = np.sum(plot_data_after_gravity['entrained_mass_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Mass Flow in Table:** {np.sum(plot_data_after_gravity['entrained_mass_flow_rate_per_dp']):.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = plot_data_adjusted['total_entrained_mass_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Mass Flow Rate (from previous stage):** {plot_data_adjusted['total_entrained_mass_flow_rate_si']:.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = np.sum(plot_data_after_gravity['entrained_volume_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Volume Flow in Table:** {np.sum(plot_data_after_gravity['entrained_volume_flow_rate_per_dp']):.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
            hrly = plot_data_adjusted['total_entrained_volume_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Volume Flow Rate (from previous stage):** {plot_data_adjusted['total_entrained_volume_flow_rate_si']:.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
           # st.info("Note: The sum of 'Entrained Flow' in the table should now precisely match the 'Total Entrained Liquid Flow Rate' from the previous stage, as the volume frequency distribution is normalized and all calculated points are displayed.")
        else:
            st.info("No data available to display in the table for gravity settling. Please check your input parameters.")

        # Data Table After Mist Extractor
        st.markdown("<h3 style='font-size: 18px;'><b>5.4. Distribution after Mist Extractor</b></h3>", unsafe_allow_html=True)
        if plot_data_after_mist_extractor['dp_values_ft'].size > 0:
            full_df_after_me = pd.DataFrame({
                "Droplet Size (µm)": plot_data_after_mist_extractor['dp_values_ft'] * FT_TO_MICRON,
                "Volume Fraction": plot_data_after_mist_extractor['volume_fraction'],
                "Cumulative Undersize": plot_data_after_mist_extractor['cumulative_volume_undersize'],
                f"Entrained Mass Flow ({mass_flow_unit})": plot_data_after_mist_extractor['entrained_mass_flow_rate_per_dp'],
                f"Entrained Volume Flow ({vol_flow_unit})": plot_data_after_mist_extractor['entrained_volume_flow_rate_per_dp']
            })
            st.dataframe(full_df_after_me.style.format({
                "Droplet Size (µm)": "{:.2f}",
                "Volume Fraction": "{:.4f}",
                "Cumulative Undersize": "{:.4f}",
                f"Entrained Mass Flow ({mass_flow_unit})": "{:.6f}",
                f"Entrained Volume Flow ({vol_flow_unit})": "{:.9f}"
            }))
            hrly = np.sum(plot_data_after_mist_extractor['entrained_mass_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Mass Flow in Table:** {np.sum(plot_data_after_mist_extractor['entrained_mass_flow_rate_per_dp']):.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = plot_data_after_gravity['total_entrained_mass_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Mass Flow Rate (from previous stage):** {plot_data_after_gravity['total_entrained_mass_flow_rate_si']:.6f} {mass_flow_unit}; ({hrly:.4f} kg/hr)")
            hrly = np.sum(plot_data_after_mist_extractor['entrained_volume_flow_rate_per_dp']) * 3600.0 # Convert to hourly
            st.markdown(f"**Sum of Entrained Volume Flow in Table:** {np.sum(plot_data_after_mist_extractor['entrained_volume_flow_rate_per_dp']):.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
            hrly = plot_data_after_gravity['total_entrained_volume_flow_rate_si'] * 3600.0 # Convert to hourly
            st.markdown(f"**Total Entrained Liquid Volume Flow Rate (from previous stage):** {plot_data_after_gravity['total_entrained_volume_flow_rate_si']:.9f} {vol_flow_unit}; ({hrly:.4f} m³/hr)")
        else:
            st.info("No data available to display in the table for mist extractor. Please check your input parameters.")
        st.markdown("---")
        if st.session_state.inputs['mist_extractor_type'] == "Mesh Pad":
            st.subheader("6. Mesh Pad Efficiency & Pressure Drop ")
            st.markdown("---")
            
            # Gather parameters
            Q_gas_mass_lb_s = st.session_state.inputs['Q_gas_mass_flow_rate_input'] * 2.20462
            rho_g_fps = to_fps(st.session_state.inputs['rho_g_input'], "density")
            Q_gas_vol_ft3_s = Q_gas_mass_lb_s / rho_g_fps
            A_installed_ft2 = st.session_state.inputs['mesh_pad_A_installed_ft2']
            mesh_pad_params = MESH_PAD_PARAMETERS[st.session_state.inputs['mesh_pad_type']]
            mesh_pad_params = mesh_pad_params.copy()
            mesh_pad_params["thickness_in"] = st.session_state.inputs['mesh_pad_thickness_in']
            K_factor = mesh_pad_params["Ks_ft_sec"]
            K_dp = st.session_state.inputs['mesh_pad_K_dp']
            rho_l_fps = to_fps(st.session_state.inputs['rho_l_input'], "density")
            mu_g_fps = to_fps(st.session_state.inputs['mu_g_input'], "viscosity")
            # Plot
            fig = plot_mesh_pad_efficiency_with_pressure(
                Q_gas_vol_ft3_s, A_installed_ft2, rho_l_fps, rho_g_fps, mu_g_fps,
                mesh_pad_params, K_factor, K_dp, results
            )
            st.markdown(f"**Actual face velocity before Mesh Pad:** {results['V_face_actual_ft_s']:.2f} ft/s; ({results['V_face_actual_ft_s']*FT_TO_M:.2f} m/s)")
            st.markdown(f"**Allowable velocity before Mesh Pad:** {results['V_allow_ft_s']:.2f} ft/s; ({results['V_allow_ft_s']*FT_TO_M:.2f} m/s)")
            st.markdown(f"**Mesh Pad Estimated Pressure Drop:** {results['dp_mesh_pad_actual']:.2f} Pa; ({results['dp_mesh_pad_actual']/1000 :.2f} kPa)")
            if results['V_face_actual_ft_s'] > results['V_allow_ft_s']:
                carryover = results['carryover_percent'] / 100.0 * plot_data_after_gravity['total_entrained_volume_flow_rate_si'] * 3600.0 # Convert to hourly
                st.markdown(f"**Estimated Carryover percent :** {results['carryover_percent']:.1f} % ")
                st.markdown(f"**Estimated Carryover :** {results['carryover_percent']:.1f} % of Inlet flow rate {plot_data_after_gravity['total_entrained_volume_flow_rate_si']*3600:.3f} m3/hr = {carryover:.3f} m3/hr")
                st.warning("Flooding risk: Operating velocity exceeds allowable face velocity for Mesh Pad.")     
            else:
                st.markdown("")
                #st.markdown(f"**Operating Region within Mesh Pad Estimated Pressure Drop:** {results['dp_mesh_pad_actual']:.2f} Pa; ({results['dp_mesh_pad_actual']/1000 :.2f} kPa)")
            st.pyplot(fig)

        # Save plots to BytesIO objects for PDF embedding
        buf_original = io.BytesIO()
        fig_original = plt.figure(figsize=(10, 6)) # Recreate figure for saving
        ax_original = fig_original.add_subplot(111)
        ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
        ax_original.plot(dp_values_microns_original, plot_data_original['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_original.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
        ax_original.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_original.tick_params(axis='y', labelcolor='black')
        ax_original.set_ylim(0, 1.05)
        ax_original.set_xlim(0, max(dp_values_microns_original) * 1.1 if dp_values_microns_original.size > 0 else 1000)
        ax2_original = ax_original.twinx()
        ax2_original.plot(dp_values_microns_original, plot_data_original['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
        ax2_original.set_ylabel('Volume/Mass Fraction', color='black', fontsize=12)
        ax2_original.tick_params(axis='y', labelcolor='black')
        max_norm_fv_original = max(plot_data_original['volume_fraction']) if plot_data_original['volume_fraction'].size > 0 else 0.1
        ax2_original.set_ylim(0, max_norm_fv_original * 1.2)
        lines_original, labels_original = ax_original.get_legend_handles_labels()
        lines2_original, labels2_original = ax2_original.get_legend_handles_labels()
        ax2_original.legend(lines_original + lines2_original, labels_original + labels2_original, loc='upper left', fontsize=10)
        plt.title('Entrainment Droplet Size Distribution (Before Inlet Device)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        fig_original.savefig(buf_original, format="png", dpi=300)
        buf_original.seek(0)
        plt.close(fig_original) # Close the plot to free memory


        buf_adjusted = io.BytesIO()
        fig_adjusted = plt.figure(figsize=(10, 6)) # Recreate figure for saving
        ax_adjusted = fig_adjusted.add_subplot(111)
        ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
        ax_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_adjusted.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
        ax_adjusted.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_adjusted.tick_params(axis='y', labelcolor='black')
        ax_adjusted.set_ylim(0, 1.05)
        ax_adjusted.set_xlim(0, max(dp_values_microns_adjusted) * 1.1 if dp_values_microns_adjusted.size > 0 else 1000)
        ax2_adjusted = ax_adjusted.twinx()
        ax2_adjusted.plot(dp_values_microns_adjusted, plot_data_adjusted['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
        ax2_adjusted.set_ylabel('Volume/Mass Fraction', color='black', fontsize=12)
        ax2_adjusted.tick_params(axis='y', labelcolor='black')
        max_norm_fv_adjusted = max(plot_data_adjusted['volume_fraction']) if plot_data_adjusted['volume_fraction'].size > 0 else 0.1
        ax2_adjusted.set_ylim(0, max_norm_fv_adjusted * 1.2)
        lines_adjusted, labels_adjusted = ax_adjusted.get_legend_handles_labels()
        lines2_adjusted, labels2_adjusted = ax2_adjusted.get_legend_handles_labels()
        ax2_adjusted.legend(lines_adjusted + lines2_adjusted, labels_adjusted + labels2_adjusted, loc='upper left', fontsize=10)
        plt.title('Entrainment Droplet Size Distribution (After Inlet Device)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        fig_adjusted.savefig(buf_adjusted, format="png", dpi=300)
        buf_adjusted.seek(0)
        plt.close(fig_adjusted) # Close the plot to free memory


        buf_after_gravity = io.BytesIO()
        fig_after_gravity = plt.figure(figsize=(10, 6)) # Recreate figure for saving
        ax_after_gravity = fig_after_gravity.add_subplot(111)
        ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
        ax_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_after_gravity.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
        ax_after_gravity.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_after_gravity.tick_params(axis='y', labelcolor='black')
        ax_after_gravity.set_ylim(0, 1.05)
        ax_after_gravity.set_xlim(0, max(dp_values_microns_after_gravity) * 1.1 if dp_values_microns_after_gravity.size > 0 else 1000)
        ax2_after_gravity = ax_after_gravity.twinx()
        ax2_after_gravity.plot(dp_values_microns_after_gravity, plot_data_after_gravity['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
        ax2_after_gravity.set_ylabel('Volume/Mass Fraction', color='black', fontsize=12)
        ax2_after_gravity.tick_params(axis='y', labelcolor='black')
        max_norm_fv_after_gravity = max(plot_data_after_gravity['volume_fraction']) if plot_data_after_gravity['volume_fraction'].size > 0 else 0.1
        ax2_after_gravity.set_ylim(0, max_norm_fv_after_gravity * 1.2)
        lines_after_gravity, labels_after_gravity = ax_after_gravity.get_legend_handles_labels()
        lines2_after_gravity, labels2_after_gravity = ax2_after_gravity.get_legend_handles_labels()
        ax2_after_gravity.legend(lines_after_gravity + lines2_after_gravity, labels_after_gravity + labels2_after_gravity, loc='upper left', fontsize=10)
        plt.title('Entrainment Droplet Size Distribution (After Gravity Settling)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        fig_after_gravity.savefig(buf_after_gravity, format="png", dpi=300)
        buf_after_gravity.seek(0)
        plt.close(fig_after_gravity) # Close the plot to free memory

        buf_after_me = io.BytesIO()
        fig_after_me = plt.figure(figsize=(10, 6)) # Recreate figure for saving
        ax_after_me = fig_after_me.add_subplot(111)
        ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_undersize'], 'o-', label='Cumulative Volume Undersize', markersize=2, color='#1f77b4')
        ax_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['cumulative_volume_oversize'], 'o-', label='Cumulative Volume Oversize', markersize=2, color='#d62728')
        ax_after_me.set_xlabel(f'Droplet Size ({micron_unit_label})', fontsize=12)
        ax_after_me.set_ylabel('Cumulative Volume Fraction', color='black', fontsize=12)
        ax_after_me.tick_params(axis='y', labelcolor='black')
        ax_after_me.set_ylim(0, 1.05)
        ax_after_me.set_xlim(0, max(dp_values_microns_after_me) * 1.1 if dp_values_microns_after_me.size > 0 else 1000)
        ax2_after_me = ax_after_me.twinx()
        ax2_after_me.plot(dp_values_microns_after_me, plot_data_after_mist_extractor['volume_fraction'], 'o-', label='Volume/Mass Fraction', markersize=2, color='#2ca02c')
        ax2_after_me.set_ylabel('Volume/Mass Fraction', color='black', fontsize=12)
        ax2_after_me.tick_params(axis='y', labelcolor='black')
        max_norm_fv_after_me = max(plot_data_after_mist_extractor['volume_fraction']) if plot_data_after_mist_extractor['volume_fraction'].size > 0 else 0.1
        ax2_after_me.set_ylim(0, max_norm_fv_after_me * 1.2)
        lines_after_me, labels_after_me = ax_after_me.get_legend_handles_labels()
        lines2_after_me, labels2_after_me = ax2_after_me.get_legend_handles_labels()
        ax2_after_me.legend(lines_after_me + lines2_after_me, labels_after_me + labels2_after_me, loc='upper left', fontsize=10)
        plt.title('Entrainment Droplet Size Distribution (After Mist Extractor)', fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.7)
        fig_after_me.savefig(buf_after_me, format="png", dpi=300)
        buf_after_me.seek(0)
        plt.close(fig_after_me) # Close the plot to free memory

    else:
        st.warning("Please go to the 'Input Parameters' page and modify inputs to trigger calculations and generate the plot data.")



# Generate Report Page
elif page == "Generate Report":
    st.header("Download PDF Report")
    st.markdown("Enter report details below. These will appear on the PDF cover page.")

    project_name = st.text_input("Project Name", value=st.session_state.get("project_name", "Default Project"))
    client_name = st.text_input("Client Name", value=st.session_state.get("client_name", "Company Name"))
    contractor_name = st.text_input("Contractor Name", value=st.session_state.get("contractor_name", "Contractor Name"))
    equipment_tag = st.text_input("Equipment Tag", value=st.session_state.get("equipment_tag", "Equipment Tag with Description"))
    prepared_by = st.text_input("Prepared By", value=st.session_state.get("prepared_by", "Prepared By Name"))
    rev_no = st.text_input("Revision No.", value=st.session_state.get("rev_no", "Rev 0"))

    # Save to session state for persistence
    st.session_state.client_name = client_name
    st.session_state.contractor_name = contractor_name
    st.session_state.equipment_tag = equipment_tag
    st.session_state.prepared_by = prepared_by
    st.session_state.rev_no = rev_no
    st.session_state.project_name = project_name

    st.markdown("Click the button below to generate and download your full calculation report as a PDF.")

    if (
        st.session_state.plot_data_original and
        st.session_state.plot_data_adjusted and
        st.session_state.plot_data_after_gravity and
        st.session_state.plot_data_after_mist_extractor and
        st.session_state.calculation_results
    ):
        # Use a placeholder to control where the download button appears
        placeholder = st.empty()
        if placeholder.button("Prepare PDF Report"):
            with st.spinner("Generating PDF report..."):
                buf_original, buf_adjusted, buf_after_gravity, buf_after_me = create_plot_buffers(
                    st.session_state.plot_data_original,
                    st.session_state.plot_data_adjusted,
                    st.session_state.plot_data_after_gravity,
                    st.session_state.plot_data_after_mist_extractor
                )
                pdf_bytes = generate_pdf_report(
                    st.session_state.inputs,
                    st.session_state.calculation_results,
                    buf_original,
                    buf_adjusted,
                    buf_after_gravity,
                    buf_after_me,
                    st.session_state.plot_data_original,
                    st.session_state.plot_data_adjusted,
                    st.session_state.plot_data_after_gravity,
                    st.session_state.plot_data_after_mist_extractor,
                    client_name,
                    contractor_name,
                    equipment_tag,
                    prepared_by,
                    rev_no,
                    project_name,
                    
                )
                placeholder.download_button(
                    label="Download Report as PDF",
                    data=pdf_bytes,
                    file_name="Liquid in Gas Carryover Report.pdf",
                    mime="application/pdf"
                )
    else:
        st.warning("Please go to the 'Carry Over Results' page and generate the plots first.")

# Summary Page
elif page == "Summary of Results":
    # Define conversion factors
    FT_TO_MICRON = 304800
    LBM_TO_KG = 0.453592
    SCF_TO_M3 = 0.0283168
    KG_TO_LB = 2.20462
    M3_TO_US_GAL = 264.172

    # --- Initialization and Calculations for the summary table ---
    # Initialize session state variables to prevent AttributeError
    if 'results' not in st.session_state:
        st.session_state['results'] = {}
    if 'calculation_results' not in st.session_state:
        st.session_state['calculation_results'] = {}
    if 'plot_data_after_gravity' not in st.session_state:
        st.session_state['plot_data_after_gravity'] = {'total_entrained_mass_flow_rate_si': 0.0}
    if 'plot_data_adjusted' not in st.session_state:
        st.session_state['plot_data_adjusted'] = {'total_entrained_mass_flow_rate_si': 0.0}
    if 'plot_data_after_mist_extractor' not in st.session_state:
        st.session_state['plot_data_after_mist_extractor'] = {'total_entrained_mass_flow_rate_si': 0.0}

    if 'results' not in st.session_state or 'calculation_results' not in st.session_state or not st.session_state['results']:
        st.warning("Please navigate to the 'Input Parameters' & 'Flow Regime' page and run the calculation first to view Summary of results.")
    else:
        # Get flow regime (assuming you have this from a previous calculation)
        # Replace with your actual flow regime calculation if available
        flow_regime_result = st.session_state.results.get('flow_regime_td_result', 'N/A')
        st.session_state.calculation_results['flow_regime'] = flow_regime_result  
        

        # Entrainment at inlet device (correctly calculated from plot_data_adjusted)
        entrainment_after_inlet_device_kg_hr = st.session_state.plot_data_adjusted.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
        st.session_state.calculation_results['entrainment_after_inlet_device_kg_hr'] = entrainment_after_inlet_device_kg_hr

        # Calculate Inlet device efficiency (%)
        inlet_device_efficiency_percent = (100 - st.session_state.calculation_results.get('E_fraction', 0.0)*100) if st.session_state.calculation_results.get('E_fraction', 0.0) is not None else 0.0   
        st.session_state.calculation_results['inlet_device_efficiency_percent'] = inlet_device_efficiency_percent

        # Gravity Settling efficiency Calculation(GE)
        entrainment_after_inlet_kg_hr = st.session_state.plot_data_adjusted.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
        entrainment_after_gasgravity_kg_hr = st.session_state.plot_data_after_gravity.get('total_entrained_mass_flow_rate_si', 0.0)*3600
        gasgravity_extractor_efficiency = (entrainment_after_inlet_kg_hr - entrainment_after_gasgravity_kg_hr) / entrainment_after_inlet_kg_hr * 100 if entrainment_after_gasgravity_kg_hr != 0 else 0
        st.session_state.calculation_results['gasgravity_extractor_efficiency'] = gasgravity_extractor_efficiency

        # Mist extractor efficiency Calculation(ME)
        entrainment_after_gravity_kg_hr = st.session_state.plot_data_after_gravity.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
        entrainment_after_me_kg_hr = st.session_state.plot_data_after_mist_extractor.get('total_entrained_mass_flow_rate_si', 0.0) * 3600
        mist_extractor_efficiency = (entrainment_after_gravity_kg_hr - entrainment_after_me_kg_hr) / entrainment_after_gravity_kg_hr * 100 if entrainment_after_gravity_kg_hr != 0 else 0
        st.session_state.calculation_results['mist_extractor_efficiency'] = mist_extractor_efficiency

        # --- Final Gas Separator Outlet Flow (US Gal / MMSCF) ---
        gas_molecular_weight = st.session_state.inputs.get('gas_molecular_weight_input', 16.0)
        gas_mass_flow_rate_kg_s = st.session_state.inputs.get('Q_gas_mass_flow_rate_input', 1.0)
        gas_mass_flow_rate_lb_hr = gas_mass_flow_rate_kg_s * 3600 * KG_TO_LB
        gas_flow_mcf_hr = (gas_mass_flow_rate_lb_hr / gas_molecular_weight) * 379.5
        gas_flow_mmscf_hr = gas_flow_mcf_hr / 1000000

        # Convert final entrainment to US Gallons/hr
        liquid_density_kg_m3 = st.session_state.inputs.get('rho_l_input', 1000.0)
        final_entrainment_volume_m3_hr = entrainment_after_me_kg_hr / liquid_density_kg_m3
        final_entrainment_volume_us_gal_hr = final_entrainment_volume_m3_hr * M3_TO_US_GAL

        # Calculate final carryover ratio
        final_carryover_gal_mmscf = final_entrainment_volume_us_gal_hr / gas_flow_mmscf_hr if gas_flow_mmscf_hr != 0 else 0
        st.session_state.results['final_carryover_gal_mmscf'] = final_carryover_gal_mmscf

        # --- Display the summary table ---
        st.markdown("### Summary of Results")
        summary_data = {
            "Parameter": [
                "Mean particle diameter dv50 at the inlet (estimated)",
                "Inlet flow regime",
                "Entrainment at inlet pipe",
                "Entrainment after Inlet Device",
                "Inlet device efficiency",
                "Entrainment after Gravity Settling Section",
                "Gravity Settling Section efficiency",
                "Entrainment after Mist Extractor Section",
                "Mist extractor efficiency",
                "Separator Performance - Liquid carryover"
            ],
            "Value": [
                f"{st.session_state.calculation_results.get('dv50_original_fps', 0.0) * FT_TO_MICRON:.2f} um",  
                flow_regime_result,
                f"{st.session_state.calculation_results.get('E_fraction', 0.0)*100:.2f} %", 
                f"{entrainment_after_inlet_device_kg_hr:.2f} kg/hr",  
                f"{inlet_device_efficiency_percent:.2f} %",
                f"{st.session_state.plot_data_after_gravity.get('total_entrained_mass_flow_rate_si', 0.0)*3600:.3f} kg/hr", 
                f"{gasgravity_extractor_efficiency:.2f} %", 
                f"{st.session_state.plot_data_after_mist_extractor.get('total_entrained_mass_flow_rate_si', 0.0)*3600:.3f} kg/hr", 
                f"{mist_extractor_efficiency:.2f} %",
                f"{final_carryover_gal_mmscf:.4f} US Gal / MMSCF"
            ]
        }

        df = pd.DataFrame(summary_data)
        st.table(df)

    
st.markdown(r"""
---
""")















