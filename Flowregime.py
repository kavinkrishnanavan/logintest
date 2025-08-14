import matplotlib.pyplot as plt
import numpy as np
import math

def interpolate_log(x_values, y_values, x_point):
    if x_point <= min(x_values) or x_point >= max(x_values):
        return np.interp(x_point, x_values, y_values, left=np.nan, right=np.nan)
    return np.exp(np.interp(np.log(x_point), np.log(x_values), np.log(y_values)))

def get_flow_regime(x, y, line_a, line_b, line_d):
    try:
        y_A_A = interpolate_log(line_a[0], line_a[1], x)
        y_B_B = interpolate_log(line_b[0], line_b[1], x)
        y_D_D = interpolate_log(line_d[0], line_d[1], x)

        if y < y_A_A:
            return 'Stratified'
        elif y < y_B_B and x < 100:
            return 'Dispersed Bubble'
        elif y < y_D_D:
            return 'Intermittent (Slug/Plug)'
        else:
            return 'Annular'
    except ValueError:
        return 'Out of bounds'

def compute_X(rho_L, rho_G, mu_L, mu_G, u_Ls, u_Gs, D, CL=0.046, CG=0.046, n=0.2, m=0.2):
    Re_L = rho_L * u_Ls * D / mu_L
    Re_G = rho_G * u_Gs * D / mu_G

    term_L = (4 * CL / D) * (Re_L ** -n) * (rho_L * u_Ls ** 2 / 2)
    term_G = (4 * CG / D) * (Re_G ** -m) * (rho_G * u_Gs ** 2 / 2)
    
    X = math.sqrt(term_L / term_G)
    return X

def compute_F(rho_L, rho_G, u_Gs, D, g=9.81, alpha=0):
    return math.sqrt((rho_G / (rho_L - rho_G)) * (u_Gs / math.sqrt(D * g * math.cos(np.radians(alpha)))))

def compute_K(F, D, u_Ls, nu_L):
    return F * math.sqrt(D * u_Ls / nu_L)

def plot_map(X, F, K, st_obj):
    """
    Plots the Taitel & Dukler map and the computed point.
    """
    # Data from digitized boundary curves
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

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    ax2 = ax.twinx()

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('Taitel & Dukler Flow Regime Map (1976)')
    ax.set_xlabel('X (Martinelli Parameter)')
    ax.set_ylabel('F')
    ax.set_xlim(0.001, 10000)
    ax.set_ylim(0.001, 10)
    ax.grid(True, which="both", ls="--")

    ax2.set_yscale('log')
    ax2.set_ylabel('K')
    ax2.set_ylim(1, 10000)

    ax.plot(line_a[0], line_a[1], 'k-', label='Stratified Transition (Line A-A)')
    ax.plot(line_b[0], line_b[1], 'k--', label='Dispersed Bubble Transition (Line B-B)')
    ax2.plot(x_c, k_c, 'b--', label='K-based Transition (Line C-C)')
    ax.plot(line_d[0], line_d[1], 'k-.', label='Annular Transition (Line D-D)')
    
    # Plot the calculated point

    if K < 6:
        ax2.plot(X, K, 'go', markersize=10, label=f'Computed Point (K={K:.2f})')
    else:    
        ax.plot(X, F, 'ro', markersize=10, label=f'Computed Point (X={X:.2f}, F={F:.2f})')
    # Plot K value if it's within a reasonable range for the K axis

    # Add text labels for regimes
    plt.text(0.01, 13, "STRATIFIED WAVY", fontsize=12, color='BLACK')
    plt.text(0.005, 4000, "ANNULAR - DISPERSED LIQUID", fontsize=12, color='BLACK')
    plt.text(100, 1000, "DISPERSED BUBBLE", fontsize=12, color='BLACK')
    plt.text(100, 50, "INTERMITTENT", fontsize=12, color='BLACK')
    plt.text(0.1, 2, "STRATIFIED SMOOTH", fontsize=12, color='BLACK')

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, loc='upper right')

    # Use the Streamlit object passed to the function to display the plot
    st_obj.pyplot(fig)
    return get_flow_regime(X, F, line_a, line_b, line_d)
