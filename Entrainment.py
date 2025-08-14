import numpy as np
import math

# The Schiller-Naumann correlation is used instead of a data table.
# This correlation is valid for Re_p < 1000.

def find_cd_from_rep(re_p):

    if re_p <= 0:
        return float('inf') # Return a very large value to handle non-positive Re_p

    # Apply the Schiller-Naumann correlation
    cd = (24.0 / re_p) * (1.0 + 0.15 * math.pow(re_p, 0.687)) + (0.42 / (1 + (42500 / re_p**1.16)))
    return cd


def run_all_calculations():
    """
    Main function to orchestrate the entire calculation process.
    """
    print("Welcome to the Complete E Value Calculator!")
    print("Please enter the values for the following parameters (SI units).")
    print("-" * 40)

    try:
        # Hardcoded A_2 as per user request
        a2 = 9e-8  #  Refer to Pan and Hanratty (2002)
        #a2 = 8.8e-5  # Large Gas velocity U_g
        
        # Get inputs from the user
        g = 9.81 # float(input("Enter acceleration due to gravity (g): "))
        rho_l = 800 #float(input("Enter liquid density (ρ_L): "))
        rho_g = 70 #float(input("Enter gas density (ρ_G): "))
        mu_g = 0.00001 #float(input("Enter gas viscosity (μ_G): "))
        mu_l = 0.0001 #float(input("Enter liquid viscosity (μ_L): "))
        d_capital = 0.6 #float(input("Enter the value for D: "))
        ug = 10 #float(input("Enter the value for U_G: "))
        sigma = 0.0115 #float(input("Enter the value for sigma (σ): "))
        w_l = 0.5 #float(input("Enter the value for W_L: "))
        
        print("-" * 40)
        
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

        # Print the final result for E
        print("-" * 40)
        print(f"The final calculated value of E is: {e_value:.6f}")

    except ValueError as ve:
        print(f"Invalid input: {ve}. Please enter numerical values.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Run the full calculation function
if __name__ == "__main__":
    run_all_calculations()
