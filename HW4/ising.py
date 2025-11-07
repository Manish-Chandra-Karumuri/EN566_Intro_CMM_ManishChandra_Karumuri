import matplotlib.pyplot as plt
from numpy import *
import argparse
import sys
import os

J = 1.5
K_B = 1.0

def get_neighbors(n, i, j):
    return [
        ((i + 1) % n, j),
        ((i - 1) % n, j),
        (i, (j + 1) % n),
        (i, (j - 1) % n)
    ]

def initialize_lattice(n):
    return random.choice([-1, 1], size=(n, n))

def calculate_total_energy(lattice, n, J):
    energy = 0
    for i in range(n):
        for j in range(n):
            spin = lattice[i, j]
            neighbor_sum = lattice[i, (j + 1) % n] + lattice[(i + 1) % n, j]
            energy += -J * spin * neighbor_sum
    return energy

def metropolis_step(lattice, n, J, T, current_energy, current_magnetization):
    i, j = random.randint(n), random.randint(n)
    spin = lattice[i, j]
    neighbor_sum = sum(lattice[ni, nj] for ni, nj in get_neighbors(n, i, j))
    delta_e = 2 * J * spin * neighbor_sum

    if delta_e <= 0 or random.rand() < exp(-delta_e / (K_B * T)):
        lattice[i, j] *= -1
        current_energy += delta_e
        current_magnetization += 2 * (-spin)
    
    return current_energy, current_magnetization

def run_simulation(n, J, T, thermal_sweeps, measure_sweeps):
    N = n * n
    lattice = initialize_lattice(n)
    E = calculate_total_energy(lattice, n, J)
    M = sum(lattice)

    for _ in range(thermal_sweeps):
        for _ in range(N):
            E, M = metropolis_step(lattice, n, J, T, E, M)

    energy_sum = 0.0
    energy_sq_sum = 0.0
    magnetization_abs_sum = 0.0

    for _ in range(measure_sweeps):
        for _ in range(N):
            E, M = metropolis_step(lattice, n, J, T, E, M)

        energy_sum += E
        energy_sq_sum += E**2
        magnetization_abs_sum += abs(M)

    total_measurements = measure_sweeps
    avg_E = energy_sum / total_measurements
    avg_E_sq = energy_sq_sum / total_measurements
    avg_M_abs = magnetization_abs_sum / total_measurements
    
    return avg_E, avg_E_sq, avg_M_abs

def run_part1(n_fixed=50):
    print(f"\n--- Running Part 1: Magnetization vs. T for n={n_fixed} ---")
    
    T_values = linspace(1.0, 5.0, 41) 
    M_per_N_values = []
    
    THERMAL_SWEEPS = 2000 
    MEASURE_SWEEPS = 5000 
    
    for T in T_values:
        if T < 1e-6: T = 1e-6
        print(f"Simulating T = {T:.2f}...")
        
        avg_E, avg_E_sq, avg_M_abs = run_simulation(n_fixed, J, T, THERMAL_SWEEPS, MEASURE_SWEEPS)
        
        M_total = avg_M_abs
        M_per_N_values.append(M_total / (n_fixed * n_fixed))
        
    plt.figure()
    plt.plot(T_values, M_per_N_values, 'o-', markersize=4, label=f'n={n_fixed} ({n_fixed*n_fixed} spins)')
    
    T_c_analytic = 2 * J / log(1 + sqrt(2))
    plt.axvline(T_c_analytic, color='r', linestyle='--', label=f'Onsager $T_C$ ({T_c_analytic:.3f})')
    
    plt.grid(True)
    plt.legend()
    plt.xlabel('Temperature (T)')
    plt.ylabel('Magnetization per Spin ($M/N$)')
    plt.title(f'2D Ising Model: Magnetization vs. Temperature (n={n_fixed})')
    plt.savefig('Ising_Magnetization_Part1.png', dpi=300)
    plt.show()

def run_part2():
    print(f"\n--- Running Part 2: Specific Heat and Finite-Size Scaling ---")

    n_values = [5, 10, 20, 30, 40, 50, 75, 100, 200, 500] 
    T_c_analytic = 2 * J / log(1 + sqrt(2))
    
    T_sweep = linspace(T_c_analytic - 0.2, T_c_analytic + 0.2, 51)
    
    THERMAL_SWEEPS = 2000
    MEASURE_SWEEPS = 5000

    C_max_N_data = []
    
    plt.figure()

    for n in n_values:
        print(f"\nSimulating Specific Heat for n={n}...")
        C_per_N_values = []
        
        if n < 30:
            T_sweep_n = linspace(1.0, 5.0, 51)
        else:
            T_sweep_n = T_sweep

        for T in T_sweep_n:
            if T < 1e-6: T = 1e-6 
            
            avg_E, avg_E_sq, avg_M_abs = run_simulation(n, J, T, THERMAL_SWEEPS, MEASURE_SWEEPS)
            N = n * n
            
            C = (avg_E_sq - avg_E**2) / (K_B * T**2)
            C_per_N_values.append(C / N)
            
        C_max_N = max(C_per_N_values)
        T_max_idx = C_per_N_values.index(C_max_N)
        T_max = T_sweep_n[T_max_idx]
        
        print(f"  n={n}: C_max/N = {C_max_N:.4f} at T = {T_max:.3f}")
        C_max_N_data.append((n, C_max_N))
        
        if n in [10, 50, 200]:
            plt.plot(T_sweep_n, C_per_N_values, label=f'C/N vs T (n={n})')

    plt.grid(True)
    plt.legend()
    plt.xlabel('Temperature (T)')
    plt.ylabel('Specific Heat per Spin ($C/N$)')
    plt.title('2D Ising Model: Specific Heat vs. Temperature Samples')
    plt.savefig('Ising_SpecificHeat_Samples_Part2.png', dpi=300)
    plt.show()

    n_for_plot = array([d[0] for d in C_max_N_data])
    C_max_N_for_plot = array([d[1] for d in C_max_N_data])

    plt.figure()
    plt.plot(log(n_for_plot), C_max_N_for_plot, 's--', label='Simulation Data')
    
    plt.grid(True)
    plt.legend()
    plt.xlabel('Log of Lattice Size ($\ln(n)$)')
    plt.ylabel('Max Specific Heat per Spin ($C_{max}/N$)')
    plt.title('2D Ising Model: Finite-Size Scaling Verification ($C_{max}/N \\sim \\ln(n)$)')
    plt.savefig('Ising_FiniteSizeScaling_Part2.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Numerically study the 2D Ising Model using the Metropolis algorithm.")
    parser.add_argument('--part', type=int, required=True, help='Which part of the bonus question to run (1 or 2).')
    
    args = parser.parse_args()
    
    if args.part == 1:
        run_part1()
    elif args.part == 2:
        run_part2()
    else:
        print("Error: Please specify --part 1 or --part 2.")
        sys.exit(1)