import argparse
import numpy as np
import matplotlib.pyplot as plt

def solve_poisson_jacobi(grid_size, a, Q, R, tolerance, max_iter=20000):
    V = np.zeros((grid_size, grid_size, grid_size))
    rho = np.zeros((grid_size, grid_size, grid_size))
    
    center = grid_size // 2
    spacing = 2 * R / (grid_size - 1)
    
    charge_pos_z = int(round(center + a / (2 * spacing)))
    charge_neg_z = int(round(center - a / (2 * spacing)))
    
    rho[center, center, charge_pos_z] = Q / spacing**3
    rho[center, center, charge_neg_z] = -Q / spacing**3
    
    iterations = 0
    error = tolerance + 1
    
    x = np.linspace(-R, R, grid_size)
    y = np.linspace(-R, R, grid_size)
    z = np.linspace(-R, R, grid_size)
    Y, X, Z = np.meshgrid(y, x, z)
    radius_sq = X**2 + Y**2 + Z**2
    boundary_mask = radius_sq >= R**2

    while error > tolerance and iterations < max_iter:
        V_new = V.copy()
        V_new[1:-1, 1:-1, 1:-1] = (1/6) * (
            V[:-2, 1:-1, 1:-1] + V[2:, 1:-1, 1:-1] +
            V[1:-1, :-2, 1:-1] + V[1:-1, 2:, 1:-1] +
            V[1:-1, 1:-1, :-2] + V[1:-1, 1:-1, 2:] +
            spacing**2 * rho[1:-1, 1:-1, 1:-1]
        )
        V_new[boundary_mask] = 0
        
        error = np.max(np.abs(V_new - V))
        V = V_new
        iterations += 1
        
    return V, iterations, center, spacing, x, z

def solve_poisson_sor(grid_size, a, Q, R, accuracy, max_iter=20000, omega=1.9):
    V = np.zeros((grid_size, grid_size, grid_size))
    rho = np.zeros((grid_size, grid_size, grid_size))
    
    center = grid_size // 2
    spacing = 2 * R / (grid_size - 1)
    
    charge_pos_z = int(round(center + a / (2 * spacing)))
    charge_neg_z = int(round(center - a / (2 * spacing)))
    
    rho[center, center, charge_pos_z] = Q / spacing**3
    rho[center, center, charge_neg_z] = -Q / spacing**3
    
    iterations = 0
    max_rel_change = accuracy + 1

    x = np.linspace(-R, R, grid_size)
    y = np.linspace(-R, R, grid_size)
    z = np.linspace(-R, R, grid_size)
    Y, X, Z = np.meshgrid(y, x, z)
    radius_sq = X**2 + Y**2 + Z**2
    boundary_mask = radius_sq >= R**2
    
    while max_rel_change > accuracy and iterations < max_iter:
        V_old_iter = V.copy()
        max_rel_change = 0.0
        for i in range(1, grid_size - 1):
            for j in range(1, grid_size - 1):
                for k in range(1, grid_size - 1):
                    if not boundary_mask[i, j, k]:
                        v_old_point = V[i,j,k]
                        term = (V[i+1,j,k] + V[i-1,j,k] + V[i,j+1,k] + V[i,j-1,k] + V[i,j,k+1] + V[i,j,k-1])
                        residue = (1/6) * (term + spacing**2 * rho[i,j,k]) - v_old_point
                        
                        V[i, j, k] += omega * residue

                        if abs(v_old_point) > 1e-9:
                            rel_change = abs(omega * residue / v_old_point)
                            if rel_change > max_rel_change:
                                max_rel_change = rel_change
        iterations += 1
        
    return iterations

def part_1(a, Q, R):
    grid_size = 51
    tolerance = 1e-5
    
    V, iters, center, spacing, x, z = solve_poisson_jacobi(grid_size, a, Q, R, tolerance)
    
    plt.figure(figsize=(10, 8))
    plt.contour(z, x, V[:, center, :], levels=np.linspace(V.min(), V.max(), 30))
    plt.title('Equipotential Lines of an Electric Dipole (x-z plane)')
    plt.xlabel('z')
    plt.ylabel('x')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.colorbar(label='Potential V')
    plt.savefig("poisson_part_1_equipotentials.png")
    plt.show()

    r = np.sqrt(x**2)
    V_r = V[center, center, :]
    
    r_fit = r[r > a]
    V_fit = V_r[r > a]
    p = Q * a
    V_analytical = p / (4 * np.pi * r_fit**2)
    
    plt.figure(figsize=(10, 6))
    plt.plot(r, V_r, 'bo-', label='Numerical Solution V(r)')
    plt.plot(r_fit, V_analytical, 'r--', label=r'Analytical $1/r^2$ Behavior')
    plt.xlabel('Distance from origin (r)')
    plt.ylabel('Potential V(r)')
    plt.title('Potential vs Distance from Origin')
    plt.legend()
    plt.grid(True)
    plt.savefig("poisson_part_1_potential_vs_r.png")
    plt.show()

def part_2(a, Q, R):
    grid_size = 31
    tolerances = np.logspace(-2, -6, 10)
    iterations_needed = []

    for tol in tolerances:
        _, iters, _, _, _, _ = solve_poisson_jacobi(grid_size, a, Q, R, tol)
        iterations_needed.append(iters)
        print(f"Tolerance: {tol:.1e}, Iterations: {iters}")

    plt.figure(figsize=(10, 6))
    plt.plot(tolerances, iterations_needed, 'o-')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Tolerance (epsilon)')
    plt.ylabel('Number of Iterations (N_iter)')
    plt.title('Iterations vs. Tolerance for Jacobi Method')
    plt.grid(True)
    plt.savefig("poisson_part_2_iters_vs_tolerance.png")
    plt.show()

def part_3(a, Q, R):
    grid_points = np.array([21, 31, 41, 51])
    accuracy = 1e-4
    iters_jacobi = []
    iters_sor = []
    
    for n in grid_points:
        print(f"Running for grid size n = {n}...")
        _, iters_j, _, _, _, _ = solve_poisson_jacobi(n, a, Q, R, tolerance=1e-7, max_iter=50000)
        iters_jacobi.append(iters_j)
        
        iters_s = solve_poisson_sor(n, a, Q, R, accuracy=accuracy)
        iters_sor.append(iters_s)

    fit_jacobi = np.polyfit(np.log(grid_points), np.log(iters_jacobi), 1)
    fit_sor = np.polyfit(np.log(grid_points), np.log(iters_sor), 1)

    print(f"Jacobi N_iter proportional to n^{fit_jacobi[0]:.2f}")
    print(f"SOR N_iter proportional to n^{fit_sor[0]:.2f}")

    plt.figure(figsize=(10, 6))
    plt.plot(grid_points, iters_jacobi, 'o-', label=f'Jacobi (slope ~ {fit_jacobi[0]:.2f})')
    plt.plot(grid_points, iters_sor, 's-', label=f'SOR (slope ~ {fit_sor[0]:.2f})')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of Grid Points (n)')
    plt.ylabel('Number of Iterations (N_iter)')
    plt.title('Iterations vs. Grid Size for Jacobi and SOR Methods')
    plt.legend()
    plt.grid(True)
    plt.savefig("poisson_part_3_iters_vs_gridsize.png")
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Solve the Poisson equation for an electric dipole.')
    parser.add_argument('--part', type=str, required=True, help='The part(s) of the assignment to run (e.g., --part=1 or --part=1,2).')
    args = parser.parse_args()

    a = 0.6
    Q_eps0 = 1.0
    R = 10.0
    
    try:
        parts_to_run = [int(p.strip()) for p in args.part.split(',')]
    except ValueError:
        print(f"Error: Invalid part number in '{args.part}'. Please provide comma-separated integers.")
        exit()

    for part_num in parts_to_run:
        if part_num == 1:
            part_1(a, Q_eps0, R)
        elif part_num == 2:
            part_2(a, Q_eps0, R)
        elif part_num == 3:
            part_3(a, Q_eps0, R)
        else:
            print(f"Warning: Part {part_num} is not a valid option and will be skipped.")

