import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def solve_diffusion(L, T, D, nx, nt):
    dx = L / (nx - 1)
    dt = T / nt
    
    alpha = D * dt / dx**2
    if alpha > 0.5:
        print(f"Warning: Stability condition not met. alpha = {alpha:.3f} > 0.5")

    x = np.linspace(-L/2, L/2, nx)
    rho = np.zeros(nx)

    box_width_sites = int(nx / 20)
    center = nx // 2
    rho[center - box_width_sites : center + box_width_sites] = 1.0
    rho /= np.sum(rho) * dx
    
    rho_history = [rho.copy()]
    time_points = [0.0]
    
    time_snapshots = [T/10, T/4, T/2, T*0.75, T]
    next_snapshot_idx = 0

    for j in range(nt):
        rho_new = rho.copy()
        for i in range(1, nx - 1):
            rho_new[i] = rho[i] + alpha * (rho[i+1] - 2*rho[i] + rho[i-1])
        rho = rho_new
        
        current_time = (j + 1) * dt
        if next_snapshot_idx < len(time_snapshots) and current_time >= time_snapshots[next_snapshot_idx]:
            rho_history.append(rho.copy())
            time_points.append(current_time)
            next_snapshot_idx += 1
            
    return x, rho_history, time_points

def part_2():
    L = 20.0
    T = 2.0
    D = 2.0
    nx = 201
    nt = 4000
    
    x, rho_history, time_points = solve_diffusion(L, T, D, nx, nt)
    
    plt.figure(figsize=(10, 7))
    print('\n-----------------------------\n')
    print("Fitting results:")
    print("-" * 40)
    print(f"{'Time (s)':<10} | {'Fit sigma':<12} | {'Theory sigma':<15}")
    print("-" * 40)

    


    for i in range(len(time_points)):
        rho = rho_history[i]
        t = time_points[i]
        
        popt, _ = curve_fit(gaussian, x, rho, p0=[np.max(rho), 0, 1])
        fit_sigma = popt[2]
        
        if t > 0:
            theory_sigma = np.sqrt(2 * D * t)
            print(f"{t:<10.3f} | {fit_sigma:<12.4f} | {theory_sigma:<15.4f}")
        else:
            print(f"{t:<10.3f} | {fit_sigma:<12.4f} | {'N/A'}")
    
      
        plt.plot(x, rho, label=f't = {t:.2f} s')
        plt.plot(x, gaussian(x, *popt), '--', color='black', alpha=0.5)
    print('\n-----------------------------\n')
    plt.title('1D Diffusion of a Sharply Peaked Profile')
    plt.xlabel('Position x')
    plt.ylabel('Density rho(x, t)')
    plt.legend()
    plt.grid(True)
    plt.savefig('diffusion_part_2_profiles.png')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Solve the 1D diffusion equation.')
    parser.add_argument('--part', type=str, required=True, help='The part(s) to run (e.g., --part=2).')
    args = parser.parse_args()

    try:
        parts_to_run = [int(p.strip()) for p in args.part.split(',')]
    except ValueError:
        print(f"Error: Invalid part number in '{args.part}'. Please provide comma-separated integers.")
        exit()

    for part_num in parts_to_run:
        if part_num == 2:
            part_2()
        else:
            print(f"Warning: Part {part_num} is not a valid option and will be skipped.")

