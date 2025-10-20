import argparse
import numpy as np
import matplotlib.pyplot as plt

def perform_walk(num_steps):
    steps = np.random.randint(0, 4, num_steps)
    x = np.zeros(num_steps + 1)
    y = np.zeros(num_steps + 1)
    
    for i in range(num_steps):
        if steps[i] == 0:
            x[i+1] = x[i] + 1
            y[i+1] = y[i]
        elif steps[i] == 1:
            x[i+1] = x[i] - 1
            y[i+1] = y[i]
        elif steps[i] == 2:
            y[i+1] = y[i] + 1
            x[i+1] = x[i]
        else: # steps[i] == 3
            y[i+1] = y[i] - 1
            x[i+1] = x[i]
            
    return x, y

def part_1(num_walks, num_steps):
    all_x = np.zeros((num_walks, num_steps + 1))
    
    for i in range(num_walks):
        x, _ = perform_walk(num_steps)
        all_x[i, :] = x

    mean_x = np.mean(all_x, axis=0)
    mean_x_sq = np.mean(all_x**2, axis=0)
    
    steps = np.arange(num_steps + 1)
    
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.plot(steps, mean_x, 'b-')
    plt.title(r'Mean Position $\langle x_n \rangle$')
    plt.xlabel('Number of Steps (n)')
    plt.ylabel(r'$\langle x_n \rangle$')
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(steps, mean_x_sq, 'r-')
    plt.title(r'Mean Square Position $\langle x_n^2 \rangle$')
    plt.xlabel('Number of Steps (n)')
    plt.ylabel(r'$\langle x_n^2 \rangle$')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('rwalk_part_1_mean_positions.png')
    plt.show()

def part_2(num_walks, num_steps):
    all_x = np.zeros((num_walks, num_steps + 1))
    all_y = np.zeros((num_walks, num_steps + 1))
    
    for i in range(num_walks):
        x, y = perform_walk(num_steps)
        all_x[i, :] = x
        all_y[i, :] = y
        
    mean_r_sq = np.mean(all_x**2 + all_y**2, axis=0)
    steps = np.arange(num_steps + 1)
    
    fit = np.polyfit(steps[1:], mean_r_sq[1:], 1)
    diffusion_constant = fit[0] / 4.0
    print('\n--------------------------\n')
    print(f"The motion is diffusive.")
    print(f"Fit for <r^2> = D_fit * n: D_fit = {fit[0]:.4f}")
    print(f"Estimated Diffusion Constant D = D_fit / 4 = {diffusion_constant:.4f}")
    print('\n--------------------------\n')
    
    plt.figure(figsize=(10, 7))
    plt.plot(steps, mean_r_sq, 'bo', label='Simulation Data')
    plt.plot(steps, fit[0] * steps + fit[1], 'r--', label=f'Linear Fit (slope={fit[0]:.3f})')
    plt.title(r'Mean Square Distance $\langle r^2 \rangle$ vs. Time (steps)')
    plt.xlabel('Number of Steps (n)')
    plt.ylabel(r'$\langle r^2 \rangle$')
    plt.grid(True)
    plt.legend()
    plt.savefig('rwalk_part_2_diffusion.png')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate a 2D random walk.')
    parser.add_argument('--part', type=str, required=True, help='The part(s) of the assignment to run (e.g., --part=1 or --part=1,2).')
    args = parser.parse_args()
    
    num_walks = 10000
    num_steps = 100

    try:
        parts_to_run = [int(p.strip()) for p in args.part.split(',')]
    except ValueError:
        print(f"Error: Invalid part number in '{args.part}'. Please provide comma-separated integers.")
        exit()

    for part_num in parts_to_run:
        if part_num == 1:
            part_1(num_walks, num_steps)
        elif part_num == 2:
            part_2(num_walks, num_steps)
        else:
            print(f"Warning: Part {part_num} is not a valid option and will be skipped.")

