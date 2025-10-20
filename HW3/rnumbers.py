import argparse
import numpy as np
import matplotlib.pyplot as plt

def generate_gaussian_box_muller(n_samples):
    u1 = np.random.rand(n_samples // 2)
    u2 = np.random.rand(n_samples // 2)
    
    z1 = np.sqrt(-2 * np.log(u1)) * np.cos(2 * np.pi * u2)
    z2 = np.sqrt(-2 * np.log(u1)) * np.sin(2 * np.pi * u2)
    
    return np.concatenate((z1, z2))

def part_1():
    num_samples_list = [1000, 1000000]
    subdivisions = [10, 20, 50, 100]

    for n_samples in num_samples_list:
        random_data = np.random.rand(n_samples)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Uniform Distribution of {n_samples} Random Numbers')
        
        for i, n_bins in enumerate(subdivisions):
            ax = axes.flatten()[i]
            ax.hist(random_data, bins=n_bins, density=True, edgecolor='black')
            ax.set_title(f'{n_bins} Subdivisions')
            ax.set_xlabel('Value')
            ax.set_ylabel('Probability Density')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'rnumbers_part_1_uniform_{n_samples}.png')
        plt.show()

def part_2():
    sigma = 1.0
    num_samples_list = [1000, 1000000]
    subdivisions = [10, 20, 50, 100]

    for n_samples in num_samples_list:
        gaussian_data = generate_gaussian_box_muller(n_samples)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Gaussian Distribution of {n_samples} Random Numbers (sigma={sigma})')
        
        x_theory = np.linspace(-4 * sigma, 4 * sigma, 200)
        y_theory = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-x_theory**2 / (2 * sigma**2))

        for i, n_bins in enumerate(subdivisions):
            ax = axes.flatten()[i]
            ax.hist(gaussian_data, bins=n_bins, density=True, range=(-4, 4), edgecolor='black')
            ax.plot(x_theory, y_theory, 'r-', lw=2, label='Theoretical Gaussian')
            ax.set_title(f'{n_bins} Subdivisions')
            ax.set_xlabel('Value')
            ax.set_ylabel('Probability Density')
            ax.legend()

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(f'rnumbers_part_2_gaussian_{n_samples}.png')
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate and analyze random numbers.')
    parser.add_argument('--part', type=str, required=True, help='The part(s) of the assignment to run (e.g., --part=1 or --part=1,2).')
    args = parser.parse_args()

    try:
        parts_to_run = [int(p.strip()) for p in args.part.split(',')]
    except ValueError:
        print(f"Error: Invalid part number in '{args.part}'. Please provide comma-separated integers.")
        exit()

    for part_num in parts_to_run:
        if part_num == 1:
            part_1()
        elif part_num == 2:
            part_2()
        else:
            print(f"Warning: Part {part_num} is not a valid option and will be skipped.")

