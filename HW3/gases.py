import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def run_simulation(grid_x, grid_y, num_steps):
    grid = np.zeros((grid_x, grid_y), dtype=int)
    
    third_x = grid_x // 3
    grid[:third_x, :] = 1  # Gas A
    grid[2*third_x:, :] = -1 # Gas B

    occupied_sites = list(zip(*np.where(grid != 0)))
    
    history = [grid.copy()]
    
    for _ in range(num_steps):
        if not occupied_sites:
            break
            
        particle_idx = np.random.randint(len(occupied_sites))
        x, y = occupied_sites[particle_idx]
        
        dx, dy = np.random.choice([-1, 1]), 0
        if np.random.rand() > 0.5:
            dx, dy = dy, dx

        nx, ny = x + dx, y + dy
        
        if 0 <= nx < grid_x and 0 <= ny < grid_y and grid[nx, ny] == 0:
            grid[nx, ny] = grid[x, y]
            grid[x, y] = 0
            occupied_sites[particle_idx] = (nx, ny)
    
    history.append(grid.copy())
    return history

def plot_densities(grids, times, filename, trials=1):
    plt.figure(figsize=(10, 7))
    for i, grid in enumerate(grids):
        n_A = np.sum(grid == 1, axis=1) / trials
        n_B = np.sum(grid == -1, axis=1) / trials
        x = np.arange(grid.shape[0])
        plt.plot(x, n_A, label=f'n_A(x) at t={times[i]}')
        plt.plot(x, n_B, '--', label=f'n_B(x) at t={times[i]}')

    plt.title(f'Linear Population Densities (Averaged over {trials} trials)')
    plt.xlabel('Position x')
    plt.ylabel('Linear Density')
    plt.legend()
    plt.grid(True)
    plt.savefig(filename)
    plt.show()

def part_1_and_2():
    grid_x, grid_y = 60, 40
    num_particles = (grid_x // 3) * grid_y * 2
    
    time_intervals = [0, num_particles * 10, num_particles * 50, num_particles * 200]
    final_grid = np.zeros((grid_x, grid_y))

    grid = np.zeros((grid_x, grid_y), dtype=int)
    third_x = grid_x // 3
    grid[:third_x, :] = 1
    grid[2*third_x:, :] = -1
    
    grids_to_plot = [grid.copy()]
    
    occupied_sites = list(zip(*np.where(grid != 0)))
    
    for step in range(1, time_intervals[-1] + 1):
        if not occupied_sites: break
        particle_idx = np.random.randint(len(occupied_sites))
        x, y = occupied_sites[particle_idx]
        
        dx, dy = np.random.choice([-1, 1]), 0
        if np.random.rand() > 0.5: dx, dy = dy, dx
        
        nx, ny = x + dx, y + dy
        
        if 0 <= nx < grid_x and 0 <= ny < grid_y and grid[nx, ny] == 0:
            grid[nx, ny] = grid[x, y]
            grid[x, y] = 0
            occupied_sites[particle_idx] = (nx, ny)
            
        if step in time_intervals[1:]:
            grids_to_plot.append(grid.copy())
            
    fig, axes = plt.subplots(1, len(grids_to_plot), figsize=(16, 4))
    for i, g in enumerate(grids_to_plot):
        axes[i].imshow(g.T, cmap='coolwarm', vmin=-1, vmax=1)
        axes[i].set_title(f'Time Step: {time_intervals[i]}')
        axes[i].set_xticks([])
        axes[i].set_yticks([])
    plt.savefig('gases_part_2_configs.png')
    plt.show()

    plot_densities(grids_to_plot, time_intervals, 'gases_part_2_densities.png')

def part_3():
    grid_x, grid_y = 60, 40
    num_trials = 100
    num_particles = (grid_x // 3) * grid_y * 2
    time_intervals = [0, num_particles * 10, num_particles * 50, num_particles * 200]
    
    avg_grids = [np.zeros((grid_x, grid_y)) for _ in time_intervals]

    for trial in range(num_trials):
        print(f"Running trial {trial + 1}/{num_trials}")
        grid = np.zeros((grid_x, grid_y), dtype=int)
        third_x = grid_x // 3
        grid[:third_x, :] = 1
        grid[2*third_x:, :] = -1
        
        avg_grids[0] += grid
        
        occupied_sites = list(zip(*np.where(grid != 0)))
        
        for step in range(1, time_intervals[-1] + 1):
            if not occupied_sites: break
            particle_idx = np.random.randint(len(occupied_sites))
            x, y = occupied_sites[particle_idx]
            
            dx, dy = np.random.choice([-1, 1]), 0
            if np.random.rand() > 0.5: dx, dy = dy, dx
            
            nx, ny = x + dx, y + dy
            
            if 0 <= nx < grid_x and 0 <= ny < grid_y and grid[nx, ny] == 0:
                grid[nx, ny] = grid[x, y]
                grid[x, y] = 0
                occupied_sites[particle_idx] = (nx, ny)
                
            if step in time_intervals[1:]:
                idx = time_intervals.index(step)
                avg_grids[idx] += grid

    plot_densities(avg_grids, time_intervals, 'gases_part_3_avg_densities.png', trials=num_trials)

def part_animate():
    print("Running animation...")
    grid_x, grid_y = 60, 40
    
    grid = np.zeros((grid_x, grid_y), dtype=int)
    third_x = grid_x // 3
    grid[:third_x, :] = 1
    grid[2*third_x:, :] = -1
    
    occupied_sites = list(zip(*np.where(grid != 0)))
    
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(grid.T, cmap='coolwarm', vmin=-1, vmax=1, animated=True)
    ax.set_title("Gas Mixing Simulation")
    ax.set_xticks([])
    ax.set_yticks([])

    steps_per_frame = 100

    def update(frame):
        for _ in range(steps_per_frame):
            if not occupied_sites:
                break
            particle_idx = np.random.randint(len(occupied_sites))
            x, y = occupied_sites[particle_idx]
            
            dx, dy = np.random.choice([-1, 1]), 0
            if np.random.rand() > 0.5:
                dx, dy = dy, dx
            
            nx, ny = x + dx, y + dy
            
            if 0 <= nx < grid_x and 0 <= ny < grid_y and grid[nx, ny] == 0:
                grid[nx, ny] = grid[x, y]
                grid[x, y] = 0
                occupied_sites[particle_idx] = (nx, ny)
        
        im.set_array(grid.T)
        return im,

    ani = animation.FuncAnimation(fig, update, frames=2000, blit=True, interval=5)

    print("Saving animation to gas_mixing.gif... This may take a moment.")
    ani.save('gas_mixing.gif', writer='pillow', fps=60)
    print("Animation saved successfully as gas_mixing.gif.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate the mixing of two gases.')
    parser.add_argument('--part', type=str, required=False, help='The part(s) to run (e.g., --part=1,2). If not specified, all parts including animation will run.')
    args = parser.parse_args()

    if args.part is None:
        print("No specific part requested. Running all parts by default.")
        part_1_and_2()
        part_3()
        part_animate()
    else:
        parts_to_run = [p.strip() for p in args.part.split(',')]
        
        animation_has_run = False
        has_run_action = False

        if any(p in ['1', '2'] for p in parts_to_run):
            part_1_and_2()
            has_run_action = True
        if '3' in parts_to_run:
            part_3()
            has_run_action = True
        if 'animate' in parts_to_run:
            part_animate()
            animation_has_run = True
            has_run_action = True

        if not has_run_action:
            print("No valid parts specified. Choose from 1, 2, 3, or animate.")
        
        if has_run_action and not animation_has_run:
            part_animate()

