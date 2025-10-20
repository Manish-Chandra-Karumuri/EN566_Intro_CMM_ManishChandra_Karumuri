import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os

def get_acceleration(theta, omega, t, Omega_D, alpha_D, gamma, g, l, is_nonlinear):
    driving_force = alpha_D * np.sin(Omega_D * t)
    damping_force = -2 * gamma * omega
    if is_nonlinear:
        restoring_force = -(g / l) * np.sin(theta)
    else:
        restoring_force = -(g / l) * theta
    return restoring_force + damping_force + driving_force

def euler_cromer_solver(Omega_D, alpha_D, gamma, g, l, is_nonlinear, t_max, dt, theta0=0.0, omega0=0.0):
    num_steps = int(t_max / dt)
    t = np.linspace(0, t_max, num_steps + 1)
    theta = np.zeros(num_steps + 1)
    omega = np.zeros(num_steps + 1)

    theta[0] = theta0
    omega[0] = omega0

    for i in range(num_steps):
        accel = get_acceleration(theta[i], omega[i], t[i], Omega_D, alpha_D, gamma, g, l, is_nonlinear)
        omega[i+1] = omega[i] + accel * dt
        theta[i+1] = theta[i] + omega[i+1] * dt
        
        if theta[i+1] > np.pi:
            theta[i+1] -= 2*np.pi
        if theta[i+1] < -np.pi:
            theta[i+1] += 2*np.pi
            
    return t, theta, omega

def rk4_solver(Omega_D, alpha_D, gamma, g, l, is_nonlinear, t_max, dt, theta0=0.0, omega0=0.0):
    num_steps = int(t_max / dt)
    t = np.linspace(0, t_max, num_steps + 1)
    theta = np.zeros(num_steps + 1)
    omega = np.zeros(num_steps + 1)

    theta[0] = theta0
    omega[0] = omega0

    for i in range(num_steps):
        ti = t[i]
        theta_i = theta[i]
        omega_i = omega[i]

        k1_theta = dt * omega_i
        k1_omega = dt * get_acceleration(theta_i, omega_i, ti, Omega_D, alpha_D, gamma, g, l, is_nonlinear)

        k2_theta = dt * (omega_i + 0.5 * k1_omega)
        k2_omega = dt * get_acceleration(theta_i + 0.5 * k1_theta, omega_i + 0.5 * k1_omega, ti + 0.5 * dt, Omega_D, alpha_D, gamma, g, l, is_nonlinear)

        k3_theta = dt * (omega_i + 0.5 * k2_omega)
        k3_omega = dt * get_acceleration(theta_i + 0.5 * k2_theta, omega_i + 0.5 * k2_omega, ti + 0.5 * dt, Omega_D, alpha_D, gamma, g, l, is_nonlinear)

        k4_theta = dt * (omega_i + k3_omega)
        k4_omega = dt * get_acceleration(theta_i + k3_theta, omega_i + k3_omega, ti + dt, Omega_D, alpha_D, gamma, g, l, is_nonlinear)

        theta[i+1] = theta_i + (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta) / 6
        omega[i+1] = omega_i + (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega) / 6
        
        if theta[i+1] > np.pi:
            theta[i+1] -= 2*np.pi
        if theta[i+1] < -np.pi:
            theta[i+1] += 2*np.pi

    return t, theta, omega

def part_1(g, l, gamma):
    print('\n-------------------------------------------------\n')
    print("Running Part 1: Analytical Resonance Frequency")
    omega_0_sq = g / l
    omega_res_sq = omega_0_sq - 2 * gamma**2
    if omega_res_sq > 0:
        omega_res = np.sqrt(omega_res_sq)
        print(f"Natural frequency (omega_0): {np.sqrt(omega_0_sq):.4f} rad/s")
        print(f"Analytical resonance frequency (Omega_D): {omega_res:.4f} rad/s")
    else:
        print("The system is overdamped, resonance does not occur in the same way.")
    print("\nThe small-angle approximation should be good as long as the amplitude is small.")
    print("However, the driving force and damping will determine the steady-state amplitude.")
    print('\n-------------------------------------------------\n')

def part_2(g, l, gamma, alpha_D):
    print('\n-------------------------------------------------\n')
    print("Running Part 2: Resonance Structure")
    
    t_max = 200.0
    dt = 0.04
    
    driving_frequencies = np.linspace(0.5, 1.5, 50)
    amplitudes_ec = []
    amplitudes_rk4 = []
    phase_shifts_rk4 = []

    for omega_d in driving_frequencies:
        # Euler-Cromer
        t_ec, theta_ec, _ = euler_cromer_solver(omega_d, alpha_D, gamma, g, l, False, t_max, dt)
        steady_state_mask_ec = t_ec > t_max * 0.6
        amplitudes_ec.append(np.max(theta_ec[steady_state_mask_ec]))

        # Runge-Kutta 4
        t_rk4, theta_rk4, _ = rk4_solver(omega_d, alpha_D, gamma, g, l, False, t_max, dt)
        steady_state_mask_rk4 = t_rk4 > t_max * 0.6
        amplitudes_rk4.append(np.max(theta_rk4[steady_state_mask_rk4]))

        # Phase Shift Calculation (using more accurate RK4)
        steady_theta = theta_rk4[steady_state_mask_rk4]
        steady_t = t_rk4[steady_state_mask_rk4]
        driving_force_signal = alpha_D * np.sin(omega_d * steady_t)

        theta_peaks, _ = find_peaks(steady_theta)
        force_peaks, _ = find_peaks(driving_force_signal)

        if len(theta_peaks) > 0 and len(force_peaks) > 0:
            last_force_peak_time_index = np.searchsorted(steady_t[force_peaks], steady_t[theta_peaks[0]]) - 1
            if last_force_peak_time_index < 0:
                 phase_shifts_rk4.append(np.nan)
                 continue
            last_force_peak_time = steady_t[force_peaks[last_force_peak_time_index]]
            time_lag = steady_t[theta_peaks[0]] - last_force_peak_time
            phase_shift = (time_lag * omega_d) % (2 * np.pi)
            if phase_shift > np.pi:
                 phase_shift -= 2 * np.pi
            phase_shifts_rk4.append(phase_shift)
        else:
            phase_shifts_rk4.append(np.nan)

    plt.figure(figsize=(12, 5))
    plt.plot(driving_frequencies, amplitudes_ec, 'bo-', label='Euler-Cromer Amplitude')
    plt.plot(driving_frequencies, amplitudes_rk4, 'rx-', label='Runge-Kutta 4 Amplitude')
    plt.title('Resonance Curve')
    plt.xlabel(r'Driving Frequency $\Omega_D$ (rad/s)')
    plt.ylabel(r'Steady-State Amplitude $\theta_0$ (rad)')
    plt.grid(True)
    plt.legend()
    plt.savefig("part_2_resonance_curve.png")
    plt.show()
    
    # FWHM Calculation (more robust method)
    amplitudes = np.array(amplitudes_rk4)
    peak_amp = np.max(amplitudes)
    peak_freq = driving_frequencies[np.argmax(amplitudes)]
    half_max = peak_amp / 2.0

    # Find crossings on either side of the peak
    try:
        # Left side
        left_idx = np.where(amplitudes[:np.argmax(amplitudes)] > half_max)[0][0]
        freq1 = np.interp(half_max, [amplitudes[left_idx-1], amplitudes[left_idx]], [driving_frequencies[left_idx-1], driving_frequencies[left_idx]])

        # Right side
        right_idx = np.where(amplitudes[np.argmax(amplitudes):] < half_max)[0][0] + np.argmax(amplitudes)
        freq2 = np.interp(half_max, [amplitudes[right_idx], amplitudes[right_idx-1]], [driving_frequencies[right_idx], driving_frequencies[right_idx-1]])
        
        fwhm = freq2 - freq1
        
        print(f"\nMaximum Amplitude: {peak_amp:.4f} rad at {peak_freq:.4f} rad/s")
        print(f"Half-Maximum Amplitude: {half_max:.4f} rad")
        print(f"Frequencies at half-max: {freq1:.4f} and {freq2:.4f} rad/s")
        print(f"Numerically extracted FWHM: {fwhm:.4f} rad/s")
        print(f"Value of gamma * 2: {2*gamma:.4f} rad/s (Expected FWHM for small damping is ~2*gamma)")
        print('\n-------------------------------------------------\n')

    except IndexError:
        print("\nCould not determine FWHM reliably. The curve may not cross the half-max threshold on both sides within the sampled frequency range.")


    plt.figure(figsize=(12, 5))
    plt.plot(driving_frequencies, phase_shifts_rk4, 'go-', label='Runge-Kutta 4 Phase Shift')
    plt.title('Phase Shift vs. Driving Frequency')
    plt.xlabel(r'Driving Frequency $\Omega_D$ (rad/s)')
    plt.ylabel(r'Phase Shift $\phi$ (rad)')
    plt.grid(True)
    plt.legend()
    plt.savefig("part_2_phase_shift.png")
    plt.show()

    resonance_freq_rk4 = driving_frequencies[np.nanargmax(amplitudes_rk4)]
    t, theta_rk4, omega_rk4 = rk4_solver(resonance_freq_rk4, alpha_D, gamma, g, l, False, t_max, dt)
    
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 1, 1)
    plt.plot(t, theta_rk4, label=fr'$\theta(t)$ at $\Omega_D = {resonance_freq_rk4:.3f}$')
    plt.title('Runge-Kutta: Motion at Resonance')
    plt.ylabel(r'Angle $\theta$ (rad)')
    plt.grid(True)
    plt.legend()
    
    plt.subplot(2, 1, 2)
    plt.plot(t, omega_rk4, label=fr'$\omega(t)$ at $\Omega_D = {resonance_freq_rk4:.3f}$', color='orange')
    plt.ylabel(r'Angular Velocity $\omega$ (rad/s)')
    plt.xlabel('Time (s)')
    plt.grid(True)
    plt.legend()
    plt.savefig("part_2_motion_at_resonance.png")
    plt.show()


def part_3(g, l, gamma, alpha_D):
    print('\n-------------------------------------------------\n')
    print("Running Part 3: Energy Analysis")
    omega_res = np.sqrt((g/l) - 2*gamma**2)
    t_max = 10 * (2 * np.pi / omega_res)
    dt = 0.01

    t, theta, omega = rk4_solver(omega_res, alpha_D, gamma, g, l, False, t_max, dt)
    
    m = 1.0 
    ke = 0.5 * m * (l * omega)**2
    pe = m * g * l * (1 - np.cos(theta))
    total_e = ke + pe
    
    plt.figure(figsize=(12, 6))
    plt.plot(t, ke, label='Kinetic Energy')
    plt.plot(t, pe, label='Potential Energy')
    plt.plot(t, total_e, label='Total Energy')
    plt.title('Energy of the Pendulum near Resonance')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (arbitrary units)')
    plt.grid(True)
    plt.legend()
    plt.savefig("part_3_energy_analysis.png")
    plt.show()
    

def part_4(g, l, gamma, alpha_D):
    print('\n-------------------------------------------------\n')
    print("Running Part 4: Non-Linear Effects")
    omega_res = np.sqrt(g/l)
    t_max = 100.0
    dt = 0.04
    
    t_lin, theta_lin, omega_lin = rk4_solver(omega_res, alpha_D, gamma, g, l, False, t_max, dt)
    t_nonlin, theta_nonlin, omega_nonlin = rk4_solver(omega_res, alpha_D, gamma, g, l, True, t_max, dt)
    
    plt.figure(figsize=(12, 6))
    plt.plot(t_lin, theta_lin, label='Linear')
    plt.plot(t_nonlin, theta_nonlin, '--', label='Non-Linear')
    plt.title(fr'Linear vs Non-Linear Pendulum ($\alpha_D = {alpha_D}$)')
    plt.xlabel('Time (s)')
    plt.ylabel(r'Angle $\theta$ (rad)')
    plt.grid(True)
    plt.legend()
    plt.savefig(f"part_4_linear_vs_nonlinear_alpha{alpha_D}.png")
    plt.show()

    alpha_D_high = 1.2
    t_high, theta_high, omega_high = rk4_solver(omega_res, alpha_D_high, gamma, g, l, True, t_max, dt)
    
    plt.figure(figsize=(12, 6))
    plt.plot(t_high, theta_high, label=fr'Non-Linear with $\alpha_D = {alpha_D_high}$')
    plt.title('Non-Linear Pendulum with High Driving Force')
    plt.xlabel('Time (s)')
    plt.ylabel(r'Angle $\theta$ (rad)')
    plt.grid(True)
    plt.legend()
    plt.savefig(f"part_4_nonlinear_high_alpha.png")
    plt.show()

def part_5(g, l, gamma):
    print('\n-------------------------------------------------\n')
    print("Running Part 5: Chaos and Lyapunov Exponent")
    Omega_D = 2.0/3.0
    t_max = 300.0
    dt = 0.04
    d_theta_0 = 0.001
    
    alphas = [0.2, 0.5, 1.2]
    
    plt.figure(figsize=(10, 7))
    for alpha_d_val in alphas:
        t, theta1, _ = rk4_solver(Omega_D, alpha_d_val, gamma, g, l, True, t_max, dt, theta0=0.2)
        t, theta2, _ = rk4_solver(Omega_D, alpha_d_val, gamma, g, l, True, t_max, dt, theta0=0.2 + d_theta_0)
        
        delta_theta = np.abs(theta1 - theta2)
        
        plt.semilogy(t, delta_theta, label=fr'$\alpha_D = {alpha_d_val}$')

    plt.title('Divergence of Trajectories (Lyapunov Exponent Estimation)')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$|\Delta\theta(t)|$ (rad)')
    plt.grid(True)
    plt.legend()
    plt.savefig("part_5_lyapunov_estimation.png")
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulate a driven, damped pendulum.')
    parser.add_argument('--part', type=str, required=True, help='The part(s) of the assignment to run (e.g., --part=2 or --part=2,3,5).')
    args = parser.parse_args()

    g = 9.8
    l = 9.8
    gamma = 0.25
    alpha_D = 0.2
    
    try:
        parts_to_run = [int(p.strip()) for p in args.part.split(',')]
    except ValueError:
        print(f"Error: Invalid part number in '{args.part}'. Please provide comma-separated integers.")
        exit()

    for part_num in parts_to_run:
        if part_num == 1:
            part_1(g, l, gamma)
        elif part_num == 2:
            part_2(g, l, gamma, alpha_D)
        elif part_num == 3:
            part_3(g, l, gamma, alpha_D)
        elif part_num == 4:
            part_4(g, l, gamma, alpha_D)
        elif part_num == 5:
            part_5(g, l, gamma)
        else:
            print(f"Warning: Part {part_num} is not a valid option and will be skipped.")

