import matplotlib.pyplot as plt
from numpy import *
import argparse


debug = False

def calculate_ideal_trajectory(u_x, u_y, delta_t, g):
    x = [0]
    y = [0]
    t = 0
    u_ynext = u_y
    t_duration = 2 * u_y / g
    if debug: print(f'Total Duration == {t_duration}')
    iter = 0
    while y[-1] >= 0:
        if t > t_duration * 1.01 and t_duration > 0: 
            break
        if t_duration == 0 and len(y) > 1:
             break
        x_next = x[-1] + (u_x * delta_t)
        y_next = y[-1] + (u_ynext * delta_t)
        if y_next < 0: break
        x.append(x_next)
        y.append(y_next)
        u_ynext = u_ynext - (g * delta_t)
        t += delta_t
        iter += 1
    if debug:  
        print(f'x values for ideal projectile - {x}') 
        print(f'y values for ideal projectile - {y}')
        print('\n\n\n')
    return x, y

def calculate_drag_trajectory(u_x, u_y, delta_t, g, rho, A, m):
    x_dr = [0]
    y_dr = [0]
    u_ydr = u_y
    u_xdr = u_x
    C = 0.5
    iter = 0
    while y_dr[-1] >= 0:
        v_mag = sqrt(u_ydr**2 + u_xdr**2)
        if debug: print(f'Iteration : {iter}  ||  Velocity Magnitude - {v_mag}')
        if v_mag == 0: break
        
        a_dr_x = -(C * rho * A * (v_mag * u_xdr)) / m
        a_dr_y = -(C * rho * A * (v_mag * u_ydr)) / m
        a_x = a_dr_x 
        a_y = -g + a_dr_y
        if debug: 
            print(f'Iteration : {iter}  ||  Acceleration with drag in x direction - {a_x}')
            print(f'Iteration : {iter}  ||  Acceleration with drag in y direction - {a_y}\n')

        x_next1 = x_dr[-1] + (u_xdr * delta_t)
        y_next1 = y_dr[-1] + (u_ydr * delta_t)
        if y_next1 < 0: break
        x_dr.append(x_next1)
        y_dr.append(y_next1)
        u_xdr += (a_x * delta_t)
        u_ydr += (a_y * delta_t)
        iter += 1
    if debug:
        print(f'x values with drag - {x_dr}')
        print(f'y values with drag - {y_dr}')
        print('\n\n\n')
    return x_dr, y_dr

def calculate_dimpled_trajectory(u_x, u_y, delta_t, g, rho, A, m):
    x_dr_dimpled = [0]
    y_dr_dimpled = [0]
    u_ydr = u_y
    u_xdr = u_x
    iter = 0
    while y_dr_dimpled[-1] >= 0:
        v_mag = sqrt(u_ydr**2 + u_xdr**2)
        if debug: print(f'Iteration : {iter}  ||  Velocity Magnitude - {v_mag}')
        if v_mag == 0: break
        
        if v_mag > 14:
            C = 7.0 / v_mag
        else:
            C = 0.5
        
        a_dr_x = -(C * rho * A * (v_mag * u_xdr)) / m
        a_dr_y = -(C * rho * A * (v_mag * u_ydr)) / m
        a_x = a_dr_x 
        a_y = -g + a_dr_y
        if debug: 
            print(f'Iteration : {iter}  ||  Acceleration with dimpled in x direction - {a_x}')
            print(f'Iteration : {iter}  ||  Acceleration with dimpled in y direction - {a_y}\n')
        x_next1 = x_dr_dimpled[-1] + (u_xdr * delta_t)
        y_next1 = y_dr_dimpled[-1] + (u_ydr * delta_t)
        if y_next1 < 0: break
        x_dr_dimpled.append(x_next1)
        y_dr_dimpled.append(y_next1)
        
        u_xdr += (a_x * delta_t)
        u_ydr += (a_y * delta_t)
        iter += 1
    if debug:
        print(f'x values with dimples - {x_dr_dimpled}')
        print(f'y values with dimples - {y_dr_dimpled}')
        print('\n\n\n')
    return x_dr_dimpled, y_dr_dimpled

def calculate_spin_trajectory(u_x, u_y, delta_t, g, rho, A, m):
    x_spin = [0]
    y_spin = [0]
    u_ydr = u_y
    u_xdr = u_x
    spin_term = 0.25 
    iter = 0
    while y_spin[-1] >= 0:
        v_mag = sqrt(u_ydr**2 + u_xdr**2)
        if debug: print(f'Iteration : {iter}  ||  Velocity Magnitude - {v_mag}')
        if v_mag == 0: break

        if v_mag > 14:
            C = 7.0 / v_mag
        else:
            C = 0.5
        
        a_dr_x = -(C * rho * A * (v_mag * u_xdr)) / m
        a_dr_y = -(C * rho * A * (v_mag * u_ydr)) / m

        a_x = a_dr_x - (spin_term * u_ydr)
        a_y = a_dr_y + (spin_term * u_xdr) - g
        if debug: 
            print(f'Iteration : {iter}  ||  Acceleration with drag and spin in x direction - {a_x}')
            print(f'Iteration : {iter}  ||  Acceleration with drag and spin in y direction - {a_y}\n')
        x_next1 = x_spin[-1] + (u_xdr * delta_t)
        y_next1 = y_spin[-1] + (u_ydr * delta_t)
        if y_next1 < 0: break
        x_spin.append(x_next1)
        y_spin.append(y_next1)
        
        u_xdr += (a_x * delta_t)
        u_ydr += (a_y * delta_t)
        iter += 1
    if debug:
        print(f'x values with drag and spin - {x_spin}')
        print(f'y values with drag and spin - {y_spin}')
        print('\n\n\n')
    return x_spin, y_spin

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and plot the trajectory of a golf ball.")
    parser.add_argument('--plot', type=str, required=True, dest='theta', help='Initial launch angle(s) in degrees, comma-separated.')
    parser.add_argument('--step', type=int, required=False, help='Number of steps to discretize x and y', default=10000)
    parser.add_argument('--mass', type=float, required=False, help='Mass of the ball', default=0.046)
    parser.add_argument('--u', type=float, required=False, help='The initial velocity of the projectile', default=70)
    parser.add_argument('--rho', type=float, required=False, help='Denisity of the air', default=1.29)
    parser.add_argument('--area', type=float, required=False, help='Frontal area of the ball', default=0.0014)
    args = parser.parse_args()

    u_ini = args.u
    g = 9.8
    rho = args.rho 
    A = args.area
    m = args.mass 

    
    angles = []

    if ',' in args.theta:
        angles = [float(a.strip()) for a in args.theta.split(',')]
        for angle in angles:
            theta_in_rad = deg2rad(angle)
            u_x = u_ini * cos(theta_in_rad)
            u_y = u_ini * sin(theta_in_rad)
            
            t_duration = 2 * u_y/g
            steps = args.step
            delta_t = t_duration / steps if t_duration > 0 else 0.01
            
            x, y = calculate_ideal_trajectory(u_x, u_y, delta_t, g)
            plt.plot(x, y, label=f"Ideal Trajectory ({angle}°)")

            x_dr, y_dr = calculate_drag_trajectory(u_x, u_y, delta_t, g, rho, A, m)
            plt.plot(x_dr, y_dr, label=f"Smooth Ball with Drag ({angle}°)")
            
            x_dimpled, y_dimpled = calculate_dimpled_trajectory(u_x, u_y, delta_t, g, rho, A, m)
            plt.plot(x_dimpled, y_dimpled, label=f"Dimpled Ball with Drag ({angle}°)")
                
            x_spin, y_spin = calculate_spin_trajectory(u_x, u_y, delta_t, g, rho, A, m)
            plt.plot(x_spin, y_spin, label=f"Dimpled Ball with Drag and Spin ({angle}°)")
    
    else:
        angle = float(args.theta)
        angles.append(angle)

        theta_in_rad = deg2rad(angle)
        u_x = u_ini * cos(theta_in_rad)
        u_y = u_ini * sin(theta_in_rad)
        
        t_duration = 2 * u_y/g
        steps = args.step
        delta_t = t_duration / steps if t_duration > 0 else 0.01
        
        x, y = calculate_ideal_trajectory(u_x, u_y, delta_t, g)
        plt.plot(x, y, label=f"Ideal Trajectory ({angle}°)")

        x_dr, y_dr = calculate_drag_trajectory(u_x, u_y, delta_t, g, rho, A, m)
        plt.plot(x_dr, y_dr, label=f"Smooth Ball with Drag ({angle}°)")
        
        x_dimpled, y_dimpled = calculate_dimpled_trajectory(u_x, u_y, delta_t, g, rho, A, m)
        plt.plot(x_dimpled, y_dimpled, label=f"Dimpled Ball with Drag ({angle}°)")
            
        x_spin, y_spin = calculate_spin_trajectory(u_x, u_y, delta_t, g, rho, A, m)
        plt.plot(x_spin, y_spin, label=f"Dimpled Ball with Drag and Spin ({angle}°)")

    plt.grid(True)
    plt.legend()
    plt.xlabel('x values calculated')
    plt.ylabel('y values calculated')
    plt.title(f'Golf ball trajectory for angles - [{args.theta}]')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.savefig('Projectile_motion.png', dpi=1000)
    plt.show()

