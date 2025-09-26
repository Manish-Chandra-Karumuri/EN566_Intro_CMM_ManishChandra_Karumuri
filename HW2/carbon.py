from numpy import *
import matplotlib.pyplot as plt
import argparse

debug = False

def calculate_N_euler(N_o, tau, time, w):
    N = [N_o]
    num_step = int(time/w)
    for i in range(num_step):
        current_N = N[i]
        N_new = current_N * (1 - (w/tau))
        N.append(N_new)
        if debug: print(f'N value for iteration {i} = {N_new}')
    return N

def calculate_R_euler(N_o, tau, time, w):
    N = [N_o]
    time_points = []
    activity = []
    num_step = int(time/w)
    for i in range(num_step):
        t = i * w
        current_N = N[i]
        act = current_N / tau
        time_points.append(t)
        activity.append(act)
        N_new = current_N * (1 - (w/tau))
        N.append(N_new)
        if debug and i < 3: # Print first 3 steps
             print(f"  Step {i}: t = {t:.1f} years, N = {current_N:.4e}, Activity = {act:.4e}")
    return time_points, activity

def calculate_N_analytical(N_o, tau, t_points):
    return N_o * exp(-t_points / tau)

def calculate_R_analytical(N_o, tau, t_points):
    N_values = calculate_N_analytical(N_o, tau, t_points)
    return N_values / tau


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Calculate and plot Activity and Number of Radioactive particles present")
    parser.add_argument('--plot', type=str, dest= 'width_str' , help='Width of the step(s), comma-separated.', required=False, default='10,100')
    parser.add_argument('--mass', type=float, help='Initial mass of the radioactive material', required=False, default=1e-12)
    parser.add_argument('--mm', type=float, help='Molar mass of the radioactive element', required=False, default=0.014)
    parser.add_argument('--time', type=float, help='Duration of simulation', required=False, default= 20000)

    args = parser.parse_args()

    initial_mass = args.mass
    molar_mass = args.mm
    avo_num = 6.022 * (10**23)
    t_half = 5700

    N_o = (initial_mass/molar_mass) * avo_num

    time = args.time
    tau = t_half/log(2)

    width = [int(w.strip()) for w in args.width_str.split(',')]

    t_analytical = linspace(0, time, 1000)
    act_analytical = calculate_R_analytical(N_o, tau, t_analytical)

    for w in width:
        time_points, activity = calculate_R_euler(N_o, tau, time, w)
        plt.plot(time_points, activity, label=f'activity values for width - {w}')

    plt.plot(t_analytical, act_analytical, label = 'Analytical')
    plt.grid(True)
    plt.xlabel('Time (in years)')
    plt.ylabel('Activity (decays/year)')
    plt.title('Activity of Carbon over time')
    plt.legend()
    plt.savefig('activity_exact_vs_num.png', dpi=1000)
    plt.show()

    avg_error = []
    for w in width:
        num_step = int(time/w)
        time_points, activity = calculate_R_euler(N_o, tau, time, w)
        
        act_analytical_at_steps = calculate_R_analytical(N_o, tau, array(time_points))
        
        error =  act_analytical_at_steps - array(activity)
        if debug: print(f'Error for width {w} = {error}')
        avg_er = sum(error)/len(error)
        avg_error.append(avg_er)
        plt.plot(time_points,error, label=f'Error for width {w}')
        plt.grid(True)
        plt.legend()
        plt.xlabel('Time (in years)')
        plt.ylabel(f'Error (Analytical - Numerical)')
        plt.title(f'Error vs time (in years)')
    plt.savefig('error_vs_time.png', dpi=1000)
    plt.show()
    plt.plot(width, avg_error, color='black', label=f'Avg. Error for Widths = {width}')
    plt.scatter(width, avg_error, marker='o')
    plt.grid(True)
    plt.xlabel('Widths ------>')
    plt.ylabel('Average Error for each width')
    plt.title('Avg. Error vs Widths')
    plt.legend()
    plt.savefig('avg_err_vs_width.png', dpi=1000)
    plt.show()

    
    # For Q1(3)
    if debug:
        widths = linspace( 1000, 4000, 100)
        deviation = []
        for w in widths:
            width_3 = w
            time_3 = 2 * t_half
            time_points_3, activity_numerical_list = calculate_R_euler(N_o, tau, time, width_3)
            index_at_2T = int(time_3 / width_3)
            activity_2t_num = activity_numerical_list[index_at_2T]
            activity_2t_exact = calculate_R_analytical(N_o, tau, time_3)
            percentage_deviation = ((activity_2t_exact - activity_2t_num) / activity_2t_exact) * 100
            deviation.append(percentage_deviation)

            print(f'\n\nPercent deviation from exact and numerical after 2 Half lives = {round(percentage_deviation, 4)} %')
            print(f"Exact Analytical Activity: {activity_2t_exact:.4e} decays/year")
            print(f"Numerical Activity (Euler): {activity_2t_num:.4e} decays/year\n\n")

        plt.plot(widths, deviation, label=f'Deviation for widths = {min(widths)} : {max(widths)}')
        plt.legend()
        plt.grid(True)
        plt.xlabel('Widths')
        plt.ylabel('Percentage deviation')
        plt.title('Percentage deviation vs Widths')
        plt.savefig('percent_deviation.png')
        plt.show()
    
    width_3 = 1000
    time_3 = 2 * t_half

    time_points_3, activity_numerical_list = calculate_R_euler(N_o, tau, time, width_3)
    index_at_2T = int(time_3 / width_3)
    activity_2t_num = activity_numerical_list[index_at_2T]
    activity_2t_exact = calculate_R_analytical(N_o, tau, time_3)
    percentage_deviation = ((activity_2t_exact - activity_2t_num) / activity_2t_exact) * 100

    print(f'\n\nPercent deviation from exact and numerical after 2 Half lives for width 1000 = {round(percentage_deviation, 4)} %')
    print(f"Exact Analytical Activity: {activity_2t_exact:.4e} decays/year")
    print(f"Numerical Activity (Euler): {activity_2t_num:.4e} decays/year\n\n")
