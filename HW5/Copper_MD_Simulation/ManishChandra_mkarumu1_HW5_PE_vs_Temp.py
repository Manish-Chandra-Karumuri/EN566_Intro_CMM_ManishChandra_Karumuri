import matplotlib.pyplot as plt
import numpy as np

heat_t = []
heat_pe = []
cool_t = []
cool_pe = []

phase = 0 

print("Reading log.lammps...")

with open('log.lammps', 'r') as f:
    for line in f:
        if "fix heat" in line:
            phase = 1
        elif "fix cool" in line:
            phase = 2
        
        parts = line.split()
        

        if len(parts) >= 4:
            try:
                step = int(parts[0])
                t = float(parts[1])
                pe = float(parts[3])
                
                
                if phase == 1:
                    heat_t.append(t)
                    heat_pe.append(pe)
                elif phase == 2:
                    cool_t.append(t)
                    cool_pe.append(pe)
            except ValueError:
                
                continue

heat_t_arr = np.array(heat_t)
heat_pe_arr = np.array(heat_pe)

if len(heat_t_arr) > 0:
    # Find melting point via max derivative
    gradients = np.gradient(heat_pe_arr, heat_t_arr)
    tm_index = np.argmax(np.abs(gradients))
    melting_point = heat_t_arr[tm_index]

    print(f"Estimated Melting Point: {melting_point:.2f} K")

    
    plt.figure(figsize=(10, 6))
    plt.plot(heat_t, heat_pe, 'r-', label='Heating')
    plt.plot(cool_t, cool_pe, 'b--', label='Cooling')
    
    
    plt.axvline(x=melting_point, color='k', linestyle=':', label=f'Tm = {melting_point:.0f} K')
    
    plt.xlabel('Temperature (K)')
    plt.ylabel('Potential Energy (eV)')
    plt.title('Copper Melting and Cooling Cycle')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig('Q2_PE_vs_Temp.png')
    print("Plot saved as Q2_PE_vs_Temp.png")
    plt.show()
else:
    print("Error: No heating data extraction. Check if 'fix heat' is in the log file.")

