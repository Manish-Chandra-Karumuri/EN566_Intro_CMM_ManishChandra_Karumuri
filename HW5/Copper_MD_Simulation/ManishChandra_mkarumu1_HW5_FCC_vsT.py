### ---------- FCC analysis part --------------

import matplotlib.pyplot as plt

# --- SIMULATION PARAMETERS (Update these to match your input file) ---
T_MIN = 293.15
T_MAX = 1700.0       
HEAT_START = 20000
RUN_STEPS = 1400000  
HOLD_STEPS = 10000


HEAT_END = HEAT_START + RUN_STEPS
COOL_START = HEAT_END + HOLD_STEPS
COOL_END = COOL_START + RUN_STEPS

def get_temperature(step):
    """Calculates temperature based on the simulation step."""
    if HEAT_START <= step <= HEAT_END:
        # Heating Phase
        progress = (step - HEAT_START) / RUN_STEPS
        return T_MIN + progress * (T_MAX - T_MIN)
    elif COOL_START <= step <= COOL_END:
        # Cooling Phase
        progress = (step - COOL_START) / RUN_STEPS
        return T_MAX - progress * (T_MAX - T_MIN)
    elif HEAT_END < step < COOL_START:
        # Holding Phase
        return T_MAX
    return None

# --- DATA PROCESSING ---
steps = []
temps = []
fcc_fractions = []

total_atoms = 4000  

print("Reading dump file...")

with open('dump.cu.lammpstrj', 'r') as f:
    while True:
        line = f.readline()
        if not line: break
        
        if "ITEM: TIMESTEP" in line:
            step = int(f.readline())
            
            
            for _ in range(7):
                f.readline()
            
            
            fcc_count = 0
            for _ in range(total_atoms):
                atom_line = f.readline().split()
                
                cna_type = float(atom_line[-1])
                if cna_type == 1.0:
                    fcc_count += 1
            
            
            T = get_temperature(step)
            if T is not None:
                steps.append(step)
                temps.append(T)
                fcc_fractions.append(fcc_count / total_atoms)


plt.figure(figsize=(10, 6))


plt.plot(temps, fcc_fractions, 'g.', markersize=2, alpha=0.6, label='FCC Fraction')



plt.xlabel('Temperature (K)', fontsize=14)
plt.ylabel('Fraction of FCC Atoms', fontsize=14)
plt.title('Structural Evolution of Copper (FCC Fraction)', fontsize=16)
plt.ylim(-0.05, 1.05)
plt.grid(True, alpha=0.3)
plt.legend()


plt.savefig('copper_fcc_fraction.png', dpi=300)
print("Plot saved as copper_fcc_fraction.png")
plt.show()