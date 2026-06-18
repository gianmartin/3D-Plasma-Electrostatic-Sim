import pandas as pd
import matplotlib.pyplot as plt

def plot_diagnostics(filename='runtime_diags.csv'):
    # 1. Load the data
    try:
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return

    # 2. Setup the figure
    # We will create 2 subplots: one for Particles, one for Energy
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # --- PLOT 1: Particle Conservation ---
    # We plot the macroparticle counts. If these drop, particles are leaving the box.
    # If they fluctuate up/down, you have a race condition in the particle list management.
    ax1.plot(df['time'], df['mp_count.O+'], label='O+ Ions', color='blue', linewidth=2)
    ax1.plot(df['time'], df['mp_count.e-'], label='Electrons', color='orange', linewidth=2)
    
    ax1.set_ylabel('Macroparticle Count')
    ax1.set_title('Particle Conservation Check')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.6)

    # Add a text annotation if counts are perfectly stable
    o_variation = df['mp_count.O+'].max() - df['mp_count.O+'].min()
    e_variation = df['mp_count.e-'].max() - df['mp_count.e-'].min()
    if o_variation == 0 and e_variation == 0:
        ax1.text(0.5, 0.5, "PERFECT CONSERVATION", transform=ax1.transAxes, 
                 ha='center', fontsize=12, color='green', fontweight='bold', 
                 bbox=dict(facecolor='white', alpha=0.8))

    # --- PLOT 2: Energy Conservation ---
    # Total KE is the sum of both species
    total_ke = df['KE.O+'] + df['KE.e-']
    
    ax2.plot(df['time'], total_ke, label='Total Kinetic Energy', color='green', linestyle='--')
    ax2.plot(df['time'], df['PE'], label='Potential Energy', color='red', linestyle='--')
    ax2.plot(df['time'], df['E_total'], label='Total System Energy', color='black', linewidth=2.5)

    ax2.set_ylabel('Energy (J)')
    ax2.set_xlabel('Time (s)')
    ax2.set_title('System Energy Stability')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.6)

    # Calculate Energy Drift (drift is the enemy of PIC codes)
    start_energy = df['E_total'].iloc[0]
    end_energy = df['E_total'].iloc[-1]
    drift_pct = ((end_energy - start_energy) / start_energy) * 100
    
    print(f"Total Energy Drift: {drift_pct:.4f}%")
    
    # Save the plot
    plt.tight_layout()
    plt.savefig('simulation_diagnostics.png')
    print("Saved plot to simulation_diagnostics.png")

if __name__ == "__main__":
    plot_diagnostics()