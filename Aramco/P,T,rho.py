import matplotlib.pyplot as plt
import numpy as np
from thermopack.cubic import cubic
from thermopack.tcPR import tcPR
# Create cubic mixture object
mix = cubic("CO2,N2", "PR")

# Define composition and initial conditions
z = np.array([0.85, 0.15])
T_range = np.linspace(-10, 25, 500) + 273.15 # Temperature range from -10 to 25 degrees Celsius
p_range = np.linspace(4e6, 10e6, 500)  # Pressure range from 4 to 7 MPa

Rhog = np.zeros((len(T_range), len(p_range)))
Rhol = np.zeros((len(T_range), len(p_range)))

# Calculate densities
for i, T in enumerate(T_range):
    for j, p in enumerate(p_range):
        flsh = mix.two_phase_tpflash(T, p, z)
        
        vg, = mix.specific_volume(T, p, z, mix.VAPPH)
        vl, = mix.specific_volume(T, p, z, mix.LIQPH)

        M1 = 44.01  # Molar mass of CO2
        M2 = 28.02  # Molar mass of N2


        if flsh.phase == 4 and mix.guess_phase(T,p,z) == 2:    
            Mg = z[0] * M1 + z[1] * M2
            Ml = 0
        elif flsh.phase == 4 and mix.guess_phase(T,p,z) == 1:    
            Mg = 0
            Ml = flsh.x[0] * M1 + flsh.x[1] * M2
        elif flsh.phase == 0:
            Mg = flsh.y[0] * M1 + flsh.y[1] * M2
            Ml = flsh.x[0] * M1 + flsh.x[1] * M2

        rhog = Mg / vg
        rhol = Ml / vl

        Rhog[i, j] = rhog / 1e3  # Convert to kg/m^3
        Rhol[i, j] = rhol / 1e3  # Convert to kg/m^3



# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Plot for Rhog (gas phase)
im1 = ax1.imshow(Rhog, extent=[p_range[0], p_range[-1], T_range[-1], T_range[0]],
                 aspect='auto')
ax1.set_xlabel('Pressure (Pa)')
ax1.set_ylabel('Temperature (°C)')
ax1.set_title('Gas Phase Density (kg/m^3)')
fig.colorbar(im1, ax=ax1, label='Density (kg/m^3)')


T1, P1 = mix.get_envelope_twophase(1.0e4, z, maximum_pressure=1.5e7)
ax1.plot( P1 , T1, color='r')
ax1.set_xlim(4e6,10e6)
ax1.set_ylim(273.15-10,273.15+25)


# Plot for Rhol (liquid phase)
im2 = ax2.imshow(Rhol, extent=[p_range[0], p_range[-1], T_range[-1], T_range[0]],
                 aspect='auto')
ax2.set_xlabel('Pressure (Pa)')
ax2.set_ylabel('Temperature (°C)')
ax2.set_title('Liquid Phase Density (kg/m^3)')
fig.colorbar(im2, ax=ax2, label='Density (kg/m^3)')

ax2.plot( P1 , T1, color='r')
ax2.set_xlim(4e6,10e6)
ax2.set_ylim(273.15-10,273.15+25)

plt.suptitle('Density of the mixture CO2/N2 as a function of T and P')
plt.tight_layout()
plt.show()


np.savetxt('Rhog.txt', Rhog)
np.savetxt('Rhol.txt', Rhol)