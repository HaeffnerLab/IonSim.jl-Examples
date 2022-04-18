using QuantumOptics
using IonSim
import PyPlot
const plt = PyPlot;

# set some plot configs
plt.matplotlib.rc("xtick", top=false)
plt.matplotlib.rc("ytick", right=false, left=false)
plt.matplotlib.rc("axes", labelsize=20, titlesize=20, grid=true)
plt.matplotlib.rc("axes", linewidth=2)
plt.matplotlib.rc("grid", alpha=0.25, linestyle="--")
plt.matplotlib.rc("font", family="Palatino", weight="medium")
plt.matplotlib.rc("figure", figsize=(8,4))
plt.matplotlib.rc("xtick.major", width=2)
plt.matplotlib.rc("ytick.major", width=2)

C = Ca40(["S-1/2", "D-1/2"])
L1 = Laser(pointing=[(1, 1.), (2, 1.)])
L2 = Laser(pointing=[(1, 1.), (2, 1.)])
chain = LinearChain(
        ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), 
        vibrational_modes=(x=[1], z=[1])
    )
T = Trap(configuration=chain, B=6e-4, Bhat=(x̂ + ẑ)/√2, lasers=[L1, L2])

axial_mode = T.configuration.vibrational_modes.z[1]
radial_mode = T.configuration.vibrational_modes.x[1]

# Set Hilbert space dimension for participating vibrational modes
axial_mode.N = 5
radial_mode.N = 5;

Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
ϵ = 40e3

L1.Δ = Δf + axial_mode.ν + ϵ 
L1.k = ẑ
L1.ϵ = x̂

L2.Δ = Δf - axial_mode.ν - ϵ
L2.k = ẑ
L2.ϵ = x̂

η = abs(get_η(axial_mode, L1, C))
Ω = √(1e3 * ϵ) / η  # This will give a 1kHz MS strength, since coupling goes like (ηΩ)^2/ϵ

Efield_from_rabi_frequency!(Ω, T, 1, 1, ("S-1/2", "D-1/2"))
Efield_from_rabi_frequency!(Ω, T, 2, 1, ("S-1/2", "D-1/2"));

h = hamiltonian(T, lamb_dicke_order=1, rwa_cutoff=1e6);
tout, sol_res = timeevolution.schroedinger_dynamic(
        0:2:3000, C["D-1/2"] ⊗ C["S-1/2"] ⊗ axial_mode[0] ⊗ radial_mode[0], h
    );

h = hamiltonian(T, lamb_dicke_order=1, rwa_cutoff=1e6);

@time tout, sol_res = timeevolution.schroedinger_dynamic(
        0:2:3000, C["D-1/2"] ⊗ C["S-1/2"] ⊗ axial_mode[0] ⊗ radial_mode[0], h
    );

SD_res = expect(ionprojector(T, "S-1/2", "D-1/2"), sol_res)
plt.plot(tout, SD_res, label="resonant")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("SD")
plt.xlabel("Time (μs)");

# We'll add the barrier term by setting a nonzero value to 
# the trap's B-field gradient:
set_gradient!(T, (1, 2), ("S-1/2", "D-1/2"), √5.25 * 1e3); 
# ΩMS = 1 kHz, so setting Δ = √5.25 kHz gives Ω' = √(ΩMS^2 + Δ^2) = 2.5 kHz

h = hamiltonian(T, lamb_dicke_order=1, rwa_cutoff=1e6);

@time tout, sol_det = timeevolution.schroedinger_dynamic(
        0:2:3000, 
        C["D-1/2"] ⊗ C["S-1/2"] ⊗ axial_mode[0] ⊗ radial_mode[0], h
    );

SD_det = expect(ionprojector(T, "S-1/2", "D-1/2"), sol_det)
plt.plot(tout, SD_res, label="resonant")
plt.plot(tout, SD_det, label="detuned")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("SD")
plt.xlabel("Time (μs)");

L3 = Laser(pointing=[(2, 1.)])
L4 = Laser(pointing=[(2, 1.)])
T.lasers = [L1, L2, L3, L4];

axial_mode = T.configuration.vibrational_modes.z[1]
radial_mode = T.configuration.vibrational_modes.x[1]

Δf = transition_frequency(T, 2, ("S-1/2", "D-1/2"))
νeff = √5.25 * 1e3

# We set the effective vibrational frequency by detuning 
# the 2-photon transition off resonance
L3.Δ = Δf + radial_mode.ν / 2 + νeff  
L3.k = x̂
L3.ϵ = ẑ

L4.Δ = Δf - radial_mode.ν / 2
L4.k = x̂
L4.ϵ = ẑ

η = abs(get_η(radial_mode, L3, C))
Ω = √((0.5e3 * radial_mode.ν) / η)  # Set κ = 0.5 kHz, since κ = ηΩ^2/ν

Efield_from_rabi_frequency!(Ω, T, 3, 2, ("S-1/2", "D-1/2"))
Efield_from_rabi_frequency!(Ω, T, 4, 2, ("S-1/2", "D-1/2"));

h = hamiltonian(T, lamb_dicke_order=1, rwa_cutoff=1.51e6);

@time tout, sol_vaet = timeevolution.schroedinger_dynamic(
        0:2:3000, 
        C["D-1/2"] ⊗ C["S-1/2"] ⊗ axial_mode[0] ⊗ radial_mode[0], h
    );

SD_vaet = expect(ionprojector(T, "S-1/2", "D-1/2"), sol_vaet)
plt.plot(tout, SD_res, label="resonant")
plt.plot(tout, SD_det, label="detuned")
plt.plot(tout, SD_vaet, label="VAET")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.xlabel("Time (μs)");
