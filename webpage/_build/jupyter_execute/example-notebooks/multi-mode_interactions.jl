using QuantumOptics
using IonSim
using DSP: periodogram
using LaTeXStrings
import PyPlot
const plt = PyPlot

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

# Setup system
C = Ca40(["S-1/2","D-1/2"])
L = Laser()
νᵣ = √14 * 1e6  # radial trap frequency
νₐ = 1.5e6  # axial trap frequency
chain = LinearChain(
        ions=[C], com_frequencies=(x=νᵣ, y=νᵣ, z=νₐ), 
        vibrational_modes=(x=[1], z=[1])
    )
T = Trap(configuration=chain, B=4e-4, Bhat=ẑ, δB=0, lasers=[L])
# set laser parameters
L.k = (x̂ + ẑ)/√2
L.ϵ = ẑ
# set carrier transition Rabi frequency to 500 kHz
Efield_from_rabi_frequency!(5e5, T, 1, 1, ("S-1/2","D-1/2"));

axial = T.configuration.vibrational_modes.z[1]
radial = T.configuration.vibrational_modes.x[1]

# set vibrational mode Hilbert space dimensions
axial.N = 8; radial.N = 8  

ρᵢ_ion = dm(C["S-1/2"])
ρᵢ_axial = thermalstate(axial, 0.5)
ρᵢ_radial = thermalstate(radial, 0.5)
# Set initial state to ρᵢ = |↓, n̄ᵣ=0.5, n̄ₐ=0.5⟩
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial;  

tspan = 0:10:400
fout(t, ρ) = real(expect(ionprojector(T, "D-1/2"), ρ))
J = (-dm(C["S-1/2"]) + dm(C["D-1/2"])) ⊗ one(radial) ⊗ one(axial)  # Collapse operator
γ = 1e4 * 1e-6 
ωc = transition_frequency(T, 1, ("S-1/2", "D-1/2"))  # Carrier transition frequency

# Scan laser detuning
Δlist = range(ωc + νₐ, stop=ωc + 2νᵣ + 1e5, length=300)
exclist = []
for Δ in Δlist
    L.Δ = Δ
    h = hamiltonian(T, rwa_cutoff=1.5e6, lamb_dicke_order=2)
    _, sol = timeevolution.master_dynamic(tspan, ρᵢ, (t, ρ) -> (h(t, ρ), [J], [J], [γ]), fout=fout)
    push!(exclist, real(sol[end]))
end 

plt.plot(@.((Δlist - ωc) / 1e6) , exclist)
plt.text(νₐ/1e6 + 0.1, 0.45, L"νₐ", fontsize=15)
plt.text(νᵣ/1e6 + 0.1, 0.45, L"νᵣ", fontsize=15)
plt.text(2νₐ/1e6, 0.2, L"2⋅νₐ", fontsize=15)
plt.text(2νᵣ/1e6 - 0.3, 0.125, L"2⋅νᵣ", fontsize=15)
plt.text((νᵣ + νₐ)/1e6 + 0.1, 0.3, L"νᵣ + νₐ", fontsize=15)
plt.text((νᵣ - νₐ)/1e6 - 0.1, 0.05, L"νᵣ - νₐ", fontsize=15)
plt.xlim((Δlist[1] - ωc) / 1e6, (Δlist[end] - ωc) / 1e6 + 0.01)
plt.ylim(-0.01, 0.5)
plt.xlabel("Detuning relative to carrier (MHz)")
plt.ylabel("Excitation")
plt.show()

# set carrier transition Rabi frequency to 1 MHz
Efield_from_rabi_frequency!(1e6, T, 1, 1, ("S-1/2","D-1/2"))
tspan = 0:0.25:800
L.Δ = ωc + axial.ν
axial.N = 10
radial.N = 10
h = hamiltonian(T, rwa_cutoff=1e4, displacement="truncated", time_dependent_eta=false)

# Set initial state to |↓, n̄ᵣ=0, n̄ₐ=1.5⟩ and solve
ρᵢ_axial = thermalstate(axial, 1.5)
ρᵢ_radial = thermalstate(radial, 0)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_ax = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

# Set initial state to |↓, n̄ᵣ=1.5, n̄ₐ=0⟩ and solve
ρᵢ_axial = thermalstate(axial, 0)
ρᵢ_radial = thermalstate(radial, 1.5)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_rad = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout);

plt.plot(tout, sol_ax, label="n̄_axial=1.5, n̄_radial=0")
plt.plot(tout, analytical.rabi_flop.(tout, 1, get_η(axial, L, C), 1.5, s=1), ls="dotted", label="n̄_axial=1.5, n̄_radial=0 (theoretical)")
plt.plot(tout, sol_rad, alpha=0.5, label="n̄_axial=0, n̄_radial=1.5", zorder=-1)
plt.legend(loc=4, framealpha=1)
plt.ylim(0, 1)
plt.xlim(0, tout[end])
plt.xlabel("Time (μs)")
plt.show()

tspan = 0:0.25:800
L.Δ = ωc + radial.ν
h = hamiltonian(T, rwa_cutoff=1e4)

# Set initial state to |↓, n̄ᵣ=0, n̄ₐ=1.5⟩ and solve
ρᵢ_axial = thermalstate(axial, 1.5)
ρᵢ_radial = thermalstate(radial, 0)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_ax = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

# Set initial state to |↓, n̄ᵣ=1.5, n̄ₐ=0⟩ and solve
ρᵢ_axial = thermalstate(axial, 0)
ρᵢ_radial = thermalstate(radial, 1.5)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_rad = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout);

plt.plot(tout, sol_rad, label="n̄_axial=0, n̄_radial=1.5")
plt.plot(tout, analytical.rabi_flop.(tout, 1, get_η(radial, L, C), 1.5, s=1), ls="dotted", label="n̄_axial=0, n̄_radial=1.5 (theoretical)")
plt.plot(tout, sol_ax, alpha=0.5, label="n̄_axial=1.5, n̄_radial=0", zorder=-1)
plt.legend(loc=4, framealpha=1)
plt.ylim(0, 1)
plt.xlim(0, tout[end])
plt.xlabel("Time (μs)")
plt.show()

fout(t, ρ) = real(expect(ionprojector(T, "D-1/2"), ρ));

tspan = 0:1:10000
axial.N = 15
radial.N = 15
L.Δ = ωc + radial.ν + axial.ν
h = hamiltonian(T, rwa_cutoff=1e4)

# Set initial state to |↓, n̄ᵣ=0, n̄ₐ=1.5⟩
ρᵢ_axial = thermalstate(axial, 1.5)
ρᵢ_radial = thermalstate(radial, 0)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_ax = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

# Set initial state to |↓, n̄ᵣ=1.5, n̄ₐ=0⟩ and solve
ρᵢ_axial = thermalstate(axial, 0)
ρᵢ_radial = thermalstate(radial, 1.5)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_rad = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

# Set initial state to |↓, n̄ᵣ=1.5, n̄ₐ=1.5⟩ and solve
ρᵢ_axial = thermalstate(axial, 1.5)
ρᵢ_radial = thermalstate(radial, 1.5)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  
tout, sol_ax_rad = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout);

plt.plot(tout, sol_rad, label="n̄_axial=0, n̄_radial=1.5")
plt.plot(tout, sol_ax, label="n̄_axial=1.5, n̄_radial=0")
plt.plot(tout, sol_ax_rad, label="n̄_axial=1.5, n̄_radial=1.5")
plt.legend(loc=1, framealpha=1)
plt.ylim(0, 1)
plt.xlim(0, tout[10000])
plt.xlabel("Time (μs)")
plt.show()

ηᵣ = get_η(radial, L, C)
ηₐ = get_η(axial, L, C)
η2 = ηᵣ * ηₐ
P(n1, n2; n̄₁=1.5, n̄₂=1.5) =  n̄₁^n1 * n̄₂^n2 / ((n̄₁ + 1)^(n1 + 1) * (n̄₂ + 1)^(n2 + 1))
sample_rate = (length(tout) - 1) / tout[end]
pwarray1 = periodogram(convert(Array{Float64}, sol_ax_rad), fs=sample_rate)
p1 = sqrt.(pwarray1.power[2:length(pwarray1.power) ÷ 2])
f1 = pwarray1.freq[2:length(pwarray1.power) ÷ 2]
pwarray2 = periodogram(convert(Array{Float64}, sol_rad), fs=sample_rate)
p2 = sqrt.(pwarray2.power[2:length(pwarray2.power) ÷ 2])
f2 = pwarray2.freq[2:length(pwarray2.power) ÷ 2]

fig, ax = plt.subplots()
plt.plot(f2, p2, label=L"n̄ₐ=0, n̄ᵣ=1.5 (simulation)")
plt.plot(f1, p1, label=L"n̄ₐ=1.5, n̄ᵣ=1.5 (simulation)")
bins = [η2 * √i for i in 1:9]

nlist = []
for n in 1:9
    templist = []
    for i in 0:9, j in 0:9
        if (i + 1) * (j + 1) == n
            push!(templist, (i, j))
        end
    end
    push!(nlist, templist)
end
double_thermal = [sum([P(n...) for n in nl]) for nl in nlist]
single_thermal = [P(n, 0) for n in 0:8]
plt.bar(
        bins,  double_thermal .* maximum(p1)/double_thermal[2], 
        width=0.0002, alpha=0.25, color="C1", 
        label=L"n̄ₐ=1.5, n̄ᵣ=1.5 (theory)"
    )
plt.bar(
        bins, single_thermal .* maximum(p2)/single_thermal[1], 
        width=0.0002, alpha=0.25, color="C0", 
        label=L"n̄ₐ=0, n̄ᵣ=1.5 (theory)"
    )
ax.set_xticks(bins)
ax.set_xticklabels([L"ηᵣ⋅ηₐ", L"√2⋅ηᵣ⋅ηₐ", L"√3⋅ηᵣ⋅ηₐ", L"2⋅ηᵣ⋅ηₐ", L"√5⋅ηᵣ⋅ηₐ", L"√6⋅ηᵣ⋅ηₐ", L"√7⋅ηᵣ⋅ηₐ", L"√8⋅ηᵣ⋅ηₐ", L"3⋅ηᵣ⋅ηₐ"], rotation=40)
ax.set_yticklabels([])
plt.xlim(0.0018, 0.006)
plt.xlabel("Frequency [2π/Ω_carrier]")
plt.ylabel("Spectral Weights [abu]")
plt.legend(loc=1)
plt.show()

tspan = 0:5:1000
ρᵢ_axial = thermalstate(axial, 0)
ρᵢ_radial = thermalstate(radial, 0)
ρᵢ = ρᵢ_ion ⊗ ρᵢ_radial ⊗ ρᵢ_axial  # Set initial state to |↓, n̄ᵣ=0, n̄ₐ=0⟩

L.Δ = ωc + 2radial.ν  # Set detuning to second order radial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[1, 2])  # First order LD approx on radial mode and second order on axial mode
tout, sol_rad = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

L.Δ = ωc + 2axial.ν  # Set detuning to second order axial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[1, 2])
tout, sol_ax = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

L.Δ = ωc + radial.ν + axial.ν  # Set detuning to axial + radial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[1, 2])
tout, sol_rad_ax = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

L.Δ = ωc + 2radial.ν  # Set detuning to second order radial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[2, 1])  # First order LD approx on radial mode and second order on axial mode
tout, sol_radr = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

L.Δ = ωc + 2axial.ν  # Set detuning to second order axial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[2, 1])
tout, sol_axr = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout)

L.Δ = ωc + radial.ν + axial.ν  # Set detuning to axial + radial sideband
h = hamiltonian(T, rwa_cutoff=1e4, lamb_dicke_order=[2, 1])
tout, sol_rad_axr = timeevolution.schroedinger_dynamic(tspan, ρᵢ, h, fout=fout);

fig, ax = plt.subplots(1, 2, sharey=true, figsize=(15, 4))
ax[1].plot(tout, sol_rad, label=L"Δ = ωc + 2⋅νᵣ")
ax[1].plot(tout, sol_ax, label=L"Δ = ωc + 2⋅νₐ")
ax[1].plot(tout, sol_rad_ax, label=L"Δ = ωc + νᵣ + νₐ")
ax[1].legend(loc=1, framealpha=1)
ax[1].set_title("lamb_dicke_order = [1, 2]")
ax[2].plot(tout, sol_radr, label=L"Δ = ωc + 2⋅νᵣ")
ax[2].plot(tout, sol_axr, label=L"Δ = ωc + 2⋅νₐ")
ax[2].plot(tout, sol_rad_axr, label=L"Δ = ωc + νᵣ + νₐ")
ax[2].set_title("lamb_dicke_order = [2, 1]")
ax[2].legend(loc=1, framealpha=1)
plt.ylim(-0.05, 1)
ax[1].set_xlim(0, tout[end]); ax[2].set_xlim(0, tout[end])
ax[1].set_xlabel("Time [μs]"); ax[2].set_xlabel("Time [μs]")
plt.show()
