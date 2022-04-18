using IonSim
using QuantumOptics: timeevolution, stochastic

using Pkg
Pkg.add("StochasticDiffEq")
Pkg.add("DSP")

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
print(typeof(C) <: Ion)

?Ca40

L = Laser()
print(L)

chain = LinearChain(
        ions=[C], com_frequencies=(x=3e6,y=3e6,z=1e6), 
        vibrational_modes=(;z=[1])
    )
print(typeof(chain) <: IonConfiguration)

chain.vibrational_modes

T = Trap(configuration=chain, B=4e-4, Bhat=ẑ, δB=0, lasers=[L]);

L.k = (x̂ + ẑ)/√2 
L.ϵ = (x̂ - ẑ)/√2;

Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L.Δ = Δf

?transition_frequency

Efield_from_pi_time!(2e-6, T, 1, 1, ("S-1/2", "D-1/2"));  # Sets pi_time to 2 μs

?Efield_from_pi_time

?hamiltonian

h = hamiltonian(T, timescale=1e-6);

tspan = 0:0.01:10
mode = T.configuration.vibrational_modes.z[1]
tout, sol = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ mode[0], h);

tspan = 0:0.01:10
mode = T.configuration.vibrational_modes.z[1]
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ mode[0], h);

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

L.Δ = Δf + 2.5e5

h = hamiltonian(T, timescale=1e-6)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

L.Δ = Δf

C.stark_shift["S-1/2"] = -1.25e5
C.stark_shift["D-1/2"] = 1.25e5

h = hamiltonian(T, timescale=1e-6)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

ψi_ion = dm(C["S-1/2"]) 
mode.N = 100
ψi_mode = thermalstate(mode, 10)
ψi = ψi_ion ⊗ ψi_mode
h = hamiltonian(T, timescale=1e-6, rwa_cutoff=Inf, time_dependent_eta=false, displacement="truncated")
tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);

# reset all artificial Stark shifts to zero
zero_stark_shift(C)  

ψi_ion = dm(C["S-1/2"]) 

# we'd like to look at a pretty hot ion, so we need to increase our mode dimension, 
# which is set to N=10 by default here we set it to 100
mode.N = 100

ψi_mode = thermalstate(mode, 10)

ψi = ψi_ion ⊗ ψi_mode

tspan = 0:0.1:50
h = hamiltonian(T, timescale=1e-6, rwa_cutoff=Inf, time_dependent_eta=false, displacement="truncated")
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

?analytical.rabi_flop

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex, label="numerical")

η = get_η(mode, L, C)
plt.plot(
        tout, analytical.rabi_flop(tout, 1/4, η, 10), 
        linestyle="--", label="analytical"
    )
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

h = hamiltonian(T, lamb_dicke_order=0)  # set lamb_dicke_order here
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex, label="numerical")
η = get_η(mode, L, C)
plt.plot(
        tout, analytical.rabi_flop(tout, 1/4, η, 10), 
        linestyle="--", label="analytical"
    )
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

tspan = 0:0.1:100
# set the motional dimension back to 10
mode.N = 10  

# tune laser frequency to blue sideband
L.Δ  = Δf + mode.ν  

h = hamiltonian(T)
@time tout, sol_blue = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

# tune laser frequency to red sideband
L.Δ  = Δf - mode.ν  

h = hamiltonian(T)
@time tout, sol_red = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

ex_blue = expect(ionprojector(T, "D-1/2"), sol_blue)
ex_red = expect(ionprojector(T, "D-1/2"), sol_red)
plt.plot(tout, ex_blue, color="blue", label="blue sideband")
plt.plot(tout, ex_red, color="red", label="red sideband")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

L.Δ  = Δf + mode.ν  # tune laser frequency to blue sideband

h = hamiltonian(T)
@time tout, sol_blue = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ fockstate(mode, 1), h)

L.Δ  = Δf - mode.ν  # tune laser frequency to red sideband

h = hamiltonian(T, rwa_cutoff=1e7)
@time tout, sol_red = timeevolution.schroedinger_dynamic(tspan, ionstate(T, "S-1/2") ⊗ fockstate(mode, 1), h)

ex_blue = expect(ionprojector(T, "D-1/2"), sol_blue)
ex_red = expect(ionprojector(T, "D-1/2"), sol_red)
step = 1
plt.plot(tout[1:step:end], ex_blue[1:step:end], color="blue", label="blue sideband")
plt.plot(tout[1:step:end], ex_red[1:step:end], color="red", label="red sideband")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time [μs]");

# tune laser frequency to blue sideband
L.Δ  = Δf + mode.ν  

h = hamiltonian(T, rwa_cutoff=1e-5) # set rwa_cutoff here
@time tout, sol_blue = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[1], h)

# tune laser frequency to red sideband
L.Δ  = Δf - mode.ν  

h = hamiltonian(T, rwa_cutoff=1e-5) # set rwa_cutoff here
@time tout, sol_red = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[1], h)

ex_blue = expect(ionprojector(T, "D-1/2"), sol_blue)
ex_red = expect(ionprojector(T, "D-1/2"), sol_red)
plt.plot(tout, ex_blue, color="blue", label="blue sideband")
plt.plot(tout, ex_red, color="red", label="red sideband")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

# add additional detuning to BSB to compensate for carrier Stark shift
L.Δ  = Δf + mode.ν - 31e3

tspan = 0:0.1:60
h = hamiltonian(T)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

# get previous value of electric field strength
E = Efield_from_pi_time(2e-6, T, 1, 1, ("S-1/2", "D-1/2")) 

# Simple amplitude ramping function
function Ω(t)
    if t < 4
        return E * sin(2π * t / 16)^2
    end
    E
end

L.E = Ω;

# add additional detuning to BSB to compensate for carrier Stark shift
L.Δ  = Δf + mode.ν - 31e3

tspan = 0:0.1:60
h = hamiltonian(T)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.plot(tout, @.(Ω(tout) / 2E), linestyle="--", label="scaled ramp")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

# To make things interesting let's also reduce the ramp time from 4 μs to 1 μs
L.E = t -> t < 1 ?  E * sin(2π * t / 4)^2 : E;

# add additional detuning to BSB to compensate for carrier Stark shift
L.Δ  = Δf + mode.ν - 31e3

tspan = 0:0.1:60
h = hamiltonian(T)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex)
plt.plot(tout, @.(Ω(tout) / 2E), linestyle="--", label="scaled ramp")
plt.xlim(tout[1], tout[end])
plt.legend(loc=1)
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

# the length of time of the frequency chirp in μs
Tp = 150

# Δϕ is equal to half the detuning range we will chirp the laser over multiplied by the timescale (1e-6)
Δϕ = 2π * 150e-3

# Set the frequency chirp
L.ϕ = t -> (-Δϕ + (2Δϕ / Tp) * t) * t;

E = Efield_from_pi_time(6.5e-6, T, 1, 1, ("S-1/2", "D-1/2"))
tr = 33
function Ω(t)
    if t < tr
        return E * sin(2π * t / 4tr)^2
    elseif tr <= t <= 150 - tr
        return E
    elseif 150 - tr < t < 150
        return E * sin(2π * (t - 150) / 4tr)^2
    else
        return 0
    end
end
L.E = Ω;

# Set the B-field to match the value in the reference
T.B = 2.9e-4
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))

L.Δ  = Δf  # set detuning back to carrier 

tspan = 0:0.1:175
h = hamiltonian(T)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)

ex = expect(ionprojector(T, "D-1/2"), sol)
plt.plot(tout, ex, lw=3)
plt.plot(
        tout, @.(L.E(tout) / 2E), 
        linestyle="--", label="scaled amplitude profile"
    )
plt.plot(
        tout, @.(L.ϕ(tout) / (2Δϕ * tout)), 
        linestyle="--", label="scaled frequency profile"
    )
plt.xlim(tout[1], tout[end])
plt.legend(loc=1)
plt.ylim(0, 1)
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

detuning_range = 2π*1e-6 .* collect(1e3:2e4:2e6)
pops = Vector{Float64}(undef, 0)
for Δϕ in detuning_range
    L.ϕ = t -> (-Δϕ/2 + (Δϕ / Tp) * t) * t
    h = hamiltonian(T)
    _, sol = timeevolution.schroedinger_dynamic(tspan, C["S-1/2"] ⊗ mode[0], h)
    push!(pops, real(expect(ionprojector(T, "D-1/2"), sol)[end]))
end
plt.plot(@.(detuning_range * 1e3), pops, lw=3, color="C1")
plt.xlim(-10, detuning_range[end] * 1e3)
plt.ylim(0, 1.01)
plt.ylabel("Final Excited State\nPopulation")
plt.xlabel("detuning range (kHz)");

# get rid of frequency chirp
L.ϕ = 0
# and recompute the transition frequency
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L.Δ = Δf;

using StochasticDiffEq

tspan = collect(0:10.:30000)
β = 1e-2
σ = √(2β)
w = StochasticDiffEq.OrnsteinUhlenbeckProcess(β, 0.0, σ, 0.0, 0.0)
w.dt = 0.1;

plt.title("O-U noise")
warray = [w(t)[1] for t in tspan]
plt.plot(tspan, warray)
plt.show()

using DSP: periodogram

sample_rate = (length(tspan) - 1)/tspan[end]
pwarray = periodogram(warray, fs=sample_rate)
freqs = pwarray.freq
abs_freqs = freqs .* 1e6
powers = pwarray.power
N = 50
for i in 1:N-1
    w = StochasticDiffEq.OrnsteinUhlenbeckProcess(β, 0.0, σ, 0.0, 0.0)
    w.dt = 0.1
    warray = [w(t)[1] for t in tspan]
    powers .+= periodogram(warray, fs=sample_rate).power
end
powers ./= powers[2]

plt.loglog(abs_freqs, powers, label="numerical")
ou_psd(f) = σ^2 / (β^2 + (2π * f)^2)
plt.loglog(abs_freqs, ou_psd.(freqs) ./ ou_psd(freqs[2]), color="C1", label="theoretical", ls="--")
plt.ylabel("relative power")
plt.xlabel("frequency [Hz]")
plt.title("O-U PSD")
plt.ylim(10e-3, 2)
plt.xlim(250, 2e4)
plt.legend()
plt.show()

# we're not paying attention to the vibrational mode here, so we set its dimension to 1,
# effectively ignoring it
mode.N = 1

# construct a zero operator
L.E = 0
T.δB = 0.1
h = hamiltonian(T)

# let's work with a Δm=2 transition
T.configuration.ions[1].selected_level_structure = ["S-1/2", "D-5/2"]

# construct noise operator
T.δB = 5e-1
hs = hamiltonian(T)
hsvec = (t,ψ) -> [hs(t, ψ)]

ψi_ion = (C["S-1/2"] + C["D-5/2"])/√2
ψi = ψi_ion ⊗ mode[0]

Ntraj = 100
ex = zero(tspan)

# iterate SDE solver Ntraj times and avergage results
w = StochasticDiffEq.OrnsteinUhlenbeckProcess(β, 0.0, σ, 0.0, 0.0)
for i in 1:Ntraj
    tout, sol = stochastic.schroedinger_dynamic(tspan, ψi, h, hsvec, noise=w,
    normalize_state=true, dt=0.1)
    ex .+= real.(expect(dm(ψi_ion) ⊗ one(mode), sol)) ./ Ntraj
end

plt.plot(tspan, ex, color="C0", label="Simulated Contrast")
plt.xlim(tspan[1], tspan[end])
plt.ylim(0, 1)
plt.legend()
plt.ylabel("Ramsey Signal")
plt.xlabel("Time (μs)")
plt.show()

γ1 = hs(1.0, 0).data[1, 1]/40
γ2 = hs(1.0, 0).data[2, 2]/40
rates = 1/4π .* abs.([γ1, γ2])
hs1 = C["S-1/2"] ⊗ C["S-1/2"]' ⊗ one(mode)
hs2 = C["D-5/2"] ⊗ C["D-5/2"]' ⊗ one(mode)
tspan = collect(0:1.:30000)
tout, sol = timeevolution.master(tspan, dm(ψi), h(1.0, 0.0), [hs1, hs2], rates=rates);

γ1 = hs(1.0, 0).data[1, 1]/40
γ2 = hs(1.0, 0).data[2, 2]/40

rates = 1/4π .* abs.([γ1, γ2])

hs1 = C["S-1/2"] ⊗ C["S-1/2"]' ⊗ one(mode)
hs2 = C["D-5/2"] ⊗ C["D-5/2"]' ⊗ one(mode)

tspan = collect(0:1.:30000)
@time tout, sol = timeevolution.master(tspan, dm(ψi), h(1.0, 0.0), [hs1, hs2], rates=rates)

ex = expect(dm(ψi_ion) ⊗ one(mode), sol)
plt.plot(tspan, ex)
plt.xlim(tspan[1], tspan[end])
plt.ylim(0, 1)
plt.ylabel("Ramsey Contrast")
plt.xlabel("Time (μs)")
plt.show()

T.configuration.ions[1].selected_level_structure = ["S-1/2", "D-1/2"]
T.B = 4e-4
T.δB = 0

# recompute the transition frequency
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L.Δ = Δf

E = Efield_from_pi_time(2e-6, T, 1, 1, ("S-1/2", "D-1/2"))
tspan = 0:0.1:60

# average over Ntraj runs
Ntraj = 1000
δE = 0.025E
ex = zero(tspan)
ψi = C["S-1/2"] ⊗ mode[0]
@time begin
    for i in 1:Ntraj
        ΔE = δE * randn()
        L.E = E + ΔE 
        h = hamiltonian(T)
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        ex .+= expect(ionprojector(T, "D-1/2"), sol) ./ Ntraj
    end
end

# compute expected τ
hz_per_E = 1 / Efield_from_rabi_frequency(1, ẑ, L, C, ("S-1/2", "D-1/2"))
τ_us = 1e6 / (2π * δE * hz_per_E)

plt.plot(tspan, ex)
plt.plot(
        tspan, @.((1 + exp(-(tspan / (√2 * τ_us))^2))/2), 
        color="C1", ls="--", label="Gaussian Envelope"
    )
plt.plot(
        tspan, @.((1 - exp(-(tspan / (√2 * τ_us))^2))/2), 
        color="C1", ls="--"
    )
plt.xlim(tspan[1], tspan[end])
plt.ylim(0, 1)
plt.legend()
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

T.configuration.ions[1].selected_level_structure = ["S-1/2", "D-1/2"]
T.B = 4e-4
T.δB = 0

# recompute the transition frequency
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L.Δ = Δf

E = Efield_from_pi_time(2e-6, T, 1, 1, ("S-1/2", "D-1/2"))
tspan = 0:0.1:60

# average over Ntraj runs
Ntraj = 1000
δE = 0.025E
ex = zero(tspan)
ψi = C["S-1/2"] ⊗ mode[0]
@time begin
    for i in 1:Ntraj
        ΔE = δE * randn()
        L.E = E + ΔE 
        h = hamiltonian(T)
        tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h)
        ex .+= expect(ionprojector(T, "D-1/2"), sol) ./ Ntraj
    end
end

# compute expected τ
hz_per_E = 1 / Efield_from_rabi_frequency(1, ẑ, L, C, ("S-1/2", "D-1/2"))
τ_us = 1e6 / (2π * δE * hz_per_E)

plt.plot(tspan, ex)
plt.plot(
        tspan, @.((1 + exp(-(tspan / (√2 * τ_us))^2))/2), 
        color="C1", ls="--", label="Gaussian Envelope"
    )
plt.plot(
        tspan, @.((1 - exp(-(tspan / (√2 * τ_us))^2))/2), 
        color="C1", ls="--"
    )
plt.xlim(tspan[1], tspan[end])
plt.ylim(0, 1)
plt.legend()
plt.ylabel("Excitation")
plt.xlabel("Time (μs)");

# Construct the system
C = Ca40(["S-1/2", "D-1/2"])
L1 = Laser(); L2 = Laser() 
chain = LinearChain(
        ions=[C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), 
        vibrational_modes=(;z=[1])
    )
T = Trap(configuration=chain, B=4e-4, Bhat=(x̂ + ẑ)/√2, lasers=[L1, L2]);

# Set the laser parameters
ϵ = 40e3
d = 80  # corrects for AC stark shift from single-photon coupling to sidebands
mode = T.configuration.vibrational_modes.z[1]
Δf = transition_frequency(T, 1, ("S-1/2", "D-1/2"))
L1.Δ = Δf + mode.ν + ϵ - d
L1.k = ẑ
L1.ϵ = x̂

L2.Δ = Δf - mode.ν - ϵ + d
L2.k = ẑ
L2.ϵ = x̂

mode.N = 5
η = abs(get_η(mode, L1, C))
Ω1 = √(1e3 * ϵ) / η  # This will give a 1kHz MS strength, since coupling goes like (ηΩ)^2/ϵ

Efield_from_rabi_frequency!(Ω1, T, 1, 1, ("S-1/2", "D-1/2"))
Efield_from_rabi_frequency!(Ω1, T, 2, 1, ("S-1/2", "D-1/2"));

ψi = C["S-1/2"] ⊗ C["S-1/2"] ⊗ mode[0];  # initial state

h = hamiltonian(T, rwa_cutoff=5e5)
tspan = 0:0.25:1000
tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);

# setup and run the Hamiltonian
h = hamiltonian(T, rwa_cutoff=5e5)
tspan = 0:0.25:1000
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);

SS = expect(ionprojector(T, "S-1/2", "S-1/2"), sol)
DD = expect(ionprojector(T, "D-1/2", "D-1/2"), sol)
SD = expect(ionprojector(T, "S-1/2", "D-1/2"), sol)
DS = expect(ionprojector(T, "D-1/2", "S-1/2"), sol)
plt.plot(tout, SS, label="SS")
plt.plot(tout, DD, label="DD")
plt.plot(tout, SD, label="SD")
plt.plot(tout, DS, label="DS")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.xlabel("Time (μs)");

?analytical.two_ion_ms

ex = analytical.two_ion_ms(tspan, 1e-6Ω1, 1e-6mode.ν, 1e-6mode.ν + 1e-6ϵ, η, 0)
plt.plot(tspan, SS)
plt.plot(tspan, DD)
plt.plot(tspan, ex[1], ls="--")
plt.plot(tspan, ex[2], ls="--")
plt.ylim(0, 1);

mode.N = 15
ρi = dm(ionstate(T, "S-1/2", "S-1/2")) ⊗ thermalstate(mode, 2)  # thermal state
h = hamiltonian(T, rwa_cutoff=5e5)
tspan = collect(0:0.25:1000)
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ρi, h)

SS = expect(ionprojector(T, "S-1/2", "S-1/2"), sol)
DD = expect(ionprojector(T, "D-1/2", "D-1/2"), sol)
SD = expect(ionprojector(T, "S-1/2", "D-1/2"), sol)
DS = expect(ionprojector(T, "D-1/2", "S-1/2"), sol)
plt.plot(tout, SS, label="SS")
plt.plot(tout, DD, label="DD")
plt.plot(tout, SD, label="SD")
plt.plot(tout, DS, label="DS")
ex = analytical.two_ion_ms(tspan, 1e-6Ω1, 1e-6mode.ν, 1e-6mode.ν + 1e-6ϵ, η, 2)
plt.plot(tspan, ex[1], ls="--")
plt.plot(tspan, ex[2], ls="--")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.xlabel("Time (μs)");

chain = LinearChain(
        ions=[C, C, C, C, C, C], com_frequencies=(x=3e6,y=3e6,z=1e6), 
        vibrational_modes=(;z=[1])
    )
T = Trap(configuration=chain, B=4e-4, Bhat=(x̂ + ẑ)/√2, lasers=[L1, L2])
global_beam!(T, L1)  # set L1 to shine equally on all ions
global_beam!(T, L2)  # set L2 to shine equally on all ions
mode = T.configuration.vibrational_modes.z[1]
mode.N = 5
η = abs(get_η(mode, L1, C))
Ω2 = √(1e3 * ϵ) / η  
Efield_from_rabi_frequency!(Ω2, T, 1, 1, ("S-1/2", "D-1/2"))
Efield_from_rabi_frequency!(Ω2, T, 2, 1, ("S-1/2", "D-1/2"))
ψi = ionstate(T, repeat(["S-1/2"], 6)...) ⊗ mode[0]

h = hamiltonian(T, rwa_cutoff=5e5)
tspan = 0:0.25:2000
@time tout, sol = timeevolution.schroedinger_dynamic(tspan, ψi, h);

S = expect(ionprojector(T, repeat(["S-1/2"], 6)...), sol)
D = expect(ionprojector(T, repeat(["D-1/2"], 6)...), sol)
plt.plot(tout, S, label="SSSSSS")
plt.plot(tout, D, label="DDDDDD")
plt.xlim(tout[1], tout[end])
plt.ylim(0, 1)
plt.legend(loc=1)
plt.title("6 Ion Molmer Sorensen")
plt.xlabel("Time (μs)");

const pc = IonSim.PhysicalConstants;

pc.ħ

pc.c

println(pc.μB / pc.ħ)
println(pc.ħ + pc.c)
println(pc.α^pc.α)
