using ControlSystems
using Plots

s = tf("s")
G = 1 / (2*s^2 + 3*s + 1)   #Replace with the desired TF

sys_poles = poles(G)
sys_zeros = tzeros(G)
wn, zeta, _ = damp(G)
K_dc = dcgain(G)[1]

println("  Poles       : ", sys_poles)
println("  Zeros       : ", isempty(sys_zeros) ? "[ none ]" : sys_zeros)
println("  ωₙ  [rad/s] : ", round.(wn,   digits=6))
println("  ζ           : ", round.(zeta, digits=6))
println("  K  (DC gain): ", round(K_dc,  digits=6))

for i in eachindex(wn)
    ζᵢ  = zeta[i]
    ωₙᵢ = wn[i]
    cls = ζᵢ < 0  ? "UNSTABLE"          :
          ζᵢ == 0 ? "UNDAMPED"          :
          ζᵢ < 1  ? "UNDERDAMPED"       :
          ζᵢ == 1 ? "CRITICALLY DAMPED" :
                    "OVERDAMPED"
    ωd_val = ζᵢ < 1 ? round(ωₙᵢ * sqrt(1 - ζᵢ^2), digits=6) : NaN
    println("  Pole[$i] class: $cls  |  ωd = $ωd_val rad/s")
end

p1 = pzmap(G,
    title  = "POLE-ZERO MAP  |  ζ=$(round.(zeta,digits=3))  K=$(round(K_dc,digits=3))",
    xlabel = "Re(s)  [σ]",
    ylabel = "Im(s)  [jω]")

p2 = plot(step(G, 15.0),
    title  = "UNIT STEP RESPONSE  |  K=$(round(K_dc,digits=4))",
    xlabel = "Time [s]",
    ylabel = "y(t)",
    label  = "y(t)",
    lw     = 2)

hline!(p2, [K_dc],
    linestyle = :dash,
    lw        = 1,
    label     = "DC Gain K=$(round(K_dc,digits=3))",
    color     = :orange)

fig = plot(p1, p2, layout = (1, 2), size = (1100, 500))
display(fig)

savefig(fig, "51a.png")
