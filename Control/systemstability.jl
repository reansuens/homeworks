using ControlSystems
using Plots
using Plots.Measures

s = tf("s")
G = 1 / (2*s^2 + 3*s + 1)   #Transfer Function defined here. and I can just change this part

sys_poles = poles(G)
sys_zeros = tzeros(G)
K_dc      = dcgain(G)[1]

den_coeffs = vec(denpoly(G)[1].coeffs)  
c_d = den_coeffs[1]   
b_d = den_coeffs[2]  
a_d = den_coeffs[3] 

omega_n = sqrt(c_d / a_d)
zeta    = b_d / (2 * a_d * omega_n)

println("  Poles       : ", sys_poles)
println("  Zeros       : ", isempty(sys_zeros) ? "[ none ]" : sys_zeros)
println("  ωₙ  [rad/s] : ", round(omega_n, digits=6))
println("  ζ           : ", round(zeta,    digits=6))
println("  K  (DC gain): ", round(K_dc,    digits=6))

cls = zeta < 0  ? "UNSTABLE"          :
      zeta == 0 ? "UNDAMPED"          :
      zeta < 1  ? "UNDERDAMPED"       :
      zeta == 1 ? "CRITICALLY DAMPED" :
                  "OVERDAMPED"

ωd_val = zeta < 1 ? round(omega_n * sqrt(1 - zeta^2), digits=6) : NaN
println("  System class: $cls  |  ωd = $ωd_val rad/s")

p1 = pzmap(G,
    title  = "POLE-ZERO MAP",
    xlabel = "Re(s)",
    ylabel = "Im(s)",
    margin = 6mm)

p2 = plot(step(G, 15.0),
    title  = "UNIT STEP RESPONSE",
    xlabel = "Time [s]",
    ylabel = "y(t)",
    label  = "y(t)",
    lw     = 2,
    margin = 6mm)

hline!(p2, [K_dc],
    linestyle = :dash,
    lw        = 1,
    label     = "K = $(round(K_dc, digits=3))",
    color     = :orange)

fig = plot(p1, p2,
    layout      = (1, 2),
    size        = (1200, 520),
    bottom_margin = 8mm,
    left_margin   = 8mm,
    top_margin    = 6mm,
    plot_title  = "ωₙ = $(round(omega_n,digits=4)) rad/s  |  ζ = $(round(zeta,digits=4))  |  K = $(round(K_dc,digits=4))  |  $(cls)")

display(fig)
savefig(fig, "51a.png") #switch the name of the file every time you run for a new hw
