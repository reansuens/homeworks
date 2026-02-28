using ControlSystems
using Plots
using Plots.Measures

s = tf("s")
G = 1 / (2*s^2 + 3*s + 1)   #Transfer Function defined here. and I can just change this part

sys_poles = poles(G)
sys_zeros = tzeros(G)
K_dc      = dcgain(G)[1]

println("  Poles       : ", sys_poles)
println("  Zeros       : ", isempty(sys_zeros) ? "[ none ]" : sys_zeros)

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
    top_margin    = 6mm)

display(fig)
savefig(fig, "52b.png") #switch the name of the file every time you run for a new hw
