using ControlSystems
using Plots
using Plots.Measures

s = tf("s")
G = (s+3)/(s^2 + 3*s + 5)

p_rlocus = rlocusplot(G, 
    title  = "G(s)", 
    xlabel = "Real [σ]",
    ylabel = "Imaginary [jω]",
    lw     = 2,
    grid   = :both,
    size   = (800, 600),
    margin = 8mm)

display(p_rlocus)

p = poles(G)
println("Open-loop Poles (k=0): ", p)
