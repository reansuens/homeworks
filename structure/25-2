struct Material
    E::Float64
end

cover = Material(20000.0)
web = Material(60000.0)

b = 150.0       
h = 100.0      
t_cover = 1.0 
t_web = 2.0  
P = 40000.0 

EA_covers = cover.E * 2 * t_cover * b
EA_webs = web.E * 2 * t_web * h
EA_total = EA_covers + EA_webs

ε = P / EA_total

F_covers = cover.E * ε * 2 * t_cover * b
F_webs = web.E * ε * 2 * t_web * h

println("Stiffness for both")
println("EA_covers: $(EA_covers) N")
println("EA_webs: $(EA_webs) N")
println("EA_total: $(EA_total) N")
println()
println("STRAIN")
println("ε: $(round(ε, sigdigits=3))")
println()
println("AXIAL FORCES")
println("Covers (both): $(round(F_covers/1000, digits=1)) kN")
println("Webs (both): $(round(F_webs/1000, digits=1)) kN")
println("Total: $(round((F_covers + F_webs)/1000, digits=1)) kN")
