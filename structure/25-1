using LinearAlgebra
struct Material
    E::Float64      
    ν::Float64     
end

polyester = Material(3000.0, 0.16)
kevlar    = Material(140000.0, 0.28)

struct Layer
    material::Material
    thickness::Float64  # mm
end

layers = [
    Layer(polyester, 15.0),
    Layer(kevlar, 5.0),
    Layer(polyester, 15.0),
    Layer(kevlar, 5.0),
    Layer(polyester, 15.0)
]

b = 100.0  
L = 1000.0 
P = 500e3 

function compute_areas(layers::Vector{Layer}, width::Float64)
    areas = Float64[]
    for layer in layers
        push!(areas, layer.thickness * width)
    end
    return areas
end

A_vec = compute_areas(layers, b)
A_total = sum(A_vec)

println("AREA")
println("#################")
for (i, (layer, A)) in enumerate(zip(layers, A_vec))
    mat_name = layer.material === polyester ? "POLYESTER" : "KEVLAR"
    println("Layer $i [$mat_name]: A = $(A) mm^2")
end
println("TOTAL AREA: A_total = $(A_total) mm^2")
println()

E_eq = sum(layers[i].material.E * A_vec[i] for i in 1:length(layers)) / A_total

println("E_equivalent = $(E_eq) N/mm^2")
println("#################")
println()

δ = (P * L) / (E_eq * A_total)


println("DISPLACEMENT")
println("#################")
println("Axial shortening: δ = $(round(δ, digits=3)) mm")
println()

ε = δ / L

println("STRAIN")
println("#################")
println("Axial strain: ε = $(round(ε, sigdigits=6))")
println()

σ_polyester = polyester.E * ε
σ_kevlar = kevlar.E * ε

println("STRESS")
println("Polyester stress: σ_p = $(round(σ_polyester, digits=2)) N/mm^2")
println("Kevlar stress:    σ_k = $(round(σ_kevlar, digits=2)) N/mm^2")
println()

ε_trans_polyester = polyester.ν * ε
ε_trans_kevlar = kevlar.ν * ε

Δt_total = sum(
    layers[i].thickness * (layers[i].material === polyester ? ε_trans_polyester : ε_trans_kevlar)
    for i in 1:length(layers)
)


println("#################")
println("Thickness increase: Δt = $((round(Δt_total, digits=4))) mm")
println()

println("FORCE EQUILIBRIUM")

F_polyester = σ_polyester * sum(A_vec[i] for i in 1:length(layers) if layers[i].material === polyester)
F_kevlar = σ_kevlar * sum(A_vec[i] for i in 1:length(layers) if layers[i].material === kevlar)
F_total = F_polyester + F_kevlar

println("Force in polyester layers: F_p = $(round(F_polyester, digits=1)) N")
println("Force in kevlar layers:    F_k = $(round(F_kevlar, digits=1)) N")
println("Total internal force:      F_t = $(round(F_total, digits=1)) N")
println("Applied load:              P   = $(P) N")
println("Residual error:            ε_f = $(round(abs(F_total - P), sigdigits=3)) N")

println("│ Axial shortening (δ)         │ $(rpad(round(δ, digits=3), 18)) │ mm     │")

println("│ Thickness increase (Δt)      │ $(rpad(round(Δt_total, digits=4), 18)) │ mm     │")
println("│ Polyester stress (σ_p)       │ $(rpad(round(σ_polyester, digits=2), 18)) │ N/mm^2  │")
println("│ Kevlar stress (σ_k)          │ $(rpad(round(σ_kevlar, digits=2), 18)) │ N/mm^2  │")
