using Printf

const BR = 150.0
const MX = 100.0
const SY = 50.0e3
const π  = pi

const Y = [
    750.0, 603.6, 250.0, -250.0, -603.6,
   -750.0, -603.6, -250.0, 250.0, 603.6
]

Ixx = 2 * BR * (
    750.0^2 +
    2 * 603.6^2 +
    2 * 250.0^2
)

σz = Float64[]

for (_, y_val) in enumerate(Y)
    push!(σz, (MX / Ixx) * y_val * 1e6)
end

C = -(SY * BR) / Ixx

qb = Float64[]

qb910 = C * Y[9]
qb101 = qb910 + C * Y[10]
qb12  = qb101 + C * Y[1]

push!(qb, qb12)

for i in 2:3
    push!(qb, qb[end] + C * Y[i])
end

# symmetry assignments
push!(qb, qb[2])
push!(qb, qb[1])
push!(qb, qb101)
push!(qb, qb910)
push!(qb, 0.0)
push!(qb, qb910)
push!(qb, qb101)

A34  = 0.5 * 500.0 * 250.0
A23  = 62500.0 +
       (45/360) * π * 500.0^2 -
       0.5 * 250.0 * 353.6

A910 = A23

A12  = 0.5 * 250.0 * 353.6 +
       (45/360) * π * 500.0^2

A101 = A12

A    = 500.0 * 1000.0 + π * 500.0^2

qso = (
    SY * 250.0
    - 2 * (
        2 * qb910 * A910 +
        2 * qb101 * A101 +
        2 * qb12  * A12  +
        2 * qb[2] * A23  +
        2 * qb[3] * A34
    )
) / (-2 * A)

q = Float64[]

for (_, val) in enumerate(qb)
    push!(q, abs(val - qso))
end

names = [
    "q12","q23","q34","q45","q56",
    "q67","q78","q89","q910","q101"
]

println("\n--- Shear Flow Results (Absolute Values) ---\n")

for (i, value) in enumerate(q)
    @printf("%-6s = %12.6f N/mm^2\n", names[i], value)
end
