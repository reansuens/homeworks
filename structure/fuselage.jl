using Printf

Br  = 150.0                 # mm^2
Mx  = 100.0                 # kN·m
Sy  = 50.0e3                # N
π   = pi

y = [750.0, 603.6, 250.0, -250.0, -603.6, -750.0,
     -603.6, -250.0, 250.0, 603.6]

Ixx = 2 * Br * (750.0^2 + 2*603.6^2 + 2*250.0^2)

σz = (Mx / Ixx) .* y .* 1e6

C = -(Sy * Br) / Ixx

qb89  = 0.0
qb910 = C * y[9]
qb78  = qb910
qb101 = qb910 + C * y[10]
qb67  = qb101
qb12  = qb101 + C * y[1]
qb56  = qb12
qb23  = qb12 + C * y[2]
qb45  = qb23
qb34  = qb23 + C * y[3]

qb = [qb12, qb23, qb34, qb45, qb56,
      qb67, qb78, qb89, qb910, qb101]

A34  = 0.5 * 500.0 * 250.0
A23  = 62500.0 + (45/360) * π * 500.0^2 - 0.5 * 250.0 * 353.6
A910 = A23
A12  = 0.5 * 250.0 * 353.6 + (45/360) * π * 500.0^2
A101 = A12
A    = 500.0*1000.0 + π*500.0^2

qso = (
    Sy*250.0
    - 2*(2*qb910*A910 +
         2*qb101*A101 +
         2*qb12*A12 +
         2*qb23*A23 +
         2*qb34*A34)
) / (-2*A)

q = abs.(qb .- qso)

names = ["q12","q23","q34","q45","q56",
         "q67","q78","q89","q910","q101"]

println("\n--- Shear Flow Results (Absolute Values) ---\n")

for i in eachindex(q)
    @printf("%-6s = %12.6f N/mm^2\n", names[i], q[i])
end
