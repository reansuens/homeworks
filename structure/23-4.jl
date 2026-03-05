using Printf
using LinearAlgebra

println("=" ^ 64)
println("  INIT: P.23.3 FOUR-CELL WING TORSIONAL STIFFNESS SOLVER")
println("=" ^ 64)

const A_I   = 161_500.0
const A_II  = 291_000.0
const A_III = 291_000.0
const A_IV  = 226_000.0

println("\n[SECTION 1] CELL AREAS (pre-specified)")
@printf("  A_I   = %9.1f mm^2\n", A_I)
@printf("  A_II  = %9.1f mm^2\n", A_II)
@printf("  A_III = %9.1f mm^2\n", A_III)
@printf("  A_IV  = %9.1f mm^2\n", A_IV)

const delta_12  = 762.0  / 0.915      # = 832.8
const delta_78  = delta_12
const delta_23  = 812.0  / 0.915      # = 887.4
const delta_67  = delta_23
const delta_34  = delta_23
const delta_56  = delta_23
const delta_45o = 1525.0 / 0.711      # = 2144.9
const delta_45i = 356.0  / 1.220      # = 291.8
const delta_36  = 406.0  / 1.625      # = 249.8
const delta_72  = 356.0  / 1.220      # = 291.8
const delta_81  = 254.0  / 0.915      # = 277.6

println("\n[SECTION 2] WALL FLEXIBILITY COEFFICIENTS delta = s/t (computed)")
@printf("  delta_12 = delta_78                  = %7.1f\n", delta_12)
@printf("  delta_23 = delta_67 = delta_34 = delta_56 = %7.1f\n", delta_23)
@printf("  delta_45deg (outer skin)             = %7.1f\n", delta_45o)
@printf("  delta_45i   (inner web, Cell I<->II) = %7.1f\n", delta_45i)
@printf("  delta_36    (web, Cell II<->III)     = %7.1f\n", delta_36)
@printf("  delta_72    (web, Cell III<->IV)     = %7.1f\n", delta_72)
@printf("  delta_81                             = %7.1f\n", delta_81)

const Sigma_I   = delta_45o + delta_45i
const Sigma_II  = delta_34  + delta_45i + delta_56 + delta_36
const Sigma_III = delta_23  + delta_36  + delta_67  + delta_72
const Sigma_IV  = delta_72  + delta_78  + delta_81  + delta_12

println("\n[SECTION 3] CELL PERIMETER SUMS Sigma_i (computed)")
@printf("  Sigma_I   = delta_45deg + delta_45i                   = %7.1f\n", Sigma_I)
@printf("  Sigma_II  = delta_34 + delta_45i + delta_56 + delta_36 = %7.1f\n", Sigma_II)
@printf("  Sigma_III = delta_23 + delta_36  + delta_67 + delta_72  = %7.1f\n", Sigma_III)
@printf("  Sigma_IV  = delta_72 + delta_78  + delta_81 + delta_12  = %7.1f\n", Sigma_IV)

row1 = [Sigma_I,     -delta_45i,   0.0,       0.0,       -2.0*A_I  ]
row2 = [-delta_45i,   Sigma_II,   -delta_36,   0.0,       -2.0*A_II ]
row3 = [0.0,         -delta_36,    Sigma_III, -delta_72,  -2.0*A_III]
row4 = [0.0,          0.0,        -delta_72,   Sigma_IV,  -2.0*A_IV ]
row5 = [2.0*A_I,      2.0*A_II,   2.0*A_III,  2.0*A_IV,   0.0      ]

A_mat = [row1'; row2'; row3'; row4'; row5']
b_vec = [0.0, 0.0, 0.0, 0.0, 1.0]

println("\n[SECTION 4] 5x5 SYSTEM MATRIX [A_mat | b_vec]")
println("  Unknowns: [q_I, q_II, q_III, q_IV, mu]")
@printf("  %-18s %-12s %-12s %-12s %-12s %-12s | RHS\n",
        "Row","q_I","q_II","q_III","q_IV","mu")
println("  " * "-"^87)
labels = ["Cell I compat.","Cell II compat.","Cell III compat.","Cell IV compat.","Torque (T=1)"]
for i in 1:5
    @printf("  %-18s", labels[i])
    for j in 1:5
        @printf(" %11.2f", A_mat[i,j])
    end
    @printf("  | %6.1f\n", b_vec[i])
end

function cramers_rule(A::Matrix{Float64}, b::Vector{Float64})::Vector{Float64}
    n = size(A, 1)
    d = det(A)
    x = Vector{Float64}(undef, n)
    for j in 1:n
        Aj       = copy(A)
        Aj[:, j] = b
        x[j]     = det(Aj) / d
    end
    return x
end

println("\n[SECTION 5] SOLVER: CRAMER'S RULE")
det_A = det(A_mat)
@printf("  det(A_mat) = %.6e\n", det_A)

x_sol = cramers_rule(A_mat, b_vec)

q_I   = x_sol[1]
q_II  = x_sol[2]
q_III = x_sol[3]
q_IV  = x_sol[4]
mu    = x_sol[5]    # G * (dtheta/dz) / T

println("\n  SOLUTION VECTOR (T = 1 N*mm normalisation):")
@printf("  q_I   = %12.6e  [N/mm per N*mm]\n", q_I)
@printf("  q_II  = %12.6e  [N/mm per N*mm]\n", q_II)
@printf("  q_III = %12.6e  [N/mm per N*mm]\n", q_III)
@printf("  q_IV  = %12.6e  [N/mm per N*mm]\n", q_IV)
@printf("  mu    = %12.6e  [G*(dtheta/dz)/T]\n", mu)

GJ_over_G    = 1.0 / mu
GJ_over_G_M6 = GJ_over_G * 1.0e-6

println("\n[SECTION 6] TORSIONAL STIFFNESS OUTPUT")
println("=" ^ 64)
@printf("  GJ_eff = T/(dtheta/dz) = %.4f x 10^6 * G  N*mm^2/rad\n", GJ_over_G_M6)
println("=" ^ 64)
