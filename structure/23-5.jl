using Printf
using LinearAlgebra

println("=" ^ 64)
const S_y  = 44_500.0
const A_I  = 232_000.0    # Cell I enclosed area [mm²]
const A_II = 258_000.0    # Cell II enclosed area [mm²]

const B_1 = 1290.0;  const B_6 = B_1
const B_2 = 1936.0;  const B_5 = B_2
const B_3 =  645.0;  const B_4 = B_3

const y_6 =  127.0;  const y_1 = -y_6
const y_5 =  203.0;  const y_2 = -y_5   # NOTE: boom 5 at web junction
const y_4 =  101.0;  const y_3 = -y_4

@printf("  S_y  = %.1f N\n", S_y)
@printf("  A_I  = %.1f mm^2\n", A_I)
@printf("  A_II = %.1f mm^2\n", A_II)
println("\n  Boom | Area [mm^2] | y [mm]")
println("  " * "-"^30)
for (lbl,B,y) in [("1",B_1,y_1),("2",B_2,y_2),("3",B_3,y_3),
                   ("4",B_4,y_4),("5",B_5,y_5),("6",B_6,y_6)]
    @printf("   %-4s  %8.1f     %8.1f\n", lbl, B, y)
end

const I_xx = B_1*y_1^2 + B_2*y_2^2 + B_3*y_3^2 +
             B_4*y_4^2 + B_5*y_5^2 + B_6*y_6^2

const C_qb = -S_y / I_xx    # shear flow coefficient

@printf("  I_xx = Σ(B_r·y_r²) = %.4e mm^4\n", I_xx)
@printf("  C_qb = -S_y/I_xx   = %.6e\n", C_qb)
@printf("  Textbook reference : 214.3 × 10^6 mm^4\n")

const d_16 = 254.0 / 1.625    # = 156.3
const d_25 = 406.0 / 2.032    # = 199.8  (shared web)
const d_34 = 202.0 / 1.220    # = 165.6
const d_12 = 647.0 / 0.915    # = 707.1
const d_56 = d_12
const d_23 = 775.0 / 0.559    # = 1386.4
const d_45 = d_23

@printf("  d_16         = %7.1f\n", d_16)
@printf("  d_25 (web)   = %7.1f\n", d_25)
@printf("  d_34         = %7.1f\n", d_34)
@printf("  d_12 = d_56  = %7.1f\n", d_12)
@printf("  d_23 = d_45  = %7.1f\n", d_23)

const qb_65 = 0.0
const qb_61 = qb_65 + C_qb * B_6 * y_6       # = -32.8 N/mm
const qb_12 = qb_61 + C_qb * B_1 * y_1       # ≈ 0 (symmetry confirmed)
const qb_25 = qb_12 + C_qb * B_2 * y_2       # = +81.7 N/mm (web 25)

const qb_54 = 0.0
const qb_43 = qb_54 + C_qb * B_4 * y_4       # raw, used for q_b,34
const qb_34 = -qb_43                           # = +13.6 N/mm (direction reversal)
const qb_23 = 0.0                              # by symmetry
const qb_45 = 0.0                              # by symmetry

println("  Wall  | q_b [N/mm]  | Notes")
println("  " * "-"^45)
@printf("  65    | %10.2f  | Cut wall (Cell I)\n",  qb_65)
@printf("  61    | %10.2f  | ref: -32.8\n",          qb_61)
@printf("  12    | %10.2f  | ~0 by symmetry\n",      qb_12)
@printf("  25    | %10.2f  | Shared web, ref: 81.7\n", qb_25)
@printf("  54    | %10.2f  | Cut wall (Cell II)\n",  qb_54)
@printf("  34    | %10.2f  | ref: 13.6\n",           qb_34)
@printf("  23    | %10.2f  | 0 by symmetry\n",       qb_23)
@printf("  45    | %10.2f  | 0 by symmetry\n",       qb_45)


const Sigma_I  = d_56 + d_16 + d_12 + d_25
const Sigma_II = d_45 + d_25 + d_23 + d_34

const qbs_I  =  qb_65*d_56 + qb_61*d_16 + qb_12*d_12 + qb_25*d_25
const qbs_II =  qb_45*d_45 - qb_25*d_25 + qb_23*d_23 + qb_34*d_34

const M_qb = qb_34 * 202.0 * 763.0 + qb_61 * 254.0 * 635.0

@printf("  Sigma_I   = %7.1f  ref: 1770.3\n", Sigma_I)
@printf("  Sigma_II  = %7.1f  ref: 3138.2\n", Sigma_II)
@printf("  qb·d (I)  = %9.2f  ref: +11197\n", qbs_I)
@printf("  qb·d (II) = %9.2f  ref: -14071.5\n", qbs_II)
@printf("  M_qb      = %12.2f  ref: -3194198\n", M_qb)

row1 = [ Sigma_I,   -d_25,    -2.0*A_I  ]
row2 = [-d_25,       Sigma_II, -2.0*A_II ]
row3 = [ 2.0*A_I,   2.0*A_II,  0.0      ]

A_mat = [row1'; row2'; row3']
b_vec = [-qbs_I, -qbs_II, -M_qb]

println("  Unknowns: [q_s0I, q_s0II, mu]")
@printf("  %-20s %12s %12s %12s | RHS\n", "Row","q_s0I","q_s0II","mu")
println("  " * "-"^72)
rlabels = ["Cell I compat.", "Cell II compat.", "Moment equil."]
for i in 1:3
    @printf("  %-20s", rlabels[i])
    for j in 1:3; @printf(" %11.3f", A_mat[i,j]); end
    @printf("  | %12.3f\n", b_vec[i])
end

function cramers_rule(A::Matrix{Float64}, b::Vector{Float64})::Vector{Float64}
    n = size(A, 1)
    d = det(A)
    x = Vector{Float64}(undef, n)
    for j in 1:n
        Aj = copy(A); Aj[:,j] = b
        x[j] = det(Aj) / d
    end
    return x
end

@printf("  det(A_mat) = %.6e\n", det(A_mat))

x_sol  = cramers_rule(A_mat, b_vec)
qs0_I  = x_sol[1]
qs0_II = x_sol[2]
mu     = x_sol[3]

@printf("  q_s0,I  = %8.4f N/mm  (ref: -1.1)\n", qs0_I)
@printf("  q_s0,II = %8.4f N/mm  (ref: +7.2)\n", qs0_II)
@printf("  mu      = %.6e\n", mu)

const q_16  = qb_61 + qs0_I          # wall 16 (= wall 61 reversed)
const q_21  = qb_65 + qs0_I          # wall 65 = wall 21 (= q_s0_I since qb_65=0)
const q_25f = qb_25 + qs0_I - qs0_II # web 25
const q_34f = qb_34 + qs0_II         # wall 34
const q_23f = qb_23 + qs0_II         # wall 23 (= q_s0_II since qb_23=0)
const q_45f = qb_45 + qs0_II         # wall 45 (= q_s0_II)

println("\n FINAL SHEAR FLOW DISTRIBUTION")
println("=" ^ 64)
println("  Wall     | q_final [N/mm] | Textbook [N/mm]")
println("  " * "-"^48)
@printf("  16       | %12.4f   | 33.9\n",  abs(q_16))
@printf("  65 = 21  | %12.4f   |  1.1\n",  abs(q_21))
@printf("  25 (web) | %12.4f   | 73.4\n",  abs(q_25f))
@printf("  34       | %12.4f   | 20.8\n",  abs(q_34f))
@printf("  23 = 45  | %12.4f   |  7.2\n",  abs(q_23f))
println("=" ^ 64)

