g          = 9.81
rho        = 1.225
W          = 238000.0
S          = 88.5
V          = 215.0
theta      = 40.0 * π / 180.0
R          = 1525.0
ddot_theta = -0.25          
l_t        = 12.2
I_yy       = 204000.0

q_dyn = 0.5 * rho * V^2

L_total = W * cos(theta) + (W / g) * (V^2 / R)
n_cg    = L_total / W

C_L = L_total / (q_dyn * S)

P = 0.0
for iter = 1:50
    global C_L, P
    M_cg   = q_dyn * S * (0.427 * C_L - 0.061)
    P_new  = (M_cg - I_yy * ddot_theta) / l_t
    L_wing = L_total - P_new
    C_L_new = L_wing / (q_dyn * S)
    err = abs(C_L_new - C_L)
    global C_L = C_L_new
    global P   = P_new
    err < 1e-10 && break
end

f = -P / W

n_tail_structural = n_cg + (abs(ddot_theta) * l_t) / g

M_cg_final = q_dyn * S * (0.427 * C_L - 0.061)

println("q_dyn  : $(round(q_dyn,              digits=1)) Pa")
println("L_total: $(round(L_total,            digits=0)) N")
println("n_CG   : $(round(n_cg,               digits=3))")
println("C_L    : $(round(C_L,                digits=4))")
println("M_CG   : $(round(M_cg_final,         digits=0)) N·m")
println("P      : $(round(P,                  digits=0)) N")
println("f      : $(round(f,                  digits=4))")
println("n_tail (structural): $(round(n_tail_structural, digits=3))")
