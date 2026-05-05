using CairoMakie
using Roots
using Printf

const X_mm  = Float64[0,25,50,75,100,125,150,175,200,225,250,275,
                      300,325,350,375,400,425,450,475,500,525,550,575,600]

const Yu_mm = Float64[10.7,24.4,37.5,47.1,58.2,66.2,67.7,63.4,59.0,55.2,
                      50.2,45.0,40.1,37.0,35.0,33.5,31.3,30.3,29.3,28.4,
                      28.5,28.0,27.0,25.4,25.4]

const YL_mm = Float64[10.7,20.7,25.0,24.5,25.5,25.5,25.5,25.5,25.5,25.5,
                      25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,25.5,
                      25.5,25.5,25.5,25.5,25.5]

const X_ho_raw  = Float64[0,50,100,150,200,250,300,350,400,450,500,550,600]
const ho_cm_raw = Float64[10.2,16.5,23.3,35.4,43.2,58.2,61.5,62.2,63.8,63.4,56.4,61.6,59.6]

const S_star = 59.2
const k      = 1.4
const P_star_Po = (2/(k+1))^(k/(k-1))   # ≈ 0.5283

S      = 152.0 .- (Yu_mm .+ YL_mm)
AAstar = S ./ S_star

function interp1(xq::Float64, xs::Vector{Float64}, ys::Vector{Float64})
    idx = searchsortedfirst(xs, xq)
    idx == 1          && return ys[1]
    idx > length(xs)  && return ys[end]
    x0, x1 = xs[idx-1], xs[idx]
    y0, y1 = ys[idx-1], ys[idx]
    return y0 + (y1-y0)*(xq-x0)/(x1-x0)
end

ho_interp = [interp1(x, X_ho_raw, ho_cm_raw) for x in X_mm]
P_Po      = 1.0 .- (ho_interp ./ 100.0) ./ 0.76

function M_from_PPo(PPo::Float64)
    (PPo <= 0.0 || PPo >= 1.0) && return NaN
    arg = PPo^(-(k-1)/k) - 1.0
    arg < 0.0 && return NaN
    return sqrt(2.0/(k-1) * arg)
end

M_exp = M_from_PPo.(P_Po)

function AAstar_func(M::Float64)
    return (1/M) * ((2/(k+1)) * (1 + (k-1)/2 * M^2))^((k+1)/(2*(k-1)))
end

const throat_idx = 7   # X=150 mm, 1-indexed

function solve_M_th(ratio::Float64, supersonic::Bool)
    ratio < 1.0 && return NaN
    f = M -> AAstar_func(M) - ratio
    if supersonic
        try; return find_zero(f, (1.0+1e-9, 10.0), Bisection()); catch; return NaN; end
    else
        try; return find_zero(f, (1e-9, 1.0-1e-9), Bisection()); catch; return NaN; end
    end
end

M_th = [i == throat_idx ? 1.0 : solve_M_th(AAstar[i], i > throat_idx)
        for i in eachindex(X_mm)]

DM = [isnan(M_exp[i]) || isnan(M_th[i]) ? NaN : M_exp[i] - M_th[i]
      for i in eachindex(X_mm)]

set_theme!(theme_dark())

const CYAN    = RGBf(0.0,  0.90, 0.95)
const ORANGE  = RGBf(1.0,  0.60, 0.10)
const GREEN   = RGBf(0.20, 0.95, 0.45)
const MAGENTA = RGBf(0.95, 0.20, 0.75)
const RED     = RGBf(0.95, 0.25, 0.25)
const YELLOW  = RGBf(0.98, 0.92, 0.20)
const BGCOL   = RGBf(0.07, 0.07, 0.10)
const GRIDCOL = RGBAf(0.35, 0.35, 0.45, 0.5)

AXIS_KW = (
    backgroundcolor = BGCOL,
    xgridcolor      = GRIDCOL,
    ygridcolor      = GRIDCOL,
    xgridwidth      = 0.6,
    ygridwidth      = 0.6,
    titlesize       = 13,
    xlabelsize      = 11,
    ylabelsize      = 11,
    xticklabelsize  = 9,
    yticklabelsize  = 9,
    titlecolor      = :white,
    xlabelcolor     = :white,
    ylabelcolor     = :white,
)

valid_exp = .!isnan.(M_exp)
valid_th  = .!isnan.(M_th)
valid_dM  = .!isnan.(DM)
X_dM      = X_mm[valid_dM]
dM_v      = DM[valid_dM]


fig1 = Figure(size=(900,420), backgroundcolor=BGCOL)

ax1a = Axis(fig1[1,1];
    title="Fig 1 — Nozzle Wall Profile  Yᵤ, Y_L  vs  X",
    xlabel="X  (mm)", ylabel="Wall offset  (mm)", AXIS_KW...)
lines!(ax1a, X_mm, Yu_mm; color=CYAN,   linewidth=2.0, label="Yᵤ  (upper wall)")
lines!(ax1a, X_mm, YL_mm; color=ORANGE, linewidth=2.0, label="Y_L  (lower wall)")
scatter!(ax1a, X_mm, Yu_mm; color=CYAN,   markersize=5)
scatter!(ax1a, X_mm, YL_mm; color=ORANGE, markersize=5)
vlines!(ax1a, [150.0]; color=RED, linestyle=:dash, linewidth=1.2, label="Throat  X=150 mm")
axislegend(ax1a; position=:lt, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

ax1b = Axis(fig1[1,2];
    title="Fig 2 — Gap Height  S  vs  X",
    xlabel="X  (mm)", ylabel="S  (mm)", AXIS_KW...)
lines!(ax1b, X_mm, S; color=GREEN, linewidth=2.2, label="S = 152 − (Yᵤ+Y_L)")
scatter!(ax1b, X_mm, S; color=GREEN, markersize=5)
hlines!(ax1b, [S_star]; color=RED,     linestyle=:dash, linewidth=1.2, label="S* = 59.2 mm")
vlines!(ax1b, [150.0];  color=MAGENTA, linestyle=:dash, linewidth=1.0, label="Throat")
axislegend(ax1b; position=:rt, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig1_NozzleProfile.png", fig1, px_per_unit=2)
println("Saved: Fig1_NozzleProfile.png")

fig2 = Figure(size=(820,400), backgroundcolor=BGCOL)
ax2 = Axis(fig2[1,1];
    title="Fig 3 — Area Ratio  A/A*  vs  X",
    xlabel="X  (mm)", ylabel="A / A*", AXIS_KW...)
lines!(ax2, X_mm, AAstar; color=CYAN, linewidth=2.2, label="A/A* = S/S*")
scatter!(ax2, X_mm, AAstar; color=CYAN, markersize=5)
hlines!(ax2, [1.0]; color=RED,     linestyle=:dash, linewidth=1.2, label="Throat  A/A* = 1")
vlines!(ax2, [150.0]; color=MAGENTA, linestyle=:dash, linewidth=1.0, label="X = 150 mm")
axislegend(ax2; position=:rt, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig2_AreaRatio.png", fig2, px_per_unit=2)
println("Saved: Fig2_AreaRatio.png")

fig3 = Figure(size=(820,440), backgroundcolor=BGCOL)
ax3 = Axis(fig3[1,1];
    title="Fig 4 — Static-to-Total Pressure Ratio  P/P₀  vs  X",
    xlabel="X  (mm)", ylabel="P / P₀", AXIS_KW...)

lines!(ax3, X_mm, P_Po; color=ORANGE, linewidth=2.4, label="P/P₀  (experimental)")
scatter!(ax3, X_mm, P_Po; color=ORANGE, markersize=6)
hlines!(ax3, [P_star_Po]; color=RED, linestyle=:dash, linewidth=1.3,
        label="P*/P₀ = $(round(P_star_Po, digits=4))")
hlines!(ax3, [1.0]; color=GRIDCOL, linestyle=:dot, linewidth=1.0)
vlines!(ax3, [150.0]; color=MAGENTA, linestyle=:dash, linewidth=1.0, label="Throat  X=150 mm")
band!(ax3, [450.0, 525.0], [0.0, 0.0], [1.0, 1.0]; color=RGBAf(0.95,0.25,0.25,0.10))
text!(ax3, 487, 0.62; text="Shock /\nsep. zone", color=RED, fontsize=9, align=(:center,:center))
axislegend(ax3; position=:lb, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig3_PressureRatio.png", fig3, px_per_unit=2)
println("Saved: Fig3_PressureRatio.png")

fig4 = Figure(size=(900,460), backgroundcolor=BGCOL)
ax4 = Axis(fig4[1,1];
    title="Fig 5 — Mach Number Distribution  M_exp  vs  M_th  vs  X",
    xlabel="X  (mm)", ylabel="Mach number  M", AXIS_KW...)

lines!(ax4, X_mm[valid_exp], M_exp[valid_exp];
       color=CYAN,   linewidth=2.4, label="M_exp  (isentropic from P/P₀)")
scatter!(ax4, X_mm[valid_exp], M_exp[valid_exp]; color=CYAN,   markersize=5)
lines!(ax4, X_mm[valid_th],  M_th[valid_th];
       color=ORANGE, linewidth=2.0, linestyle=:dash, label="M_th  (A/A* isentropic)")
scatter!(ax4, X_mm[valid_th],  M_th[valid_th];  color=ORANGE, markersize=4)
hlines!(ax4, [1.0];   color=RED,     linestyle=:dash, linewidth=1.2, label="M = 1  (sonic)")
vlines!(ax4, [150.0]; color=MAGENTA, linestyle=:dash, linewidth=1.0, label="Throat")
band!(ax4, [450.0, 525.0], [0.0, 0.0], [2.5, 2.5]; color=RGBAf(0.95,0.25,0.25,0.10))
text!(ax4, 487, 0.80; text="Shock /\nsep. zone", color=RED, fontsize=9, align=(:center,:center))
axislegend(ax4; position=:lt, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig4_MachDistribution.png", fig4, px_per_unit=2)
println("Saved: Fig4_MachDistribution.png")

fig5 = Figure(size=(900,440), backgroundcolor=BGCOL)
ax5 = Axis(fig5[1,1];
    title="Fig 6 — Mach Deviation  ΔM = M_exp − M_th  vs  X",
    xlabel="X  (mm)", ylabel="ΔM  =  M_exp − M_th", AXIS_KW...)

hlines!(ax5, [0.0]; color=:white, linewidth=0.8, linestyle=:dot)
band!(ax5, X_dM, min.(dM_v, 0.0), max.(dM_v, 0.0);
      color=[v >= 0 ? RGBAf(0.20,0.95,0.45,0.25) : RGBAf(0.95,0.25,0.25,0.20)
             for v in dM_v])
lines!(ax5, X_dM, dM_v;   color=GREEN, linewidth=2.4, label="ΔM = M_exp − M_th")
scatter!(ax5, X_dM, dM_v; color=GREEN, markersize=6)
vlines!(ax5, [150.0]; color=RED, linestyle=:dash, linewidth=1.2, label="Throat  X=150 mm")
band!(ax5, [450.0, 525.0], [-0.55, -0.55], [0.0, 0.0]; color=RGBAf(0.95,0.92,0.20,0.10))
text!(ax5, 75,   0.14; text="ΔM > 0\nBL displacement\n(convergent)",
      color=GREEN,  fontsize=9, align=(:center,:center))
text!(ax5, 340, -0.18; text="ΔM < 0\nViscous loss\n(divergent)",
      color=RED,    fontsize=9, align=(:center,:center))
text!(ax5, 490, -0.51; text="Shock /\nseparation",
      color=YELLOW, fontsize=9, align=(:center,:center))
axislegend(ax5; position=:rb, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig5_MachDeviation.png", fig5, px_per_unit=2)
println("Saved: Fig5_MachDeviation.png")

fig6 = Figure(size=(900,460), backgroundcolor=BGCOL)
ax6 = Axis(fig6[1,1];
    title="Fig 7 — Back-Pressure Regimes and Choked Mass Flow  (Theoretical)",
    xlabel="Axial position  x/L  (normalised)",
    ylabel="P / P₀", AXIS_KW...)

x_norm   = collect(range(0.0, 1.0, length=400))
throat_x = 0.25

function lineB(x)
    x <= throat_x ?
        1.0 - 0.18*(x/throat_x) :
        0.82 + 0.10*((x-throat_x)/(1-throat_x))
end
function lineC(x)
    if x <= throat_x
        return 1.0 - 0.47*(x/throat_x)
    elseif x <= 0.65
        t = (x - throat_x)/(0.65 - throat_x)
        return 0.53 - 0.25*t
    else
        t = (x - 0.65)/(1 - 0.65)
        return 0.28 + 0.42*t
    end
end
function lineD(x)
    x <= throat_x ?
        1.0 - 0.47*(x/throat_x) :
        0.53 - 0.35*((x - throat_x)/(1 - throat_x))^0.65
end

lB = lineB.(x_norm)
lC = lineC.(x_norm)
lD = lineD.(x_norm)

lines!(ax6, x_norm, lB; color=CYAN,    linewidth=2.0, label="Line B — subsonic, unchoked")
lines!(ax6, x_norm, lC; color=YELLOW,  linewidth=2.0, label="Line C — choked, normal shock in divergent")
lines!(ax6, x_norm, lD; color=MAGENTA, linewidth=2.0, label="Line D — fully supersonic (experimental condition)")
hlines!(ax6, [P_star_Po]; color=RED, linestyle=:dash, linewidth=1.2,
        label="P*/P₀ = $(round(P_star_Po, digits=4))")
vlines!(ax6, [throat_x]; color=:white, linestyle=:dot, linewidth=1.0, label="Throat")
text!(ax6, 0.5, 0.95;
      text="ṁ_max  (choked, Lines C & D  — mass flow independent of Pb)",
      color=:white, fontsize=9, align=(:center,:center))
axislegend(ax6; position=:lb, backgroundcolor=RGBAf(0,0,0,0.6), labelsize=10)

save("Fig6_BackPressureRegimes.png", fig6, px_per_unit=2)
println("Saved: Fig6_BackPressureRegimes.png")


fig7 = Figure(size=(1400,900), backgroundcolor=BGCOL)

Label(fig7[0, 1:3],
      "Experiment 4 — Convergent–Divergent Nozzle: Pressure Distribution Analysis";
      fontsize=15, color=:white, font=:bold)

axA = Axis(fig7[1,1]; title="(A)  S  vs  X", xlabel="X (mm)", ylabel="S (mm)", AXIS_KW...)
lines!(axA, X_mm, S; color=GREEN, linewidth=2.0)
hlines!(axA, [S_star]; color=RED,     linestyle=:dash, linewidth=1.0)
vlines!(axA, [150.0];  color=MAGENTA, linestyle=:dash, linewidth=0.8)

axB = Axis(fig7[1,2]; title="(B)  P/P₀  vs  X", xlabel="X (mm)", ylabel="P/P₀", AXIS_KW...)
lines!(axB, X_mm, P_Po; color=ORANGE, linewidth=2.0)
scatter!(axB, X_mm, P_Po; color=ORANGE, markersize=4)
hlines!(axB, [P_star_Po]; color=RED,     linestyle=:dash, linewidth=1.0)
vlines!(axB, [150.0];     color=MAGENTA, linestyle=:dash, linewidth=0.8)
band!(axB, [450.0, 525.0], [0.0, 0.0], [1.0, 1.0]; color=RGBAf(0.95,0.25,0.25,0.10))

axC = Axis(fig7[1,3]; title="(C)  M_exp & M_th  vs  X", xlabel="X (mm)", ylabel="M", AXIS_KW...)
lines!(axC, X_mm[valid_exp], M_exp[valid_exp]; color=CYAN,   linewidth=2.0, label="M_exp")
lines!(axC, X_mm[valid_th],  M_th[valid_th];  color=ORANGE, linewidth=2.0, linestyle=:dash, label="M_th")
hlines!(axC, [1.0]; color=RED,     linestyle=:dash, linewidth=1.0)
vlines!(axC, [150.0]; color=MAGENTA, linestyle=:dash, linewidth=0.8)
axislegend(axC; position=:lt, labelsize=8, backgroundcolor=RGBAf(0,0,0,0.5))

axD = Axis(fig7[2,1]; title="(D)  ΔM  vs  X", xlabel="X (mm)", ylabel="ΔM", AXIS_KW...)
hlines!(axD, [0.0]; color=:white, linewidth=0.7, linestyle=:dot)
band!(axD, X_dM, min.(dM_v,0.0), max.(dM_v,0.0);
      color=[v>=0 ? RGBAf(0.20,0.95,0.45,0.25) : RGBAf(0.95,0.25,0.25,0.20)
             for v in dM_v])
lines!(axD, X_dM, dM_v; color=GREEN, linewidth=2.0)
scatter!(axD, X_dM, dM_v; color=GREEN, markersize=4)
vlines!(axD, [150.0]; color=RED, linestyle=:dash, linewidth=1.0)
band!(axD, [450.0, 525.0], [-0.6, -0.6], [0.0, 0.0]; color=RGBAf(0.95,0.92,0.20,0.10))

# E: Area ratio
axE = Axis(fig7[2,2]; title="(E)  A/A*  vs  X", xlabel="X (mm)", ylabel="A/A*", AXIS_KW...)
lines!(axE, X_mm, AAstar; color=CYAN, linewidth=2.0)
hlines!(axE, [1.0]; color=RED,     linestyle=:dash, linewidth=1.0)
vlines!(axE, [150.0]; color=MAGENTA, linestyle=:dash, linewidth=0.8)

# F: Back-pressure regimes
axF = Axis(fig7[2,3]; title="(F)  Back-Pressure Regimes", xlabel="x/L", ylabel="P/P₀", AXIS_KW...)
lines!(axF, x_norm, lB; color=CYAN,    linewidth=1.8, label="Line B")
lines!(axF, x_norm, lC; color=YELLOW,  linewidth=1.8, label="Line C")
lines!(axF, x_norm, lD; color=MAGENTA, linewidth=1.8, label="Line D")
hlines!(axF, [P_star_Po]; color=RED, linestyle=:dash, linewidth=1.0)
axislegend(axF; position=:lb, labelsize=8, backgroundcolor=RGBAf(0,0,0,0.5))

save("Fig7_CombinedPanel.png", fig7, px_per_unit=2)
println("Saved: Fig7_CombinedPanel.png")


println("\n" * "="^72)
println("  EXPERIMENT 4 — C-D NOZZLE  |  COMPUTED RESULTS")
println("="^72)
@printf("  %-6s  %-7s  %-6s  %-7s  %-7s  %-7s  %-7s  %-8s\n",
        "X(mm)", "S(mm)", "A/A*", "ho(cm)", "P/P0", "M_exp", "M_th", "DM")
println("  " * "-"^70)
for i in eachindex(X_mm)
    @printf("  %6.0f  %7.1f  %6.3f  %7.2f  %7.4f  %7.4f  %7.4f  %+8.4f\n",
            X_mm[i], S[i], AAstar[i], ho_interp[i],
            P_Po[i], M_exp[i], M_th[i], DM[i])
end
println("="^72)
println()
println("  Throat:  X = 150 mm  |  S* = $(S_star) mm  |  P*/P₀ = $(round(P_star_Po,digits=4))")
println()
println("  Output files:")
for f in ["Fig1_NozzleProfile.png","Fig2_AreaRatio.png","Fig3_PressureRatio.png",
          "Fig4_MachDistribution.png","Fig5_MachDeviation.png",
          "Fig6_BackPressureRegimes.png","Fig7_CombinedPanel.png"]
    println("    $f")
end
println()
