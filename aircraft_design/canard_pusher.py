import math
from dataclasses import dataclass, replace

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, Slider

DEG, G0 = math.pi / 180, 9.80665
INK, MUTED, RULE, GRID = "#0e1a25", "#566574", "#c7d0d8", "#dce3e9"
BG, PANEL = "#edf0f3", "#f8fafb"
CAN, WING, GOOD, WARN, BAD = "#1f5fcf", "#b36a00", "#1b7f47", "#9a6100", "#c2301d"

plt.rcParams.update({"mathtext.fontset": "cm", "text.color": INK, "axes.labelcolor": MUTED, "font.size": 9})


@dataclass
class P:
    m: float; Iy: float; Sw: float; Sc: float; cbar: float      # mass, inertia, areas, mean chord
    Xc: float; Xw: float; Xp: float; zT: float; Sp: float        # stations (m aft of nose), thrust height, disc area
    Aw: float; Ac: float; clw: float; clc: float; kdw: float     # lift slopes, CLmax, canard downwash on wing
    Pav: float; Tst: float; dcMax: float                         # power (W), static thrust (N), canard travel (rad)
    Xcg: float = 0.0
    cm0: float = -0.03; cd0: float = 0.025; kiw: float = 0.04; kic: float = 0.06
    cnp: float = 0.2; iw: float = 0.02; ic: float = 0.0; kup: float = 0.10


PRESETS = [
    dict(name="Starship", V0=125, h0=3000, Kn0=0.15,
         p=dict(m=6760, Iy=38000, Sw=26, Sc=5.2, cbar=1.6, Xc=2.9, Xw=9.4, Xp=12.2, zT=0.15, Sp=9.8,
                Aw=4.6, Ac=4.0, clw=1.45, clc=1.6, kdw=0.15, Pav=1.43e6, Tst=26e3, dcMax=0.26)),
    dict(name="Avtek 400A", V0=100, h0=3000, Kn0=0.12,
         p=dict(m=3900, Iy=11000, Sw=15.5, Sc=3.2, cbar=1.45, Xc=2.2, Xw=7.2, Xp=8.9, zT=1.0, Sp=5.85,
                Aw=4.9, Ac=4.2, clw=1.5, clc=1.6, kdw=0.08, Pav=0.85e6, Tst=14e3, dcMax=0.26)),
    dict(name="XP-55 Ascender", V0=120, h0=3000, Kn0=0.05,
         p=dict(m=3300, Iy=6600, Sw=21.8, Sc=2.8, cbar=1.8, Xc=1.8, Xw=6.3, Xp=8.7, zT=0.1, Sp=10.5,
                Aw=3.8, Ac=3.8, clw=0.85, clc=1.6, kdw=0.12, Pav=0.76e6, Tst=12e3, dcMax=0.17)),
]


def rho_at(h):
    return 1.225 * math.exp(-max(h, 0.0) / 9300)


def thrust_max(p, V):
    return min(p.Tst, p.Pav / max(V, 1.0))


def cl(a, slope, clmax):
    """Linear lift to stall, then a drop to 40 percent of CLmax over 0.15 rad."""
    a_s = clmax / slope
    x = abs(a)
    if x <= a_s:
        return slope * a
    return math.copysign(clmax * (1 - 0.6 * min(1.0, (x - a_s) / 0.15)), a)


def depth(a, slope, clmax):
    return min(1.0, max(0.0, (abs(a) - clmax / slope) / 0.15))


def aero(p, V, al, q, dc, T, rho):
    """Returns L, D, M, Lc, Lw, canard local alpha, wing local alpha."""
    qb = 0.5 * rho * V * V
    lc, lw, lp = p.Xcg - p.Xc, p.Xw - p.Xcg, p.Xp - p.Xcg
    ac = (1 + p.kup) * al + p.ic + dc - q * lc / V
    CLc = cl(ac, p.Ac, p.clc)
    aw = al + p.iw - p.kdw * CLc / p.Ac + q * lw / V
    CLw = cl(aw, p.Aw, p.clw)
    sc, sw = math.sin(ac), math.sin(aw)
    CDw = p.cd0 + p.kiw * CLw**2 + 1.1 * sw * sw * depth(aw, p.Aw, p.clw)
    CDc = p.kic * CLc**2 + 1.1 * sc * sc * depth(ac, p.Ac, p.clc)
    Lc, Lw = qb * p.Sc * CLc, qb * p.Sw * CLw
    D = qb * (p.Sw * CDw + p.Sc * CDc)
    M = Lc * lc - Lw * lw + qb * p.Sw * p.cbar * p.cm0 - T * p.zT - qb * p.Sp * p.cnp * al * lp
    return Lc + Lw, D, M, Lc, Lw, ac, aw


def deriv(p, s, dc, tau):
    V, al, q, th, h = max(s[0], 8.0), s[1], s[2], s[3], s[4]
    T = tau * thrust_max(p, V)
    L, D, M, *_ = aero(p, V, al, q, dc, T, rho_at(h))
    g = th - al
    return np.array([(T * math.cos(al) - D) / p.m - G0 * math.sin(g),
                     q - (L + T * math.sin(al)) / (p.m * V) + G0 * math.cos(g) / V,
                     M / p.Iy, q, V * math.sin(g)])


def rk4(p, s, dc, tau, dt):
    k1 = deriv(p, s, dc, tau)
    k2 = deriv(p, s + dt / 2 * k1, dc, tau)
    k3 = deriv(p, s + dt / 2 * k2, dc, tau)
    k4 = deriv(p, s + dt * k3, dc, tau)
    return s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def trim_at(p, V, h, guess):
    """Level flight by Newton. Unknowns: alpha, canard deflection, thrust. Returns (converged, z)."""
    rho, z = rho_at(h), np.array(guess, float)

    def res(z):
        L, D, M, *_ = aero(p, V, z[0], 0.0, z[1], z[2], rho)
        return np.array([z[2] * math.cos(z[0]) - D, L + z[2] * math.sin(z[0]) - p.m * G0, M])

    st = (1e-6, 1e-6, 1e-3)
    for _ in range(80):
        r = res(z)
        if np.linalg.norm(r) < 1e-3:
            break
        J = np.empty((3, 3))
        for j in range(3):
            d = np.zeros(3)
            d[j] = st[j]
            J[:, j] = (res(z + d) - r) / st[j]
        try:
            dz = np.linalg.solve(J, r)
        except np.linalg.LinAlgError:
            break
        sc = min(1.0, 0.1 / (abs(dz[0]) + 1e-12), 0.1 / (abs(dz[1]) + 1e-12), 8000 / (abs(dz[2]) + 1e-9))
        z = z - sc * dz
    ok = bool(np.all(np.isfinite(z)) and np.linalg.norm(res(z)) < 1.0)
    return ok, z


def trim_any(p, V, h):
    first = None
    for g in ([0.05, 0.05, 5000.0], [0.1, 0.0, 9000.0], [0.02, -0.05, 3000.0], [0.2, 0.1, 12000.0]):
        ok, z = trim_at(p, V, h, g)
        if ok:
            return True, z
        first = first or (ok, z)
    return first


def neutral_point(p):
    kc = p.Sc * p.Ac * (1 + p.kup)
    kw = p.Sw * p.Aw * (1 - p.kdw * (1 + p.kup))
    kp = p.Sp * p.cnp
    return (kc * p.Xc + kw * p.Xw + kp * p.Xp) / (kc + kw + kp), kc, kw, kp


def linearize(p, V, h, z):
    al, dc, T = z
    tau = T / thrust_max(p, V)
    s0 = np.array([V, al, 0.0, al, h])
    st = (1e-3, 1e-6, 1e-6, 1e-6)
    A = np.zeros((4, 4))
    for j in range(4):
        d = np.zeros(5)
        d[j] = st[j]
        A[:, j] = (deriv(p, s0 + d, dc, tau)[:4] - deriv(p, s0 - d, dc, tau)[:4]) / (2 * st[j])
    return A


def stall_sweep(p, V0, h, z):
    """Trim at falling speed and record how close each surface is to its stall angle."""
    rows, end, rho = [], "slow", rho_at(h)
    for V in np.arange(V0, 17.9, -1.0):
        ok, zn = trim_at(p, V, h, z)
        if not ok:
            end = "lost"
            break
        z = zn
        _, _, _, _, _, ac, aw = aero(p, V, z[0], 0.0, z[1], z[2], rho)
        fc, fw = abs(ac) * p.Ac / p.clc, abs(aw) * p.Aw / p.clw
        rows.append((V, fc, fw))
        if max(fc, fw) >= 1:
            end = "stall"
            break
        if abs(z[1]) > p.dcMax:
            end = "authority"
            break
    return np.array(rows), end


def simulate(p, V0, h0, z, kind, dt=0.025):
    al, dc0, T0 = z
    tau0 = min(1.0, T0 / thrust_max(p, V0))
    t_end = (20.0, 40.0, 45.0)[kind]
    s = np.array([V0, al, 0.0, al, h0])
    out, note = [(0.0, *s)], ""
    for i in range(int(t_end / dt)):
        t = i * dt
        dc, tau = dc0, tau0
        if kind == 0 and 1.0 <= t < 2.0:
            dc += DEG
        elif kind == 1 and t >= 1.0:
            tau = min(1.0, 1.3 * tau0)
        elif kind == 2:
            tau, dc = 0.0, dc0 + min(5 * DEG, 0.2 * DEG * t)
        s = rk4(p, s, float(np.clip(dc, -p.dcMax, p.dcMax)), tau, dt)
        if not np.all(np.isfinite(s)):
            note = "numerical failure"
            break
        if s[4] <= 0:
            note = "ground contact"
            break
        if abs(s[1]) > 1.2:
            note = "angle of attack out of range"
            break
        out.append((t + dt, *s))
    return np.array(out), note      # columns: t, V, alpha, q, theta, h


KINDS = ["1 deg canard pulse", "+30% power step", "idle power, slow pull"]
BOUNDS = {"Xcg": lambda b: (b.Xc + 0.5, b.Xw + 0.4), "Sc": lambda b: (0.4 * b.Sc, 2.2 * b.Sc),
          "zT": lambda b: (-1.0, 1.8), "Xp": lambda b: (0.5, b.Xp + 3.0)}


def make(pr):
    p = P(**pr["p"])
    return replace(p, Xcg=neutral_point(p)[0] - pr["Kn0"] * p.cbar)


fig = plt.figure(figsize=(13.2, 8.2), num="Canard flight bench", facecolor=BG)
fig.text(0.06, 0.985, "Canard flight bench", fontsize=15, weight="bold", va="top")
ax_bal = fig.add_axes([0.06, 0.855, 0.55, 0.085], facecolor=BG)
ax_cm = fig.add_axes([0.06, 0.50, 0.26, 0.29])
ax_st = fig.add_axes([0.35, 0.50, 0.26, 0.29])
ax_at = fig.add_axes([0.06, 0.08, 0.26, 0.31])
ax_sp = fig.add_axes([0.35, 0.08, 0.26, 0.31])


def heading(y, s):
    fig.text(0.67, y, s, fontsize=8, color=MUTED, weight="bold")


def radio(rect, labels, cb):
    ax = fig.add_axes(rect, facecolor=BG)
    for sp in ax.spines.values():
        sp.set_visible(False)
    rb = RadioButtons(ax, labels, active=0, activecolor=INK)
    for t in rb.labels:
        t.set_fontsize(9)
    rb.on_clicked(cb)
    return rb


state = dict(pi=0, kind=0, dirty=False)
base = lambda: make(PRESETS[state["pi"]])
params = lambda: replace(base(), **{k: s.val for k, s in SL.items()})

heading(0.958, "AIRCRAFT")
heading(0.846, "INPUT")
heading(0.722, "BALANCE AND PROPULSION, APPLIED ON RELEASE")
heading(0.300, "MODEL")

SL, b0 = {}, make(PRESETS[0])
for i, (key, name, step) in enumerate([("Xcg", "CG station", 0.01), ("Sc", "Canard area", 0.05),
                                       ("zT", "Thrust above CG", 0.05), ("Xp", "Prop station", 0.05)]):
    ax = fig.add_axes([0.78, 0.675 - 0.04 * i, 0.17, 0.022], facecolor=GRID)
    lo, hi = BOUNDS[key](b0)
    s = Slider(ax, name, lo, hi, valinit=getattr(b0, key), valstep=step, valfmt="%.2f", initcolor="none",
               handle_style={"facecolor": INK, "edgecolor": INK, "size": 9})
    s.poly.set_visible(False)
    s.label.set_fontsize(9)
    s.valtext.set_fontsize(9)
    s.on_changed(lambda v: state.update(dirty=True))
    SL[key] = s

verdict = fig.text(0.67, 0.515, "", fontsize=12, weight="bold", va="top")
detail = fig.text(0.67, 0.485, "", fontsize=8, family="monospace", va="top", linespacing=1.5)

EQ = [(0.255, r"$X_{np}=\frac{k_cX_c+k_wX_w+k_pX_p}{k_c+k_w+k_p}$"),
      (0.213, r"$k_c=S_ca_c(1+\epsilon_\alpha),\quad k_p=S_pC_{N\alpha}$"),
      (0.183, r"$k_w=S_wa_w\left[1-\epsilon_c(1+\epsilon_\alpha)\right]$"),
      (0.150, r"$K_n=\frac{X_{np}-X_{cg}}{\bar{c}}$"),
      (0.110, r"$M=L_cl_c-L_wl_w+\bar{q}S_w\bar{c}C_{m0}-Tz_T-\bar{q}S_pC_{N\alpha}\,\alpha\,l_p$"),
      (0.075, r"$l_c=X_{cg}-X_c,\quad l_w=X_w-X_{cg},\quad l_p=X_p-X_{cg}$"),
      (0.040, r"$\omega_{sp}\approx\sqrt{\bar{q}\,(k_c+k_w+k_p)(X_{np}-X_{cg})/I_y}$")]
for y, s in EQ:
    fig.text(0.67, y, s, fontsize=10, va="center")


def prep(ax, title, xl):
    ax.clear()
    ax.set_facecolor(PANEL)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(RULE)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel(xl, fontsize=8)


def ruler(p, Xnp, col):
    ax = ax_bal
    ax.clear()
    ax.set_facecolor(BG)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    ax.set_yticks([])
    ax.set_xlim(0, max(p.Xp, p.Xw) + 1.5)
    ax.set_ylim(-1.6, 1.5)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.axhline(0, color=INK, lw=1.2)
    x0, x1 = sorted((Xnp, p.Xcg))
    ax.fill_between([x0, x1], -0.5, 0.5, color=col, alpha=0.35, lw=0)
    for n, X, c in (("canard", p.Xc, CAN), ("wing", p.Xw, WING), ("prop", p.Xp, INK)):
        ax.plot([X, X], [-0.5, 0.5], color=c, lw=3)
        ax.text(X, 0.7, n, color=c, ha="center", fontsize=9)
    st = Xnp > p.Xcg
    ax.plot(p.Xcg, 0, "o", mfc=PANEL, mec=INK, ms=8, zorder=3)
    ax.plot(Xnp, 0, "^", color=INK, ms=7, zorder=3)
    ax.text(p.Xcg + (-0.15 if st else 0.15), -0.95, "CG", ha="right" if st else "left", fontsize=9)
    ax.text(Xnp + (0.15 if st else -0.15), -0.95, "NP", ha="left" if st else "right", fontsize=9)


def refresh(*_):
    pr = PRESETS[state["pi"]]
    p, V0, h0 = params(), pr["V0"], pr["h0"]
    rho = rho_at(h0)
    Xnp, kc, kw, kp = neutral_point(p)
    Kn = (Xnp - p.Xcg) / p.cbar
    col, word = (GOOD, "Stable") if Kn >= 0.05 else (WARN, "Barely stable") if Kn >= 0 else (BAD, "Unstable")
    ok, z = trim_any(p, V0, h0)

    ruler(p, Xnp, col)
    prep(ax_cm, "pitching moment coefficient against angle of attack", "angle of attack (deg)")
    prep(ax_st, "approach to stall at 1 g, % of stall angle", "speed (m/s)")
    prep(ax_at, "pitch attitude (solid) and angle of attack (dashed), deg", "time (s)")
    prep(ax_sp, "speed, m/s", "time (s)")
    info = [f"NP {Xnp:.2f} m   CG {p.Xcg:.2f} m   margin {Kn * p.cbar:.2f} m"]

    if not ok:
        info.append(f"no level-flight trim at {V0:.0f} m/s")
    else:
        al, dc, T = z
        tau = T / thrust_max(p, V0)
        qb = 0.5 * rho * V0**2
        ag = np.arange(-4, 24.01, 0.25)
        A = [aero(p, V0, a * DEG, 0.0, dc, T, rho) for a in ag]
        cm = np.array([x[2] for x in A]) / (qb * p.Sw * p.cbar)
        fc = np.abs([x[5] for x in A]) * p.Ac / p.clc
        fw = np.abs([x[6] for x in A]) * p.Aw / p.clw
        ax_cm.plot(ag, cm, color=INK, lw=1.6)
        ax_cm.axhline(0, color=MUTED, lw=0.8)
        for f_, c_, n_ in ((fc, CAN, "canard stall"), (fw, WING, "wing stall")):
            i = int(np.argmax(f_ >= 1))
            if f_[i] >= 1:
                ax_cm.axvline(ag[i], color=c_, ls="--", lw=1.3, label=n_)
        ax_cm.plot(al / DEG, 0, "o", color=INK, ms=6)
        if ax_cm.get_legend_handles_labels()[0]:
            ax_cm.legend(frameon=False, fontsize=8)

        _, _, _, Lc, Lw, _, _ = aero(p, V0, al, 0.0, dc, T, rho)
        info += [f"trim {V0:.0f} m/s: alpha {al / DEG:.1f} deg",
                 f"canard {dc / DEG:+.1f} deg of +-{p.dcMax / DEG:.0f}",
                 f"power {tau * 100:.0f}%" + ("  ABOVE AVAILABLE" if tau > 1 else ""),
                 f"canard carries {100 * Lc / max(Lc + Lw, 1):.0f}% of lift"]

        ev = np.linalg.eigvals(linearize(p, V0, h0, z))
        cx = sorted((e for e in ev if e.imag > 1e-6), key=lambda e: -abs(e))
        names = ("short period", "phugoid") if len(cx) == 2 else ("oscillation",)
        for n_, e in zip(names, cx):
            info.append(f"{n_}: wn {abs(e):.2f} rad/s, zeta {-e.real / abs(e):.2f}")
        re_ = sorted((e.real for e in ev if abs(e.imag) <= 1e-6), reverse=True)
        if re_:
            info.append("real roots: " + ", ".join(f"{r:+.2f}" for r in re_) + " /s")
        if ev.real.max() > 1e-9:
            info.append(f"divergent, doubling time {math.log(2) / ev.real.max():.1f} s")

        rows, end = stall_sweep(p, V0, h0, z)
        if len(rows):
            Vs, fcl, fwl = rows[-1]
            if end == "stall":
                first, other = ("canard", fwl) if fcl >= fwl else ("wing", fcl)
                info.append(f"{first} stalls first near {Vs:.0f} m/s, the other at {other * 100:.0f}%")
            elif end == "authority":
                info.append(f"canard travel runs out near {Vs:.0f} m/s, before a stall")
            else:
                info.append(f"no stall reached, trim ends near {Vs:.0f} m/s")
            ax_st.plot(rows[:, 0], rows[:, 1] * 100, color=CAN, lw=1.8, label="canard")
            ax_st.plot(rows[:, 0], rows[:, 2] * 100, color=WING, lw=1.8, label="wing")
            ax_st.axhline(100, color=BAD, ls="--", lw=1.2)
            ax_st.set_ylim(0, 115)
            ax_st.invert_xaxis()
            ax_st.legend(frameon=False, fontsize=8, loc="lower left")

        S, note = simulate(p, V0, h0, z, state["kind"])
        ax_at.plot(S[:, 0], S[:, 4] / DEG, color=INK, lw=1.6)
        ax_at.plot(S[:, 0], S[:, 2] / DEG, color=INK, ls="--", lw=1.3)
        ax_sp.plot(S[:, 0], S[:, 1], color=INK, lw=1.6)
        if note:
            for a_ in (ax_at, ax_sp):
                a_.text(0.98, 0.04, "ended: " + note, transform=a_.transAxes, ha="right", color=BAD, fontsize=9)

    verdict.set_text(f"{word}, {Kn * 100:.1f}% MAC margin")
    verdict.set_color(col)
    detail.set_text("\n".join(info))
    fig.canvas.draw_idle()


def on_aircraft(label):
    state["pi"] = [pr["name"] for pr in PRESETS].index(label)
    b = base()
    for k, s in SL.items():
        lo, hi = BOUNDS[k](b)
        s.eventson = False
        s.valmin, s.valmax = lo, hi
        s.ax.set_xlim(lo, hi)
        s.set_val(getattr(b, k))
        s.eventson = True
    refresh()


def on_kind(label):
    state["kind"] = KINDS.index(label)
    refresh()


def on_release(_):
    if state["dirty"]:
        state["dirty"] = False
        refresh()


rb_ac = radio([0.665, 0.865, 0.30, 0.088], [pr["name"] for pr in PRESETS], on_aircraft)
rb_in = radio([0.665, 0.752, 0.30, 0.088], KINDS, on_kind)
fig.canvas.mpl_connect("button_release_event", on_release)
refresh()
plt.show()
