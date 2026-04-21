# 01 — Vasicek analítico e sensibilidades
#
# Reproduz a fórmula-limite de Vasicek (2002) e explora dependências paramétricas
# que Bolder não tabela explicitamente: densidade de perda para diferentes ρ,
# WCDR em função de α, e capital IRB numa matriz 4×4 de PD × ρ.
#
# Uso: julia --project=. examples/01_vasicek_analitico.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "CapitalEconomico.jl"))
using .CapitalEconomico

using Plots
using Printf

const FIGDIR = joinpath(@__DIR__, "figs")
mkpath(FIGDIR)

println("=== 01 — Vasicek analítico e sensibilidades ===\n")

# -----------------------------------------------------------------------------
# Plot 1 — Densidade de perda para ρ ∈ {0.08, 0.12, 0.20, 0.24}, PD fixo em 1%
# -----------------------------------------------------------------------------
pd = 0.01
rhos = [0.08, 0.12, 0.20, 0.24]
xs = 0.001:0.001:0.15

plt_density = plot(
    xlabel = "Taxa de perda L",
    ylabel = "Densidade",
    title  = "Distribuição-limite de Vasicek  (PD = 1%)",
    legend = :topright,
    size = (800, 500),
)
for ρ in rhos
    ys = [vasicek_pdf(x, pd, ρ) for x in xs]
    plot!(plt_density, xs, ys, label = "ρ = $(Int(round(ρ*100)))%", linewidth = 2)
end
savefig(plt_density, joinpath(FIGDIR, "vasicek_densidade.png"))
println("✓ figs/vasicek_densidade.png")

# -----------------------------------------------------------------------------
# Plot 2 — WCDR em função de α, para três níveis de PD, ρ_corporate(PD)
# -----------------------------------------------------------------------------
αs = 0.990:0.0005:0.9995
pds_to_plot = [0.005, 0.01, 0.03]

plt_wcdr = plot(
    xlabel = "Nível de confiança α",
    ylabel = "WCDR (taxa de default no cenário α)",
    title  = "Worst-case default rate em função de α",
    legend = :topleft,
    size = (800, 500),
)
for pd_i in pds_to_plot
    ρ_i = rho_corporate(pd_i)
    ys = [wcdr(α, pd_i, ρ_i) for α in αs]
    plot!(plt_wcdr, αs, ys, label = "PD = $(round(pd_i*100, digits=1))%", linewidth = 2)
end
savefig(plt_wcdr, joinpath(FIGDIR, "vasicek_wcdr_alpha.png"))
println("✓ figs/vasicek_wcdr_alpha.png\n")

# -----------------------------------------------------------------------------
# Tabela 4×4 — Capital IRB (%) para PD × ρ, LGD 40%, α 99.9%, M = 2.5
# -----------------------------------------------------------------------------
pd_grid  = [0.003, 0.010, 0.030, 0.070]
rho_grid = [0.08, 0.12, 0.18, 0.24]
lgd = 0.40
α   = 0.999
M   = 2.5

cap_matrix = Array{Float64}(undef, length(pd_grid), length(rho_grid))
for (i, pd_i) in enumerate(pd_grid), (j, ρ_j) in enumerate(rho_grid)
    cap_matrix[i, j] = 100 * irb_capital(pd_i, lgd, ρ_j, α; M=M)
end

println("Capital IRB (% do EAD) — LGD 40%, α 99.9%, M = 2.5\n")
print("  ", rpad("PD \\ ρ", 12))
for ρ in rho_grid
    print(rpad(@sprintf("ρ=%.0f%%", ρ*100), 12))
end
println()
println("  ", "─"^(12 + 12*length(rho_grid)))
for (i, pd_i) in enumerate(pd_grid)
    print("  ", rpad(@sprintf("PD=%.1f%%", pd_i*100), 12))
    for j in eachindex(rho_grid)
        print(rpad(@sprintf("%.2f%%", cap_matrix[i,j]), 12))
    end
    println()
end
println()

# Sanity checks
println("── Sanity checks")
@printf("  WCDR(0.999, 0.01, 0.12)            = %.5f  (~ 9%% esperado)\n",
        wcdr(0.999, 0.01, 0.12))
@printf("  IRB capital(PD=0.5%%, LGD=40%%, ρ=12%%, α=99.9%%, M=1) = %.3f%%\n",
        100 * irb_capital(0.005, 0.40, 0.12, 0.999; M=1.0))
@printf("  IRB capital(PD=1%%, LGD=45%%, ρ=corp(1%%), α=99.9%%, M=2.5) = %.3f%%  (deve estar em 5–10%%)\n",
        100 * irb_capital(0.01, 0.45, rho_corporate(0.01), 0.999; M=2.5))

println("\n✓ Experimento 01 concluído.")
