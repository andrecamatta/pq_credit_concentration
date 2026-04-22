# 03 — Multi-fator vs ASRF via Monte Carlo
#
# Simula um modelo multi-fator gaussiano com 3 setores e correlações intra e
# inter conhecidas, e mede o capital (VaR 99,9% menos perda esperada) por
# simulação direta. Compara com a fórmula IRB de um fator (ASRF) sob a mesma
# ρ_intra para identificar o direcionamento de Pykhtin (2004):
#
#   - Em carteira balanceada, a estrutura multi-fator captura diversificação
#     entre setores que o ASRF ignora → capital cai.
#   - Em carteira concentrada num setor, a diversificação entre setores deixa
#     de proteger → capital sobe aproximando-se ou ultrapassando o ASRF.
#
# Uso: julia --project=. examples/03_pykhtin_multifator.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "CapitalEconomico.jl"))
using .CapitalEconomico

using Random
using Plots
using Printf
using DataFrames
using Statistics

const FIGDIR = joinpath(@__DIR__, "figs")
mkpath(FIGDIR)

Random.seed!(20260421)

println("=== 03 — Multi-fator vs ASRF via Monte Carlo ===\n")

pd        = 0.015
lgd       = 0.40
alpha     = 0.999
rho_intra = 0.20
rho_inter = 0.08
n_sim     = 40_000
n_per     = 400

sector_names = ["Indústria", "Comércio/Serviços", "Imóveis"]

scenarios = [
    ("Balanceada",      [1/3, 1/3, 1/3]),
    ("Moderada",        [0.50, 0.30, 0.20]),
    ("Concentrada",     [0.70, 0.20, 0.10]),
    ("Muito concentrada", [0.85, 0.10, 0.05]),
]

println(@sprintf("PD=%.1f%%, LGD=%.0f%%, α=%.2f%%, ρ_intra=%.0f%%, ρ_inter=%.0f%%",
                 100pd, 100lgd, 100alpha, 100rho_intra, 100rho_inter))
println(@sprintf("Monte Carlo: %d simulações, %d devedores por setor\n", n_sim, n_per))

# ASRF baseline: o que a fórmula IRB de um fator devolveria assumindo ρ_intra uniforme
port_ref = build_sector_portfolio([1/3, 1/3, 1/3]; pd = pd, lgd = lgd)
k_asrf = asrf_capital(port_ref, rho_intra, alpha)
@printf("Capital ASRF (1-fator, ρ=%.0f%% uniforme) = %.3f%% do EAD\n\n",
        100rho_intra, 100 * k_asrf)

fields = [("Cenário", 22, nothing), ("Pesos", 22, nothing),
          ("MC capital (%)", 16, "%.3f"), ("vs ASRF", 12, "%+.1f%%")]
print("  ")
for (name, w, _) in fields
    print(rpad(name, w))
end
println()
println("  ", "─"^sum(f[2] for f in fields))

results = DataFrame(name = String[], weights = String[],
                    mc_capital_pct = Float64[], delta_pct = Float64[])

for (name, weights) in scenarios
    port = build_sector_portfolio(weights; pd = pd, lgd = lgd)
    mc   = multifactor_mc_capital(port; rho_intra = rho_intra,
                                   rho_inter = rho_inter,
                                   alpha = alpha, n_sim = n_sim,
                                   n_per_sector = n_per)
    delta = 100 * (mc.capital / k_asrf - 1)
    push!(results, (name,
                    join(round.(weights .* 100; digits=0), "/") * "%",
                    100 * mc.capital, delta))
    @printf("  %-22s%-22s%-16s%-12s\n",
            name,
            join(round.(weights .* 100; digits=0), "/") * "%",
            @sprintf("%.3f", 100 * mc.capital),
            @sprintf("%+.1f%%", delta))
end

# -----------------------------------------------------------------------------
# Sensibilidade ao ρ_inter na carteira concentrada
# -----------------------------------------------------------------------------
println("\nSensibilidade ao ρ_inter (carteira concentrada 70/20/10):")
port_conc = build_sector_portfolio([0.70, 0.20, 0.10]; pd = pd, lgd = lgd)
rho_inters = [0.00, 0.05, 0.08, 0.12, 0.18]

rows_inter = DataFrame(rho_inter = Float64[], mc_capital_pct = Float64[],
                        delta_pct = Float64[])

for ρi in rho_inters
    mc = multifactor_mc_capital(port_conc; rho_intra = rho_intra,
                                 rho_inter = ρi, alpha = alpha,
                                 n_sim = n_sim, n_per_sector = n_per)
    delta = 100 * (mc.capital / k_asrf - 1)
    push!(rows_inter, (ρi, 100 * mc.capital, delta))
    @printf("  ρ_inter = %.2f  →  capital MC = %.3f%%  (vs ASRF %+.1f%%)\n",
            ρi, 100 * mc.capital, delta)
end

# -----------------------------------------------------------------------------
# Gráficos
# -----------------------------------------------------------------------------
plt_conc = bar(
    results.name, results.mc_capital_pct,
    xlabel = "Cenário de concentração setorial",
    ylabel = "Capital (% do EAD)",
    title  = "Capital multi-fator (MC) vs ASRF por nível de concentração",
    label  = "Multi-fator MC",
    size = (900, 520),
    fillcolor = :orange,
    fillalpha = 0.75,
    legend = :bottomright,
    ylim = (0, 1.15 * 100 * k_asrf),
)
hline!(plt_conc, [100 * k_asrf],
       label = "ASRF (correl. uniforme 20%)", linestyle = :dash, color = :navy, linewidth = 2)
savefig(plt_conc, joinpath(FIGDIR, "pykhtin_comparacao.png"))
println("\n✓ figs/pykhtin_comparacao.png")

plt_inter = plot(
    rho_inters, rows_inter.delta_pct,
    xlabel = "ρ inter-setor",
    ylabel = "Delta multi-fator MC vs ASRF (+%)",
    title  = "Carteira 70/20/10 — sensibilidade ao ρ inter-setor",
    linewidth = 2, marker = :circle,
    legend = false,
    size = (900, 500),
)
hline!(plt_inter, [0.0], linestyle = :dash, color = :gray)
savefig(plt_inter, joinpath(FIGDIR, "pykhtin_sensibilidade_inter.png"))
println("✓ figs/pykhtin_sensibilidade_inter.png")

println("\n✓ Experimento 03 concluído.")
