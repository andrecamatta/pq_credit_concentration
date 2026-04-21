# 02 — HHI e ajuste de granularidade em carteira sintética
#
# Constrói uma carteira de 300 nomes com exposições Pareto-distribuídas
# (mimetizando o perfil corporate de um banco médio brasileiro) e aplica o
# ajuste de granularidade de Gordy-Lütkebohmert (2013). Em seguida, varre
# tamanhos de carteira (50, 100, 300, 1000, 5000) para mostrar como o GA
# escala com a granularidade.
#
# Uso: julia --project=. examples/02_hhi_granularity.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "CapitalEconomico.jl"))
using .CapitalEconomico

using Random
using Printf
using Plots
using DataFrames
using Statistics: mean

const FIGDIR = joinpath(@__DIR__, "figs")
mkpath(FIGDIR)

Random.seed!(20260421)

println("=== 02 — HHI e granularity adjustment em carteira sintética ===\n")

# Parâmetros comuns
pd_uniform  = 0.015   # 1,5% PD média corporate
lgd_uniform = 0.40
alpha       = 0.999

# -----------------------------------------------------------------------------
# Perfil-âncora: carteira de 300 exposições Pareto (perfil S2/S3 corporate)
# -----------------------------------------------------------------------------
n_anchor = 300
exposures_anchor = pareto_portfolio(n_anchor, 1000.0; shape = 1.5)
pds_anchor  = fill(pd_uniform, n_anchor)
lgds_anchor = fill(lgd_uniform, n_anchor)

hhi_anchor = hhi(exposures_anchor)
neff_anchor = effective_n(exposures_anchor)
top_share = maximum(exposures_anchor) / sum(exposures_anchor)

println("Carteira-âncora (300 nomes, Pareto shape=1.5):")
@printf("  HHI              = %.4f\n", hhi_anchor)
@printf("  N efetivo (1/HHI) = %.1f\n", neff_anchor)
@printf("  Maior exposição  = %.2f%% do total\n", 100 * top_share)

ga_anchor = granularity_adjustment(exposures_anchor, pds_anchor, lgds_anchor; alpha = alpha)
@printf("  IRB assintótico  = %.3f%% do EAD\n", 100 * ga_anchor.irb_asymptotic)
@printf("  Granularity adj. = %.3f%% do EAD\n", 100 * ga_anchor.ga_absolute)
@printf("  GA como %% do IRB = +%.2f%%\n\n",  100 * ga_anchor.ga_absolute / ga_anchor.irb_asymptotic)

# Contraste com carteira homogênea de mesmo N
exposures_homog = fill(1000.0 / n_anchor, n_anchor)
ga_homog = granularity_adjustment(exposures_homog, pds_anchor, lgds_anchor; alpha = alpha)
@printf("Carteira homogênea de mesmo N = 300 (comparação):\n")
@printf("  HHI              = %.4f   (= 1/N)\n", hhi(exposures_homog))
@printf("  GA como %% do IRB = +%.2f%%\n\n", 100 * ga_homog.ga_absolute / ga_homog.irb_asymptotic)

# -----------------------------------------------------------------------------
# Varredura: mesmas Pareto shape=1.5, varia N ∈ {50, 100, 300, 1000, 5000}
# -----------------------------------------------------------------------------
println("Varredura de granularidade (Pareto shape=1.5, PD=1.5%, LGD=40%, α=99.9%):\n")

ns = [50, 100, 300, 1000, 5000]
n_seeds_table = 100
rows = DataFrame(
    N         = Int[],
    HHI       = Float64[],
    N_eff     = Float64[],
    top_pct   = Float64[],
    IRB_pct   = Float64[],
    GA_abs_pct = Float64[],
    GA_rel_pct = Float64[],
)
for n in ns
    hhi_vals    = Float64[]
    neff_vals   = Float64[]
    top_vals    = Float64[]
    ga_abs_vals = Float64[]
    ga_rel_vals = Float64[]
    local k_irb = 0.0
    for seed in 1:n_seeds_table
        Random.seed!(20260421 + n * 1000 + seed)
        exp_n = pareto_portfolio(n, 1000.0; shape = 1.5)
        pds_n  = fill(pd_uniform,  n)
        lgds_n = fill(lgd_uniform, n)
        ga     = granularity_adjustment(exp_n, pds_n, lgds_n; alpha = alpha)
        push!(hhi_vals,    hhi(exp_n))
        push!(neff_vals,   effective_n(exp_n))
        push!(top_vals,    100 * maximum(exp_n) / sum(exp_n))
        push!(ga_abs_vals, 100 * ga.ga_absolute)
        push!(ga_rel_vals, 100 * ga.ga_absolute / ga.irb_asymptotic)
        k_irb = 100 * ga.irb_asymptotic
    end
    push!(rows, (
        n,
        mean(hhi_vals),
        mean(neff_vals),
        mean(top_vals),
        k_irb,
        mean(ga_abs_vals),
        mean(ga_rel_vals),
    ))
end

# Imprime tabela manualmente
fields = [("N", 6, "%d"), ("HHI", 10, "%.4f"), ("N_ef", 8, "%.1f"),
          ("Maior (%)", 12, "%.2f"), ("IRB (%)", 10, "%.3f"),
          ("GA abs (%)", 12, "%.4f"), ("GA / IRB (+%)", 14, "%.2f")]
print("  ")
for (name, w, _) in fields
    print(rpad(name, w))
end
println()
println("  ", "─"^sum(f[2] for f in fields))
for row in eachrow(rows)
    vals = [row.N, row.HHI, row.N_eff, row.top_pct, row.IRB_pct, row.GA_abs_pct, row.GA_rel_pct]
    print("  ")
    for ((name, w, fmt), v) in zip(fields, vals)
        s = fmt == "%d" ? @sprintf("%d", v) :
            fmt == "%.4f" ? @sprintf("%.4f", v) :
            fmt == "%.1f" ? @sprintf("%.1f", v) :
            fmt == "%.2f" ? @sprintf("%.2f", v) :
            fmt == "%.3f" ? @sprintf("%.3f", v) : string(v)
        print(rpad(s, w))
    end
    println()
end

# Salva CSV para reuso pela tabela do artigo
using CSV
CSV.write(joinpath(FIGDIR, "hhi_sweep.csv"), rows)
println("\n✓ figs/hhi_sweep.csv")

# -----------------------------------------------------------------------------
# Plot — GA como % do IRB em função de N, para shape Pareto ∈ {1.2, 1.5, 2.0}
# -----------------------------------------------------------------------------
ns_plot = [50, 100, 200, 300, 500, 1000, 2000, 5000]
shapes  = [1.3, 1.5, 2.0]
labels  = ["Pareto α=1.3 (cauda pesada)", "Pareto α=1.5 (S2/S3 típico)", "Pareto α=2.0 (mais disperso)"]
n_seeds = 40

plt = plot(
    xlabel = "Número de exposições (N)",
    ylabel = "Granularity Adjustment (+% sobre IRB)",
    title  = "GA em função de N — média de $n_seeds sementes por ponto",
    xscale = :log10,
    legend = :topright,
    size = (900, 550),
)
for (shape, label) in zip(shapes, labels)
    ys = Float64[]
    for n in ns_plot
        samples = Float64[]
        for seed in 1:n_seeds
            Random.seed!(1000 * n + seed)
            e = pareto_portfolio(n, 1000.0; shape = shape)
            pds_e = fill(pd_uniform, n)
            lgds_e = fill(lgd_uniform, n)
            ga = granularity_adjustment(e, pds_e, lgds_e; alpha = alpha)
            push!(samples, 100 * ga.ga_absolute / ga.irb_asymptotic)
        end
        push!(ys, mean(samples))
    end
    plot!(plt, ns_plot, ys, label = label, linewidth = 2, marker = :circle)
end
savefig(plt, joinpath(FIGDIR, "hhi_granularity_sweep.png"))
println("✓ figs/hhi_granularity_sweep.png")

# Plot complementar: Lorenz-like curve (distribuição de exposições) da carteira-âncora
plt2 = plot(
    title  = "Distribuição de exposições — carteira-âncora (N=300, Pareto α=1.5)",
    xlabel = "Nome (ordenado por tamanho decrescente)",
    ylabel = "Exposição (% do total)",
    legend = false,
    size = (900, 500),
)
sorted_exp = sort(exposures_anchor, rev = true)
bar!(plt2, 1:length(sorted_exp), 100 .* sorted_exp ./ sum(sorted_exp),
     linewidth = 0, fillalpha = 0.8)
savefig(plt2, joinpath(FIGDIR, "pareto_portfolio.png"))
println("✓ figs/pareto_portfolio.png")

println("\n✓ Experimento 02 concluído.")
