module MultiFactorPykhtin

using Distributions: Normal, cdf, quantile
using LinearAlgebra: cholesky, Symmetric
using Random: AbstractRNG, default_rng, randn
using Statistics: quantile as stat_quantile

using ..VasicekModel: wcdr

export multifactor_mc_capital, asrf_capital, build_sector_portfolio, SectorPortfolio

const STANDARD_NORMAL = Normal()

"""
    SectorPortfolio(weights, pd, lgd)

Multi-sector portfolio descriptor. `weights` is the fraction of total EAD in
each sector (must sum to 1). `pd` and `lgd` are assumed common across sectors
for simplicity; the framework generalises trivially to heterogeneous values.
"""
struct SectorPortfolio
    weights::Vector{Float64}
    pd::Float64
    lgd::Float64
end

function build_sector_portfolio(weights::AbstractVector{<:Real}; pd::Real=0.01,
                                 lgd::Real=0.4)
    isapprox(sum(weights), 1.0; atol=1e-8) || error("weights must sum to 1")
    return SectorPortfolio(collect(float.(weights)), float(pd), float(lgd))
end

"""
    asrf_capital(portfolio, rho, alpha=0.999)

Single-factor ASRF capital (Basileia II IRB, unexpected loss per EAD) assuming
a uniform asset correlation ρ across all sectors. This is the regulatory
baseline that ignores sectoral structure.
"""
function asrf_capital(portfolio::SectorPortfolio, rho::Real, alpha::Real=0.999)
    w = wcdr(alpha, portfolio.pd, rho)
    return portfolio.lgd * (w - portfolio.pd)
end

"""
    multifactor_mc_capital(portfolio; rho_intra, rho_inter, alpha, n_sim, n_per_sector, rng)

Monte Carlo estimate of the portfolio unexpected-loss capital under a
multi-factor Gaussian model with intra-sector correlation `rho_intra` and
inter-sector correlation `rho_inter`. Simulates `n_sim` scenarios, each with
`n_per_sector` obligors per sector, and returns the α-quantile of the
portfolio loss minus the expected loss.

Setup per scenario s:
  - Sample K sector factors Y_k with Cov(Y_j, Y_k) = rho_inter (and Var = 1).
  - For each obligor i in sector k: X_{ki} = √ρ_intra Y_k + √(1-ρ_intra) ε_{ki}.
  - Default if X_{ki} < Φ⁻¹(pd). Contribution to loss: w_k/n_per_sector * lgd.
  - Portfolio loss L_s = Σ_k w_k * (defaults_k / n_per_sector) * lgd.

The α-quantile of {L_s} minus the expected loss (pd * lgd) is the unexpected
loss capital requirement under the multi-factor model.
"""
function multifactor_mc_capital(portfolio::SectorPortfolio;
                                 rho_intra::Real, rho_inter::Real,
                                 alpha::Real=0.999, n_sim::Integer=50_000,
                                 n_per_sector::Integer=500,
                                 rng::AbstractRNG=default_rng())
    K = length(portfolio.weights)
    threshold = quantile(STANDARD_NORMAL, portfolio.pd)
    sqrt_rho = sqrt(rho_intra)
    sqrt_1mrho = sqrt(1 - rho_intra)

    # Inter-sector factor covariance
    Σ = fill(float(rho_inter), K, K)
    for k in 1:K
        Σ[k, k] = 1.0
    end
    L = cholesky(Symmetric(Σ)).L

    losses = zeros(Float64, n_sim)

    @inbounds for s in 1:n_sim
        # sample correlated sector factors
        Z  = randn(rng, K)
        Y  = L * Z
        total_loss = 0.0
        for k in 1:K
            # conditional default probability in sector k given Y_k
            cond_pd = cdf(STANDARD_NORMAL,
                          (threshold - sqrt_rho * Y[k]) / sqrt_1mrho)
            # simulate binomial defaults among n_per_sector obligors
            defaults = 0
            for _ in 1:n_per_sector
                if rand(rng) < cond_pd
                    defaults += 1
                end
            end
            total_loss += portfolio.weights[k] * (defaults / n_per_sector) * portfolio.lgd
        end
        losses[s] = total_loss
    end

    q = stat_quantile(losses, alpha)
    expected_loss = portfolio.pd * portfolio.lgd
    return (capital = q - expected_loss,
            var_quantile = q,
            expected_loss = expected_loss,
            losses = losses)
end

end # module
