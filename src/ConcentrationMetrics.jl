module ConcentrationMetrics

using Random: AbstractRNG, default_rng
using Distributions: Pareto, rand
using Statistics: mean

using ..VasicekModel: irb_capital, wcdr, rho_corporate

export hhi, effective_n, pareto_portfolio, granularity_adjustment, ga_fraction_of_irb

"""
    hhi(exposures)

Herfindahl-Hirschman index of a vector of exposures, computed as
HHI = Σ (EAD_i / Σ EAD_j)². Returns a value in [0, 1]; 1/N for a homogeneous
portfolio of N exposures; 1 for a single-name concentration.
"""
function hhi(exposures::AbstractVector{<:Real})
    total = sum(exposures)
    total > 0 || error("portfolio has zero exposure")
    s = exposures ./ total
    return sum(abs2, s)
end

"""
    effective_n(exposures)

Effective number of obligors (1/HHI).
"""
effective_n(exposures::AbstractVector{<:Real}) = 1 / hhi(exposures)

"""
    pareto_portfolio(n, total; shape=1.5, rng=default_rng())

Synthetic portfolio of `n` exposures, total notional `total`, with sizes drawn
from a Pareto distribution of shape parameter `shape`. Lower shape → heavier
concentration. Default `shape=1.5` produces a realistic medium-concentration
profile where the top-1 exposure is roughly 2–5% of the total in N = 300.
"""
function pareto_portfolio(n::Integer, total::Real; shape::Real=1.5,
                          rng::AbstractRNG=default_rng())
    raw = rand(rng, Pareto(shape), n)
    scale = total / sum(raw)
    return raw .* scale
end

"""
    granularity_adjustment(exposures, pds, lgds; alpha=0.999, rho_fn=VasicekModel.rho_corporate)

First-order granularity adjustment to the ASRF capital, following
Gordy & Lütkebohmert (2013). The contribution per exposure is proportional to
the squared weight sᵢ² times a function of PD, LGD and ρ; sum normalised by
total exposure. The formula used here is the standard closed form approximation:

    GA ≈ (1 / 2K_IRB) * Σ sᵢ² · V_i(pd, lgd, ρ)

where V_i is the variance of the conditional loss and K_IRB is the asymptotic
IRB capital. Returns the absolute GA (fraction of total EAD).
"""
function granularity_adjustment(exposures::AbstractVector{<:Real},
                                 pds::AbstractVector{<:Real},
                                 lgds::AbstractVector{<:Real};
                                 alpha::Real=0.999,
                                 rho_fn=nothing)
    length(exposures) == length(pds) == length(lgds) ||
        error("exposures, pds and lgds must have equal length")
    rho_fn = rho_fn === nothing ? rho_corporate : rho_fn

    total = sum(exposures)
    weights = exposures ./ total

    # asymptotic IRB capital (portfolio-level, weighted average)
    k_asymp = sum(weights[i] * irb_capital(pds[i], lgds[i], rho_fn(pds[i]), alpha)
                  for i in eachindex(pds))

    # unconditional loss variance per exposure (Bernoulli with loss lgd)
    # V_i = lgd² · PD · (1 - PD)
    variances = @. lgds^2 * pds * (1 - pds)

    ga_raw = sum(weights[i]^2 * variances[i] for i in eachindex(weights)) / (2 * max(k_asymp, 1e-12))
    return (ga_absolute = ga_raw, irb_asymptotic = k_asymp)
end

"""
    ga_fraction_of_irb(exposures, pds, lgds; alpha=0.999)

Granularity adjustment expressed as a fraction of the asymptotic IRB capital.
Useful for reporting the GA as "+X% on top of IRB".
"""
function ga_fraction_of_irb(exposures::AbstractVector{<:Real},
                             pds::AbstractVector{<:Real},
                             lgds::AbstractVector{<:Real};
                             alpha::Real=0.999)
    res = granularity_adjustment(exposures, pds, lgds; alpha=alpha)
    return res.ga_absolute / res.irb_asymptotic
end

end # module
