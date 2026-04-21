module VasicekModel

using Distributions: Normal, cdf, quantile

export vasicek_cdf, vasicek_pdf, wcdr, irb_capital, rho_corporate, maturity_adjustment

const STANDARD_NORMAL = Normal()

"""
    vasicek_cdf(x, pd, rho)

Vasicek (2002) limit CDF of the loss rate for a homogeneous, infinitely granular
portfolio under a one-factor Gaussian model. Returns P(L ≤ x).
"""
function vasicek_cdf(x::Real, pd::Real, rho::Real)
    0.0 < x < 1.0 || return x <= 0.0 ? 0.0 : 1.0
    0.0 < rho < 1.0 || error("rho must be in (0,1)")
    num = sqrt(1 - rho) * quantile(STANDARD_NORMAL, x) - quantile(STANDARD_NORMAL, pd)
    return cdf(STANDARD_NORMAL, num / sqrt(rho))
end

"""
    vasicek_pdf(x, pd, rho)

Density of the Vasicek limit loss distribution.
"""
function vasicek_pdf(x::Real, pd::Real, rho::Real)
    0.0 < x < 1.0 || return 0.0
    Φinv = quantile(STANDARD_NORMAL, x)
    a = quantile(STANDARD_NORMAL, pd)
    k = sqrt(1 - rho) / sqrt(rho)
    z = (sqrt(1 - rho) * Φinv - a) / sqrt(rho)
    φ_Φinv = exp(-0.5 * Φinv^2) / sqrt(2π)
    φ_z    = exp(-0.5 * z^2) / sqrt(2π)
    return (k * φ_z) / φ_Φinv
end

"""
    wcdr(alpha, pd, rho)

Worst-case default rate at confidence α: the α-quantile of the loss distribution
conditional on an extreme systemic outcome. This is the core of the IRB formula.
"""
function wcdr(alpha::Real, pd::Real, rho::Real)
    num = quantile(STANDARD_NORMAL, pd) + sqrt(rho) * quantile(STANDARD_NORMAL, alpha)
    return cdf(STANDARD_NORMAL, num / sqrt(1 - rho))
end

"""
    rho_corporate(pd)

Asset correlation for corporate, sovereign and bank exposures under Basel II IRB
(BCBS 128, para. 272). Varies smoothly from 0.24 for very low PD to 0.12 for high PD.
"""
function rho_corporate(pd::Real)
    0.0 < pd < 1.0 || error("pd must be in (0,1)")
    w = (1 - exp(-50 * pd)) / (1 - exp(-50.0))
    return 0.12 * w + 0.24 * (1 - w)
end

"""
    maturity_adjustment(pd, M)

Maturity adjustment factor (1 + (M - 2.5) b(PD)) / (1 - 1.5 b(PD))
from BCBS 128, para. 272.
"""
function maturity_adjustment(pd::Real, M::Real)
    b = (0.11852 - 0.05478 * log(pd))^2
    return (1 + (M - 2.5) * b) / (1 - 1.5 * b)
end

"""
    irb_capital(pd, lgd, rho, alpha; M=1.0)

Basel II IRB capital requirement (unexpected loss) per unit of EAD.
Default: α = 0.999, M = 1.0. Returns the fraction of EAD to hold as capital.
"""
function irb_capital(pd::Real, lgd::Real, rho::Real, alpha::Real=0.999; M::Real=1.0)
    k = lgd * (wcdr(alpha, pd, rho) - pd)
    return k * maturity_adjustment(pd, M)
end

end # module
