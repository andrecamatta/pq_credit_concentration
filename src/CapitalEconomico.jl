module CapitalEconomico

include("VasicekModel.jl")
include("ConcentrationMetrics.jl")
include("MultiFactorPykhtin.jl")

using .VasicekModel
using .ConcentrationMetrics
using .MultiFactorPykhtin

# re-export the public API
export vasicek_cdf, vasicek_pdf, wcdr, irb_capital, rho_corporate, maturity_adjustment
export hhi, effective_n, pareto_portfolio, granularity_adjustment, ga_fraction_of_irb
export multifactor_mc_capital, asrf_capital, build_sector_portfolio, SectorPortfolio

end # module
