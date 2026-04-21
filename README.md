# pq_credit_concentration

Experimentos numéricos em Julia sobre risco de concentração em carteiras de crédito. Código de apoio ao artigo *Concentração em crédito: por que diversificar não basta* (Pílulas de Quant, André Camatta).

## Quick start

```bash
git clone https://github.com/andrecamatta/pq_credit_concentration.git
cd pq_credit_concentration
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/01_vasicek_analitico.jl
julia --project=. examples/02_hhi_granularity.jl
julia --project=. examples/03_pykhtin_multifator.jl
```

Requer Julia ≥ 1.10. As figuras geradas vão para `examples/figs/`; as tabelas numéricas vão para stdout e também como arquivo CSV em `examples/figs/`.

## Experimentos

### 01 — Vasicek analítico e sensibilidades

Implementa a distribuição-limite de Vasicek (2002) para uma carteira homogênea infinitamente granular e explora sensibilidades que Bolder não tabela:

- densidade de perda L para quatro valores de ρ (0,08 / 0,12 / 0,20 / 0,24);
- *worst-case default rate* em função do nível de confiança α, de 99,0% a 99,97%;
- matriz 4×4 de capital IRB (%) em função de PD × ρ, com LGD de 40% e maturidade 2,5.

### 02 — HHI e ajuste de granularidade em carteira sintética

Carteira com exposições Pareto-distribuídas (α = 1,5) mimetizando o perfil *corporate* de um banco de porte médio. Computa HHI, número efetivo de nomes (1/HHI) e ajuste de granularidade de primeira ordem de Gordy e Lütkebohmert (2013). Varre carteiras de 50, 100, 300, 1.000 e 5.000 nomes (médias sobre 100 sementes por ponto para estabilidade) e plota a dependência do GA em três perfis de concentração (α = 1,3 / 1,5 / 2,0).

### 03 — Multi-fator via Monte Carlo

Simula um modelo Gaussiano multi-fator com três setores (pesos 1/3, 50/30/20, 70/20/10, 85/10/5), ρ intra = 0,20 e ρ inter = 0,08, e mede o VaR 99,9% do portfólio por simulação direta (40.000 cenários, 400 devedores por setor). Compara com a fórmula IRB de um fator (ASRF) para identificar como a estrutura multi-fator muda o capital à medida que a concentração setorial aumenta.

## Estrutura

```
pq_credit_concentration/
├── Project.toml
├── Manifest.toml
├── README.md
├── src/
│   ├── CapitalEconomico.jl
│   ├── VasicekModel.jl
│   ├── ConcentrationMetrics.jl
│   └── MultiFactorPykhtin.jl
└── examples/
    ├── 01_vasicek_analitico.jl
    ├── 02_hhi_granularity.jl
    ├── 03_pykhtin_multifator.jl
    └── figs/              (output — gerado em execução)
```

## Referências

Bolder, D. J. (2022). *Modelling Economic Capital*. Springer.
Gordy, M. B.; Lütkebohmert, E. (2013). Granularity adjustment for regulatory capital assessment. *IJCB* 9(3):33–71.
Pykhtin, M. (2004). Multi-factor adjustment. *Risk* 17(3):85–90.
Vasicek, O. A. (2002). The distribution of loan portfolio value. *Risk* 15(12):160–162.
