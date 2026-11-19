# Beyond Tractor Ownership: Modeling Stochastic Farm Power Thresholds, Custom-Hiring Transitions, and Tillage Timeliness in Tanzanian Smallholder Agriculture

[![Python 3.12](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![Author](https://img.shields.io/badge/Author-Kadeghe%20G.%20Fue-darkblue.svg)]()
[![Target Journal](https://img.shields.io/badge/Target%20Journal-Cogent%20Food%20%26%20Agriculture%20(T%26F)-orange.svg)](https://www.tandfonline.com/toc/oafa20/current)
[![Pipeline Status](https://img.shields.io/badge/Pipeline-Complete%20%26%20Verified-success.svg)]()
[![Microdata Sample](https://img.shields.io/badge/Microdata-N%20%3D%2046%2C493%20Holdings%20(3%20Waves)-purple.svg)]()

---

> [!IMPORTANT]
> **Author & Principal Investigator:** This research methodology, stochastic engineering simulation suite, and microeconometric pipeline were conceived, developed, and conducted by **Kadeghe G. Fue** ([Department of Agricultural Engineering, School of Engineering and Technology](https://www.sua.ac.tz), Sokoine University of Agriculture (SUA), Morogoro, Tanzania; Email: `kadefue@sua.ac.tz`).

---

## Executive Summary

This repository contains the reproducible computational research pipeline (`modeling.py`) for the study:

> **"Beyond Tractor Ownership: Modeling Stochastic Farm Power Thresholds, Custom-Hiring Transitions, and Tillage Timeliness in Tanzanian Smallholder Agriculture"**  
> *Author:* Kadeghe G. Fue  
> *Target Journal:* **Cogent Food & Agriculture** (Taylor & Francis, Section: *Land & Farm Management*)
> *Submission Status:* Prepared for submission to *Cogent Food & Agriculture*

The study pioneers a bidirectional coupling between **agricultural machinery systems engineering** ([ASABE Standards D497.7 and EP496.3](https://elibrary.asabe.org)) and **microeconometric panel econometrics** ([Deaton 1985](https://doi.org/10.1016/0304-4076(85)90059-3); [Papke & Wooldridge 1996](https://doi.org/10.1002/(SICI)1099-1255(199611)11:6<619::AID-JAE418>3.0.CO;2-1); [Mundlak 1978](https://doi.org/10.2307/1913649); [Hansen 2000](https://doi.org/10.1111/1468-0262.00124)).

Using three consecutive waves of the Tanzania Annual Agricultural Sample Survey (**AASS 2022/23, 2023/24, and 2024/25**) spanning **46,493 cultivating smallholder holdings** (`46,489` analytical holdings after area QA filtering) across all **31 administrative regions** (26 Mainland + 5 Zanzibar), the framework resolves the stochastic mismatch between operational farm motive power, narrow rainfall-dependent tillage windows, and heavy edaphic draft resistance.

```text
       ┌─────────────────────────────────────────────────────────────────────────────────────────┐
       │                           SYSTEMS MODELING ARCHITECTURE                                 │
       ├───────────────────────────┬───────────────────────────────┬─────────────────────────────┤
       │  1. INPUTS & CONSTRAINTS  │   2. COUPLED ENGINE           │  3. OUTCOMES & POLICY       │
       ├───────────────────────────┼───────────────────────────────┼─────────────────────────────┤
       │ • 3-Wave AASS (N=46,493)  │ • 1,000-Draw Monte Carlo      │ • Median EPD: 0.315 kW/ha   │
       │ • Masika/Vuli Rain Windows│ • ASABE D497.7 / EP496.3      │ • 62.8% Below 0.5 kW/ha     │
       │ • Edaphic Draft Regimes   │ • Deaton Synthetic Cohorts    │ • Tipping Point: 0.20–0.45  │
       │ • MRDI Breakdown & Downtime• Papke-Wooldridge Frac Logit  │ • +0.1405 ASR per kW/ha     │
       │ • Custom-Hiring Queues    │ • Hansen Threshold Estimation │ • Multi-Village Hub Dispatch│
       └───────────────────────────┴───────────────────────────────┴─────────────────────────────┘
```

---

## Abstract

### 1. Context
Agricultural mechanization across Sub-Saharan African farming systems is severely constrained by the stochastic mismatch between operational farm power capacity, narrow rainfall-dependent tillage windows, and workability limits on heavy soils. While tractor adoption is rising through commercial custom-hiring markets, equipment breakdowns and service queueing frequently delay primary land preparation, leading to seasonal crop abandonment.

### 2. Objective
This study quantifies Effective Power Density ($\text{EPD}$, $\text{kW}\cdot\text{ha}^{-1}$) across diverse farm power portfolios in Tanzania, identifies critical mechanization tipping points required to relieve tillage timeliness bottlenecks, and evaluates institutional machinery dispatch interventions to close seasonal harvesting gaps.

### 3. Methods
We develop a coupled engineering–econometric framework integrating ASABE machinery management standards with a 1,000-draw Monte Carlo simulation engine and microeconometric pseudo-panel validation. Using three consecutive waves of the Tanzania Annual Agricultural Sample Survey (AASS 2022/23–2024/25; $N = 46,493$ harmonized holdings; $46,489$ analytical holdings across 31 regions), we simulate $\text{EPD}$ across human, draught animal, two-wheel tractor (2WT), and four-wheel tractor (4WT) tiers, adjusting for equipment breakdown downtime (MRDI) and custom-hiring queueing allocation ($\phi$). We introduce the Area Survival Ratio ($\text{ASR}$, harvested-to-planted area) to measure timeliness bottleneck relief, validating system elasticities and power thresholds via Deaton synthetic cohort fixed-effects, Papke-Wooldridge fractional logit, Mundlak correlated random effects, and Hansen threshold regressions.

### 4. Results and Conclusions
National weighted mean $\text{EPD}$ stands at $2.054\ \text{kW}\cdot\text{ha}^{-1}$ (median $0.315\ \text{kW}\cdot\text{ha}^{-1}$, SD $7.726\ \text{kW}\cdot\text{ha}^{-1}$), leaving $62.8\%$ of holdings ($60.9\%\text{--}71.7\%$ across smallholder tiers) below the $0.5\ \text{kW}\cdot\text{ha}^{-1}$ mechanization threshold, while commercial custom-hiring supplies $69.5\%\text{--}81.4\%$ of mechanical power in commercialized zones. Econometric validation confirms a highly significant $\text{EPD}$ elasticity of $+0.1405$ on $\text{ASR}$ ($p < 0.001$) in the Deaton pseudo-panel ($N = 379$ cohort-wave observations) and reveals structural mechanization tipping points at $\gamma^* = 0.20\text{--}0.45\ \text{kW}\cdot\text{ha}^{-1}$. Heavy vertisols exhibit an acute $1.304\ \text{kW}\cdot\text{ha}^{-1}$ power deficit ($86.8\%$ of holdings in deficit), where draught animal and 2WT capacities prove insufficient to avert timeliness-induced abandonment during abbreviated planting windows.

### 5. Significance
Coupling ASABE engineering simulations with microeconometric pseudo-panels demonstrates that custom-hiring markets serve as the primary engine of Sub-Saharan mechanization transitions. Achieving food system resilience requires moving beyond uncoordinated individual tractor ownership toward multi-village machinery dispatch hubs and specialized implement financing to break timeliness bottlenecks on difficult soils.

---

## Theoretical & Methodological Overview

Rather than reproducing full mathematical derivations here, this section summarizes the conceptual logic, engineering priors, and econometric estimators implemented in [`modeling.py`](modeling.py). Readers seeking exact equations and implementation details can inspect the referenced code sections and primary literature below.

### 1. Stochastic Farm Power Density (EPD) & Operational Readiness

For each cultivating holding, **Effective Power Density (EPD, in kW/ha)** measures total operational motive power available across four energy tiers divided by effective cultivated area (floored at `0.05 ha` to prevent small-plot division artifacts). Uncertainty in machine rated power and animal draft output is propagated through a **1,000-draw Monte Carlo simulation** (`StochasticPowerEngine` in [`modeling.py`](modeling.py#L760-L828)):

| Power Tier | Source & Ownership | Rated Power Distribution | Tractive / Work Efficiency (`ETA`) | Hiring Allocation (`PHI`) or Breakdown Rate (`P_BROKEN`) |
|:---|:---|:---|:---|:---|
| **Human Labor** | Working-age adults (15–64 yrs) active in agriculture | `0.075 kW` per worker | `0.85` | Floor of 1 active worker per cultivating holding |
| **Draught Animals** | Owned ox-pairs (up to 3 pairs) | Normal `N(0.75, 0.04²)` kW/pair | `0.70` | `6%` breakdown/health loss (`0.30` allocation if hired) |
| **Two-Wheel Tractor (2WT)** | Owned or custom-hired power tillers | Uniform `U(7.5, 11.0)` kW | `0.68` | `28%` breakdown rate if owned; `0.22` allocation if hired |
| **Four-Wheel Tractor (4WT)** | Owned or custom-hired tractors | Uniform `U(35.0, 60.0)` kW | `0.72` | `22%` breakdown rate if owned; `0.15` allocation if hired |

- **Machinery Reliability Degradation Index (MRDI):** Captures the fraction of owned machinery or draught animals that are operational during the tillage window (`working / total owned`).
- **Custom-Hiring Allocation (`PHI`):** Discounts hired machinery power to reflect service rationing, road transit time, and queuing delays during synchronized planting dates.
- **Further Reading:** [ASABE Standards D497.7 & EP496.3](https://elibrary.asabe.org) (Agricultural Machinery Management); [`StochasticPowerEngine`](modeling.py#L760-L828) and global priors [`P_RANGE`, `ETA`, `PHI`, `P_BROKEN`](modeling.py#L117-L130).

---

### 2. Edaphic Draft Resistance & Specific Tillage Deficit

Following **ASABE Standards D497.7 / EP496.3**, holdings are classified into three edaphic draft regimes using a crop-by-region rule (`assign_structure()` in [`modeling.py`](modeling.py#L637-L669)):

| Regime | Soil Typology | Draft Resistance (`kN/m²`) | Required Power `P_req` (`kW/ha`) |
|:---|:---|:---|:---|
| **Regime 1** | Vertisol / Mbuga clay | `65.0–95.0` | `1.85` |
| **Regime 2** | Highland loam | `40.0–65.0` | `1.35` |
| **Regime 3** | Upland sandy-loam | `20.0–40.0` | `0.85` |

A holding's **Specific Engineering Power Deficit** (`kW/ha`) is the positive gap between its regime-specific power requirement (`P_req`) and its posterior mean EPD. Holdings with EPD below `P_req` cannot complete primary tillage within optimal soil moisture windows.

---

### 3. Area Survival Ratio (ASR) and Timeliness Bottlenecks

The **Area Survival Ratio (ASR)** is defined as the ratio of harvested area to planted area (bounded between `0` and `1`) for the *Masika* long-rains, *Vuli* short-rains, and pooled annual cropping seasons (`compute_asr()` in [`modeling.py`](modeling.py#L832-L857)). Inadequate motive power forces protracted land preparation past optimal soil moisture windows, resulting in poor seedbed tilth, delayed germination, weed competition, and seasonal crop abandonment.

---

### 4. Microeconometric Validation Suite

To test whether engineering EPD gains translate into higher crop area survival, [`modeling.py`](modeling.py#L863-L1064) estimates four complementary econometric frameworks:

1. **Deaton (1985) Synthetic Pseudo-Panel Fixed Effects:**  
   Aggregates repeated cross-sections into `379` synthetic cohort-wave cells (`Region × Farm Size Tier × Edaphic Draft Regime`, minimum 10 observations per cell) and estimates Weighted Least Squares (WLS) with cohort and wave fixed effects (`build_pseudo_panel()`, `estimate_pseudo_panel_fe()`).  
   *Reference:* Deaton, A. (1985), *Journal of Econometrics*, [10.1016/0304-4076(85)90059-3](https://doi.org/10.1016/0304-4076(85)90059-3).

2. **Papke–Wooldridge (1996) Fractional Logit GLM:**  
   Models bounded `ASR ∈ [0, 1]` using a binomial logit quasi-likelihood with standard errors clustered at the `Region × Wave` level, including edaphic interactions, owned-vs-hired power decompositions, and *Masika*-vs-*Vuli* seasonal splits (`estimate_household_models()`).  
   *Reference:* Papke, L. E., & Wooldridge, J. M. (1996), *Journal of Applied Econometrics*, [10.1002/(SICI)1099-1255(199611)11:6<619::AID-JAE418>3.0.CO;2-1](https://doi.org/10.1002/(SICI)1099-1255(199611)11:6<619::AID-JAE418>3.0.CO;2-1).

3. **Mundlak (1978) Correlated Random Effects (CRE):**  
   Augments the fractional logit specification with cohort-level means of time-varying covariates (EPD, 2WT/4WT power shares, MRDI, and log area) to control for unobserved cohort-level heterogeneity.  
   *Reference:* Mundlak, Y. (1978), *Econometrica*, [10.2307/1913649](https://doi.org/10.2307/1913649).

4. **Hansen (2000) Non-Linear Threshold Estimation:**  
   Performs a 47-point grid search (`0.20–2.50 kW/ha`) using piecewise linear splines stratified by edaphic draft regime to identify the mechanization tipping point `gamma*` where marginal survival gains shift from steep bottleneck relief to diminishing returns (`compute_hansen_thresholds()`).  
   *Reference:* Hansen, B. E. (2000), *Econometrica*, [10.1111/1468-0262.00124](https://doi.org/10.1111/1468-0262.00124).

---

## Architecture of `modeling.py`

[`modeling.py`](modeling.py) (2,217 lines) is a self-contained, single-command reproducible pipeline organized into 10 sequential sections:

| Section | Lines | Core Classes / Functions | Responsibility |
|:---|:---|:---|:---|
| **0. Global Configuration** | `49–159` | `SEED`, `REGIONS`, `CORRIDOR_OF_REGION`, `P_RANGE`, `ETA`, `PHI`, `audit()` | Sets deterministic RNG seed (`202605`), `N_DRAWS = 1000`, region/zone/corridor mappings, ASABE engineering priors, Okabe-Ito color palettes, and logging/audit trail. |
| **1. Normalization & Statistical Helpers** | `160–261` | `wmean()`, `wquantile()`, `wratio()`, `winsor()`, `decode_numeric()`, `crop_group()` | Survey-weighted mean/quantile/ratio utilities, 99th-percentile winsorization, NBS label decoders, and canonical crop grouping. |
| **2. Metadata & Smart Schema Resolver** | `262–633` | `MetadataDictionary`, `WaveSource`, `resolve_wave_dir()`, `read_season()`, `harmonize_wave()`, `soil_validation()` | Parses `data/AASS_metadata/*.csv` dictionaries; discovers `.dta`/`.sav`/`.csv`/`.parquet`/`.xlsx` files per wave; reconciles `holding_id`/`Holding_ID` shifts, acre-to-hectare conversions, household labor rosters, and observed mechanization proxies. |
| **3. Structural Attributes & Imputation** | `634–759` | `assign_structure()`, `eb_rate()`, `StructuralImputer` | Filters `0 < area <= 100 ha` (clipping EPD denominator at `0.05 ha`), assigns farm-size tiers, edaphic draft regimes, and corridors; imputes absent implement rosters via hot-deck donor resampling (ox-pairs), Empirical-Bayes shrunk propensities (custom hiring), literature priors with 2023/24 shed tilt (owned 4WT/2WT), and binomial breakdown rates. |
| **4. Module A — Stochastic Power Engine** | `760–828` | `StochasticPowerEngine` | Runs chunked (`chunk=4000`) 1,000-draw Monte Carlo simulation per holding; computes posterior `epd_mean`, `epd_sd`, quantiles (`q025`–`q975`), deficit probabilities (`p_lt05`, `p_lt15`, `p_lt_req`), analytic source decompositions (`c_hum`, `c_ox_*`, `c_2wt_*`, `c_4wt_*`), dormant power, and MRDI ratios. |
| **5. Module C — Area Survival Ratio** | `829–859` | `compute_asr()` | Computes Masika, Vuli, and pooled `ASR`, `abandoned_ha`, winsorized EPD components, standardized draft resistance `R_c`, and `Region x Wave` cluster IDs. |
| **6. Module D — Econometric Suite** | `860–1066` | `build_pseudo_panel()`, `estimate_pseudo_panel_fe()`, `estimate_household_models()`, `compute_hansen_thresholds()` | Builds `379`-cell Deaton pseudo-panel (`n_obs >= 10`); fits Deaton WLS FE (linear & log-EPD), Pooled OLS, Papke-Wooldridge Fractional Logit (baseline, edaphic interactions, owned vs. hired, Masika vs. Vuli), Mundlak CRE, and Hansen regime-stratified threshold splines. |
| **7. Export Engine — Publication Tables** | `1067–1538` | `save_table()`, `generate_tables()`, `generate_results_json()` | Exports 10 tables (`table_00`–`table_09`) simultaneously to `.csv`, `.tex` (LaTeX `booktabs`), and `.md`, plus the structured summary payload `outputs/modeling_results.json`. |
| **8. Visualization Engine — Figures** | `1539–1912` | `setup_mpl_theme()`, `save_figure()`, `generate_figures()` | Renders 9 publication figures (`fig_01`–`fig_09`) at 600 DPI `.png` and vector `.pdf` using Okabe-Ito palettes. |
| **9. Cartographic Engine — Maps** | `1913–2116` | `load_tanzania_geography()`, `add_cartographic_elements()`, `generate_maps()` | Loads shapefiles from `data/tanzania_maps/`, reprojects to **UTM Zone 36S (`EPSG:32736`)**, and renders 4 maps (`map_01`–`map_04`) with scalebars, north arrows, and Kilombero/Usangu insets. |
| **10. Main Execution & Verification Gate** | `2117–2217` | `main()` | Orchestrates Stages 1–8, verifies deliverable counts (`>= 9` tables, `>= 8` figures, `>= 4` maps) and non-zero file sizes, and logs total execution time (`outputs/modeling.log`). |

---

## Empirical Microdata Architecture

The study ingests and harmonizes three consecutive waves of the Tanzania Annual Agricultural Sample Survey (**AASS**), administered by the National Bureau of Statistics (NBS) and Ministry of Agriculture:

| Survey Wave | Agricultural Year | Harmonized Holdings ($N$) | Analytical Holdings | Regional Coverage | Mechanization Modules |
|:---|:---|:---|:---|:---|:---|
| **Wave 1** | 2022/2023 | 13,973 | 13,973 | 31 Regions | Custom-Hire Services, Buildings/Shed, Crops |
| **Wave 2** | 2023/2024 | 16,599 | 16,595 | 31 Regions | Buildings/Machinery Shed, Crops, Labor |
| **Wave 3** | 2024/2025 | 15,921 | 15,921 | 31 Regions | Livestock/Bulls, Buildings, Crops, Labor |
| **Total Pooled** | **2022–2025** | **46,493** | **46,489** | **31 Regions (All Tanzania)** | **Harmonized 3-Wave Multi-Tier Sample** |

### Spatial Stratification & Corridors
Holdings are stratified into 5 agricultural mechanization corridors (`CORRIDOR_OF_REGION` in `modeling.py`):
1. **C1 Southern Highlands 4WT:** Iringa, Njombe, Mbeya (non-paddy), Songwe, Ruvuma, Rukwa.
2. **C2 Northern Commercial 4WT:** Arusha, Kilimanjaro, Manyara, Tanga.
3. **C3 Lowland Paddy 2WT:** Morogoro (Kilombero Basin), Pwani, Katavi, and Mbeya paddy holdings (Usangu Plains).
4. **C4 Central-Western DAP:** Dodoma, Singida, Tabora, Shinyanga, Simiyu, Mwanza, Geita, Mara.
5. **C5 Hand-Hoe Coastal/Insular:** Kagera, Kigoma, Lindi, Mtwara, Dar es Salaam, and all 5 Zanzibar regions.

---

## Repository Directory Tree

```text
.
├── data/
│   ├── AASS 2022/ & AASS 2022_2023/    # Wave 1 Raw Microdata (13,973 holdings)
│   ├── AASS 2023/ & AASS 2023_2024/    # Wave 2 Raw Microdata (16,599 holdings)
│   ├── AASS 2024/ & AASS 2024_2025/    # Wave 3 Raw Microdata (15,921 holdings)
│   ├── AASS_metadata/                  # NBS Survey Questionnaires, Codebooks & Dictionaries
│   ├── literature/                     # Reference literature & priors
│   └── tanzania_maps/                  # Shapefiles: Regions, Districts, AEZs, Basins (EPSG:32736)
├── outputs/
│   ├── tables/                         # 10 Publication Tables (.csv, .tex [booktabs], .md)
│   ├── figures/                        # 9 Publication Figures (600 DPI .png, vector .pdf)
│   ├── maps/                           # 4 Cartographic Maps (600 DPI .png, vector .pdf)
│   ├── modeling_results.json           # Validated empirical results summary payload
│   └── modeling.log                    # Full execution & verification audit log
├── modeling.py                         # End-to-end simulation & econometric pipeline (2,217 lines)
├── requirements.txt                    # Pinned Python 3.12 dependencies
├── .gitignore                          # Virtualenv & cache exclusions
└── README.md                           # Documentation & architectural blueprint
```

---

## Deliverables & Scientific Catalog

### 1. Publication Tables (`outputs/tables/`)

All 10 tables are exported simultaneously in `.csv`, `.md`, and `.tex` (LaTeX `booktabs`):

| File | Title & Contents | Key Empirical Result |
|:---|:---|:---|
| [`table_00_schema_harmonization_audit`](outputs/tables/table_00_schema_harmonization_audit.md) | 3-Wave cross-walk, schema reconciliation & imputation audit | Documents identifier shifts (`holding_id`, `Holding_ID`), acre-to-ha conversions, and imputation anchors. |
| [`table_01_sample_mechanization_inventory`](outputs/tables/table_01_sample_mechanization_inventory.md) | Sample architecture and machinery inventory across 5 corridors | Owned 4WT prevalence is $0.0\%\text{--}0.8\%$; custom-hired 4WT reaches $27.8\%$ in C2 Northern and $22.8\%$ in C1 Southern Highlands. |
| [`table_02_stochastic_epd_distributions`](outputs/tables/table_02_stochastic_epd_distributions.md) | Monte Carlo $\text{EPD}$ ($\text{kW}\cdot\text{ha}^{-1}$) posterior statistics | National weighted mean is $2.054\ \text{kW}\cdot\text{ha}^{-1}$, median $0.315\ \text{kW}\cdot\text{ha}^{-1}$; $62.8\%$ fall below $0.5\ \text{kW}\cdot\text{ha}^{-1}$. |
| [`table_03_power_source_decomposition`](outputs/tables/table_03_power_source_decomposition.md) | Power supply decomposition by corridor and energy source | Custom-hired 4WT delivers $69.5\%\text{--}81.4\%$ of mechanical power in commercial zones; draught animals supply $26.0\%$ in C4. |
| [`table_04_machinery_reliability_downtime`](outputs/tables/table_04_machinery_reliability_downtime.md) | Machinery Operational Readiness (MRDI) and downtime losses | Breakdown downtime destroys $10.9\%\text{--}31.2\%$ of motive power during critical planting windows. |
| [`table_05_edaphic_draft_tillage_deficit`](outputs/tables/table_05_edaphic_draft_tillage_deficit.md) | Edaphic draft resistance ($R_{\text{draft}}$) vs. specific power deficit | Vertisols require $1.85\ \text{kW}\cdot\text{ha}^{-1}$; observed mean $\text{EPD}$ is $1.116$, yielding a $1.304\ \text{kW}\cdot\text{ha}^{-1}$ deficit across $86.8\%$ of holdings. |
| [`table_06_seasonal_asr_masika_vuli`](outputs/tables/table_06_seasonal_asr_masika_vuli.md) | Planting-to-harvested Area Survival Ratio by season & power tier | Masika long-rains display lower survival ($0.656$) than Vuli short-rains ($0.722$), a $12.8\%$ abandonment gap. |
| [`table_07_pseudopanel_cohort_regressions`](outputs/tables/table_07_pseudopanel_cohort_regressions.md) | Deaton (1985) synthetic cohort WLS regressions ($N=379$) | Linear $\text{EPD}$ elasticity: $\beta = +0.1405$ ($p < 0.001$); draft interaction: $\beta = -0.2364$ ($p < 0.001$). |
| [`table_08_fractional_logit_asr_models`](outputs/tables/table_08_fractional_logit_asr_models.md) | Papke-Wooldridge fractional logit & Mundlak CRE models | 2WT hiring yields $\text{AME} = +0.0030$ ($p < 0.05$); vertisol interaction dampens survival ($\beta = -0.0206$, $p < 0.10$). |
| [`table_09_threshold_transition_matrix`](outputs/tables/table_09_threshold_transition_matrix.md) | Markov cohort transitions & Hansen threshold parameters | $80.8\%$ of cohorts in Tier 1 ($<0.5\ \text{kW}\cdot\text{ha}^{-1}$) remain trapped; $0.0\%$ transition directly to Tier 4. |

---

### 2. Publication Figures (`outputs/figures/`)

Rendered at **600 DPI PNG** and vector **PDF** with Okabe-Ito colorblind-safe palettes:

1. [`fig_01_monte_carlo_epd_raincloud`](outputs/figures/fig_01_monte_carlo_epd_raincloud.png): Empirical distributions of $\text{EPD}$ across 5 corridors with FAO ($1.5\ \text{kW}\cdot\text{ha}^{-1}$) and SSA ($0.5\ \text{kW}\cdot\text{ha}^{-1}$) benchmark lines.
2. [`fig_02_power_stack_transition_alluvial`](outputs/figures/fig_02_power_stack_transition_alluvial.png): Dynamic power stack transitions from manual labor to custom-hired mechanization over 2022–2025.
3. [`fig_03_edaphic_draft_vs_asr_response`](outputs/figures/fig_03_edaphic_draft_vs_asr_response.png): Non-parametric kernel regression showing the collapse in $\text{ASR}$ on high-draft vertisols when power $< 0.5\ \text{kW}\cdot\text{ha}^{-1}$.
4. [`fig_04_masika_vs_vuli_survival_gap`](outputs/figures/fig_04_masika_vs_vuli_survival_gap.png): Seasonal survival disparity showing how mechanized tiers close the $12.8\%$ Masika abandonment gap.
5. [`fig_05_pseudopanel_cohort_trajectories`](outputs/figures/fig_05_pseudopanel_cohort_trajectories.png): Longitudinal trajectory vectors of synthetic cohorts tracing $\text{EPD}$ growth against $\text{ASR}$ gains over 3 waves.
6. [`fig_06_owned_vs_hired_efficiency_frontier`](outputs/figures/fig_06_owned_vs_hired_efficiency_frontier.png): Dual-panel tillage efficiency frontiers comparing private ownership vs. custom-hiring across farm sizes ($0.2\text{--}12.0\ \text{ha}$) and the 14-day soil moisture window.
7. [`fig_07_econometric_forest_marginal_effects`](outputs/figures/fig_07_econometric_forest_marginal_effects.png): Forest plot comparing $\ln(\text{EPD})$ marginal effects across Pooled OLS, Fractional Logit, Mundlak CRE, and Deaton FE.
8. [`fig_08_monte_carlo_sensitivity_tornado`](outputs/figures/fig_08_monte_carlo_sensitivity_tornado.png): Parameter sensitivity decomposition across rated power, tractive efficiency, and custom-hiring allocation priors.
9. [`fig_09_engineering_capacity_reliability_frontier`](outputs/figures/fig_09_engineering_capacity_reliability_frontier.png): Engineering capacity utilization scatter plot displaying $\text{EPD}$ versus $\text{ASR}$ colored by the Machinery Reliability Degradation Index (MRDI).

---

### 3. Cartographic Maps (`outputs/maps/`)

Projected to **UTM Zone 36S (`EPSG:32736`)** with publication scalebars, north arrows, and administrative boundary overlays:

1. [`map_01_effective_power_density_evolution`](outputs/maps/map_01_effective_power_density_evolution.png): Regional $\text{EPD}$ in 2024/25 and net structural change $\Delta\text{EPD}$ (2022–2025) with proportional custom-hiring symbols.
2. [`map_02_dominant_mechanization_typology`](outputs/maps/map_02_dominant_mechanization_typology.png): Spatial typology classifying 31 regions into 4WT commercial, 2WT paddy, draught ox, and hand-hoe zones.
3. [`map_03_bivariate_edaphic_power_deficit`](outputs/maps/map_03_bivariate_edaphic_power_deficit.png): $3 \times 3$ bivariate choropleth identifying acute vulnerability hotspots where heavy vertisol draft collides with low power.
4. [`map_04_area_survival_ratio_hotspots`](outputs/maps/map_04_area_survival_ratio_hotspots.png): Regional crop survival choropleth with detailed insets for the Kilombero River Basin and Usangu Plains.

---

## Reproducibility & Execution Blueprint

### 1. Environment Configuration

```bash
# Create and activate Python 3.12 virtual environment
python3.12 -m venv .venv
source .venv/bin/activate

# Install exact pinned dependencies
pip install -r requirements.txt
```

### 2. Run the Full Modeling & Econometric Pipeline

```bash
./.venv/bin/python modeling.py
```

*Autonomously executes:*
- Ingestion & schema harmonization across `46,493` holdings (Waves 1–3).
- Structural attribute assignment & defensive imputation (`46,489` analytical holdings).
- 1,000-draw Monte Carlo machinery simulation engine.
- Area Survival Ratio ($\text{ASR}$) derivation for Masika and Vuli cropping seasons.
- Microeconometric models: Deaton pseudo-panel WLS, Pooled OLS, Fractional Logit GLM, Mundlak CRE, and Hansen threshold search.
- Compilation of 10 tables, 9 figures, 4 spatial maps, `outputs/modeling_results.json`, and `outputs/modeling.log`.

---

## Verification Audit

Every execution of `modeling.py` concludes with an automated verification gate checking target deliverable counts and non-zero file sizes:

```text
╒════════════════════════╤════════════════╤═════════════════════════╤═══════════════════════╕
│ Deliverable Category   │   Unique Items │   Total Files Generated │ Verification Status   │
╞════════════════════════╪════════════════╪═════════════════════════╪═══════════════════════╡
│ Tables (Target >= 9)   │             10 │                      30 │ PASSED                │
├────────────────────────┼────────────────┼─────────────────────────┼───────────────────────┤
│ Figures (Target >= 8)  │              9 │                      18 │ PASSED                │
├────────────────────────┼────────────────┼─────────────────────────┼───────────────────────┤
│ Maps (Target >= 4)     │              4 │                       8 │ PASSED                │
╘════════════════════════╧════════════════╧═════════════════════════╧═══════════════════════╛
All 56 analytical files verified. Pipeline finished with ZERO unhandled exceptions.
```

---

## Data Availability Statement

This study is prepared for submission to **Cogent Food & Agriculture** (Taylor & Francis) under the title *"Beyond Tractor Ownership: Modeling Stochastic Farm Power Thresholds, Custom-Hiring Transitions, and Tillage Timeliness in Tanzanian Smallholder Agriculture"*.

- **Primary Microdata & Shapefiles:** The empirical analysis uses three consecutive waves of the Tanzania Annual Agricultural Sample Survey (**AASS 2022/23, 2023/24, and 2024/25**; $N = 46,493$ harmonized holdings; $46,489$ analytical holdings) administered by the National Bureau of Statistics (NBS) Tanzania and the Ministry of Agriculture, together with official regional/basin vector shapefiles (`data/tanzania_maps/`). Raw AASS microdata are archived by the NBS Tanzania Microdata Catalog (`[INSERT NBS MICRODATA CATALOG URL / ACCESSION NUMBER]`).
- **Code, Summary Datasets & Reproducibility Archive:** The complete simulation and econometric pipeline ([modeling.py](modeling.py)), pinned dependencies ([requirements.txt](requirements.txt)), metadata dictionaries (`data/AASS_metadata/`), and all derived analytical outputs are archived at `[INSERT REPOSITORY URL, e.g., Zenodo / GitHub]` under Persistent Identifier `[INSERT DOI / ACCESSION NUMBER, e.g., https://doi.org/10.5281/zenodo.XXXXXXX]`.
- **Values Behind Reported Measures (Means, SDs, Regressions):** All summary statistics, distribution quantiles, regression coefficients, standard errors, and transition matrices reported in the manuscript are provided directly in [`outputs/tables/`](outputs/tables/) (`table_00`–`table_09` in `.csv`, `.md`, and `.tex` formats) and [`outputs/modeling_results.json`](outputs/modeling_results.json).
- **Values Used to Build Figures & Maps:** All underlying values plotted in Figures 1–9 ([`outputs/figures/`](outputs/figures/)) and Maps 1–4 ([`outputs/maps/`](outputs/maps/)) are deterministically reproduced from the input microdata and shapefiles by running `modeling.py` with fixed random seed `SEED = 202605`.
- **Points Extracted from Images:** Not applicable; no quantitative data points were extracted from external images for analysis. All spatial analyses use original vector shapefiles projected to UTM Zone 36S (`EPSG:32736`).
- **Variable Descriptions & Codebook:** Complete variable names, types, and value-label codes (including `holding_id`, `weight`, `region_code`, household member age `15–64` and agricultural work status, seasonal planted/harvested areas, crop types, draught bulls, custom-hiring indicators, and machinery shed ownership) are documented in [`data/AASS_metadata/data_dictionary_catalog_52.csv`](data/AASS_metadata/data_dictionary_catalog_52.csv) and audited in [`outputs/tables/table_00_schema_harmonization_audit.md`](outputs/tables/table_00_schema_harmonization_audit.md).

---

## Citation & Working Paper Reference

This study was conceived, developed, and authored by **Kadeghe G. Fue** as a single-author investigation.

When referencing this methodology, empirical dataset, or codebase, please cite:

```bibtex
@article{fue2026mechanization,
  title={Beyond Tractor Ownership: Modeling Stochastic Farm Power Thresholds, Custom-Hiring Transitions, and Tillage Timeliness in Tanzanian Smallholder Agriculture},
  author={Fue, Kadeghe G.},
  year={2026},
  journal={Cogent Food \& Agriculture (Submitted / Preprint)},
  note={Target journal: Cogent Food \& Agriculture (Taylor \& Francis, Land \& Farm Management).}
}
```

---

*Authored and maintained by **Kadeghe G. Fue** ([Department of Agricultural Engineering, School of Engineering and Technology](https://www.sua.ac.tz), Sokoine University of Agriculture (SUA), Morogoro, Tanzania).*
