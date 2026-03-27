# Causal Inference for Functional Data: Literature Review

A curated bibliography of methods at the intersection of causal inference
and functional data analysis (FDA), organized by topic.

---

## 1. Treatment Effects with Functional Outcomes (FATE)

**Ecker, K., de Luna, X., & Schelin, L. (2024).**
Causal inference with a functional outcome.
*Journal of the Royal Statistical Society Series C: Applied Statistics*, 73(1), 221--240.
DOI: [10.1093/jrsssc/qlad092](https://doi.org/10.1093/jrsssc/qlad092)

> Defines the Functional Average Treatment Effect (FATE) -- the expected
> difference between potential outcome functions for treated vs. control
> units. Develops an outcome regression estimator with simultaneous
> confidence bands controlling coverage over the full functional domain.
> Applied to urban vs. rural residence effects on lifetime income
> trajectories in Sweden.

**Testa, L., Boschi, T., Chiaramonte, F., Kennedy, E. H., & Reimherr, M. (2025).**
Doubly-robust functional average treatment effect estimation.
*arXiv preprint*, 2501.06024.
DOI: [10.48550/arXiv.2501.06024](https://doi.org/10.48550/arXiv.2501.06024)

> DR-FoS: a doubly robust estimator for FATE that is consistent even if
> either the outcome or treatment assignment model is misspecified.
> Converges to a Gaussian process, enabling valid simultaneous confidence
> bands. Applied to the SHARE health dataset.

**Salmaso, F., Testa, L., & Chiaramonte, F. (2026).**
A doubly robust machine learning approach for disentangling treatment
effect heterogeneity with functional outcomes.
*arXiv preprint*, 2602.11118.

> FOCaL (Functional Outcome Causal Learning): a doubly robust meta-learner
> for functional conditional average treatment effects (F-CATE).
> Enables personalized treatment effect curves.

**Fang, C. & Liebl, D. (2025).**
Making event study plots honest: a functional data approach to causal
inference.
*arXiv preprint*, 2512.06804.

> Functional Difference-in-Differences estimator converging to a Gaussian
> process. Provides simultaneous confidence bands for event study plots
> with equivalence testing (pre-treatment) and relevance testing
> (post-treatment).

**Ieva, F., Fontana, N., Pivato, C. A., Di Angelantonio, E., & Secchi, P. (2025).**
Enhancing causal inference in functional data: a method for estimating
time-varying causal treatment effects.
In *New Trends in Functional Statistics and Related Fields (IWFOS 2025)*,
Contributions to Statistics. Springer.
DOI: [10.1007/978-3-031-92383-8_35](https://doi.org/10.1007/978-3-031-92383-8_35)

> Evaluates causal effect of binary treatment on functional outcomes with
> sparse and irregularly measured data using a weighting approach combined
> with intervalwise testing.

---

## 2. Functional Treatments (Dose-Response)

**Tan, R., Huang, W., Zhang, Z., & Yin, G. (2025).**
Causal effect of functional treatment.
*Journal of Machine Learning Research*, 26(91), 1--39.
[Link](http://jmlr.org/papers/v26/23-0381.html)

> Three estimators for the average dose-response functional (ADRF) when
> the treatment itself is a curve: functional stabilized weight, outcome
> regression, and doubly robust. Resolves the fundamental problem that
> generalized propensity scores lack a density for functional treatments.
> Applied to EEG data.

**Jiang, Z., Cui, E., & Huling, J. D. (2026).**
Estimating causal effects of functional treatments with modified
functional treatment policies.
*arXiv preprint*, 2602.09145.

> Modified Functional Treatment Policy (MFTP): estimates the average
> potential outcome when each individual slightly modifies their treatment
> trajectory. Uses FPCA decomposition to define population averages over
> functional variables. Applied to NHANES accelerometer data examining
> nighttime disruptions and mortality.

**Zhang, X., Xue, W., & Wang, Q. (2021).**
Covariate balancing functional propensity score for functional treatments
in cross-sectional observational studies.
*Computational Statistics & Data Analysis*, 163, 107303.
DOI: [10.1016/j.csda.2021.107303](https://doi.org/10.1016/j.csda.2021.107303)

> Defines propensity score through top FPC scores (since functional
> densities do not exist). Two covariate balancing methods minimizing
> correlation between FPC scores and covariates. Applied to body shape
> and visceral adipose tissue.

---

## 3. Causal Mediation with Functional Data

**Lindquist, M. A. (2012).**
Functional causal mediation analysis with an application to brain
connectivity.
*Journal of the American Statistical Association*, 107(500), 1297--1309.
DOI: [10.1080/01621459.2012.695640](https://doi.org/10.1080/01621459.2012.695640)

> **Foundational paper.** Extends structural equation models to
> functional mediators via a linear functional SEM (lfSEM) with basis
> function representations and penalized regression. Identifies
> conditions for causal interpretability and introduces an IV approach
> for relaxing sequential ignorability. Applied to fMRI thermal pain
> data, revealing *when* mediation effects occur.

**Coffman, D. L., Dziak, J. J., Litson, K., Chakraborti, Y., Piper, M. E., & Li, R. (2023).**
A causal approach to functional mediation analysis with application to a
smoking cessation intervention.
*Multivariate Behavioral Research*, 58(5), 859--876.
DOI: [10.1080/00273171.2022.2149449](https://doi.org/10.1080/00273171.2022.2149449)

> Combines time-varying effect modeling (TVEM) with functional regression
> within a potential outcomes framework. Extends functional mediation to
> binary outcomes. R package **`funmediation`** on CRAN.

**Zhao, Y., Luo, X., Sobel, M. E., Lindquist, M. A., & Caffo, B. S. (2025).**
Causal functional mediation analysis with an application to functional
magnetic resonance imaging data.
*Biostatistics*, 26(1), kxaf019.
DOI: [10.1093/biostatistics/kxaf019](https://doi.org/10.1093/biostatistics/kxaf019)

> Fully functional mediation: treatment, mediator, AND outcome are all
> continuous functions. Semiparametric functional linear SEMs for
> individual effect curves. R package **`cfma`** on CRAN. Applied to
> task-based fMRI for dynamic brain connectivity.

**Zeng, S., Rosenbaum, S., Alberts, S. C., Archie, E. A., & Li, F. (2021).**
Causal mediation analysis for sparse and irregular longitudinal data.
*Annals of Applied Statistics*, 15(2).
DOI: [10.1214/20-AOAS1427](https://doi.org/10.1214/20-AOAS1427)

> Extends causal mediation to sparse and irregular longitudinal data
> treated as realizations of smooth stochastic processes. Uses FPCA for
> dimension reduction in structural equation models.

**Zhao, Y. & Luo, X. (2019).**
Granger mediation analysis of multiple time series with an application
to fMRI.
*Biometrics*, 75(3), 788--798.
DOI: [10.1111/biom.13056](https://doi.org/10.1111/biom.13056)

> Integrates causal mediation analysis with vector autoregressive (VAR)
> models for time series. Extends to multilevel data. Applied to fMRI
> brain pathway analysis.

---

## 4. Causal Discovery (Functional DAGs and SEMs)

**Lee, K.-Y. & Li, L. (2022).**
Functional structural equation model.
*Journal of the Royal Statistical Society Series B*, 84(2), 600--629.
DOI: [10.1111/rssb.12471](https://doi.org/10.1111/rssb.12471)

> Functional SEM for estimating directional (causal) relations from
> multivariate functional data. Two-layer RKHS for nonlinear
> function-on-function relationships. Score function at the linear
> operator level recovers the true directional order without requiring
> the faithfulness condition. Applied to brain connectivity.

**Lee, K.-Y., Li, L., & Li, B. (2024).**
Functional directed acyclic graphs.
*Journal of Machine Learning Research*, 25(78), 1--48.
[Link](https://www.jmlr.org/papers/v25/22-1038.html)

> Estimates DAGs from multivariate functional data using conditional
> covariance and partial correlation operators for conditional
> independence testing. Adapts the PC-algorithm to functional directed
> graphs. Proves uniform convergence rates and consistency even as graph
> size diverges. Applied to proteomic time-course data.

**Zhou, F., He, K., Wang, K., Xu, Y., & Ni, Y. (2023).**
Functional Bayesian networks for discovering causality from multivariate
functional data.
*Biometrics*, 79(4), 3279--3293.
DOI: [10.1111/biom.13922](https://doi.org/10.1111/biom.13922)

> FLiNG-BN: encodes causal structures with a DAG using basis coefficient
> modeling. Non-Gaussianity enables unique causal identification even
> with noisy measurements. Fully Bayesian inference (MCMC) for
> simultaneous basis function and DAG estimation. Applied to EEG brain
> connectivity.

**Roy, S., Wong, R. K. W., & Ni, Y. (2023).**
Directed cyclic graph for causal discovery from multivariate functional
data.
*Advances in Neural Information Processing Systems (NeurIPS)*, 36,
42762--42774.

> FENCE: first causal model for functional data with cycles (feedback
> loops). Operator-based non-recursive linear SEM with a
> low-dimensional causal embedded space preserving all relevant causal
> information. Fully Bayesian with uncertainty quantification. Applied
> to EEG brain connectivity in alcoholism studies.

**Lan, T., Li, Z., Lin, J., Li, Z., Bai, L., Li, M., Tsung, F., Zhao, R., & Zhang, C. (2024).**
MultiFun-DAG: multivariate functional directed acyclic graph.
*arXiv preprint*, 2404.13836.

> First DAG model for nodes representing multivariate functional data.
> Hidden bilinear multivariate function-to-function regression with
> low-rank decomposition. EM algorithm with acyclicity constraints.
> Applied to urban traffic congestion.

**Li, T., Fan, E., Li, T., & Zhu, H. (2026).**
Causal inference in biomedical imaging via functional linear structural
equation models.
*arXiv preprint*, 2601.20610.

> FLSEM for uncovering causal links among genetic, imaging, and clinical
> data under unobserved confounding. L0-penalized three-step estimator.
> Uses scalar instrumental variables for identification. Applied to UK
> Biobank imaging data.

**Lee, K.-Y. & Li, L. (2024).**
Functional SEM with latent variables.
*arXiv preprint*, 2412.19242.

> Extends FSEMs to incorporate latent variables modeled as Gaussian
> Processes. Supports concurrent, historical, linear, and smooth effects
> with penalized likelihood estimation.

**Kennerberg, P. & Wit, E. C. (2025).**
Functional SEMs with out-of-sample guarantees.
*arXiv preprint*, 2503.20072.

> Risk minimization framework for functional SEMs using linear,
> potentially unbounded operators. Novel worst-risk decomposition
> theorem for distributional shifts.

---

## 5. Granger Causality for Functional Time Series

**Sen, R., Majumdar, A., & Sikaria, S. (2022).**
Bayesian testing of Granger causality in functional time series.
*Journal of Quantitative Economics*, 20(1).
DOI: [10.1007/s40953-022-00306-x](https://doi.org/10.1007/s40953-022-00306-x)

> Multivariate functional autoregressive model (MFAR) with Bayes Factor
> for Granger causality detection. Applied to yield curve causality
> between countries and meteorological factors causing pollution.

**Shang, H. L., Ji, K., & Beyaztas, U. (2021).**
Granger causality of bivariate stationary curve time series.
*Journal of Forecasting*, 40(4), 626--635.
DOI: [10.1002/for.2732](https://doi.org/10.1002/for.2732)

> Tests Granger causality using generalized measures of correlation for
> bivariate curve time series. Applied to sea surface temperature and
> atmospheric pressure.

**Beyaztas, U. & Shang, H. L. (2021).**
Dynamic functional principal components for testing causality.
*Foundations*, 2(2), 22.

> Projects functional time series onto dynamic FPCs, then applies
> standard multivariate Granger tests to the resulting score vectors.

**Saumard, M. (2017).**
Linear causality in the sense of Granger with stationary functional
time series.
In *Functional Statistics and Related Fields*, pp. 225--231. Springer.
DOI: [10.1007/978-3-319-55846-2_30](https://doi.org/10.1007/978-3-319-55846-2_30)

> Formal definition of Granger causality for functional time series via
> ordering of covariance operators.

**Hormann, S., Kidzinski, L., & Hallin, M. (2015).**
Dynamic functional principal components.
*Journal of the Royal Statistical Society Series B*, 77(2), 319--348.
DOI: [10.1111/rssb.12076](https://doi.org/10.1111/rssb.12076)

> Foundational work on frequency-domain dynamic FPCA for functional time
> series. Provides the dimension reduction framework that subsequent
> Granger causality tests build upon.

---

## 6. Instrumental Variables for Functional Data

**Florens, J.-P. & Van Bellegem, S. (2015).**
Instrumental variable estimation in functional linear models.
*Journal of Econometrics*, 186(2), 465--476.
DOI: [10.1016/j.jeconom.2015.02.020](https://doi.org/10.1016/j.jeconom.2015.02.020)

> Foundational IV paper for functional linear models. Both covariates
> and instruments are functional. Shows that estimation constitutes an
> ill-posed inverse problem; uses Tikhonov regularization. Introduces
> "instrument strength" in the functional setting.

**Seong, D. & Seo, W.-K. (2025).**
Functional instrumental variable regression with an application to
estimating the impact of immigration on native wages.
*Econometric Theory*, 41, 1248--1283.
DOI: [10.1017/S0266466624000252](https://doi.org/10.1017/S0266466624000252)

> FPCA-based IV estimators with full distributional theory for practical
> inference. Applied to immigration and native worker wages.

**Petrovich, J., Taoufik, B., & Davis, Z. G. (2023).**
Instrumental variable estimation for functional concurrent regression
models.
*Journal of Applied Statistics*, 51(8), 1570--1589.
DOI: [10.1080/02664763.2023.2229968](https://doi.org/10.1080/02664763.2023.2229968)

> Two-stage least squares for functional concurrent regression with
> endogenous predictors. Adapted for sparse functional data using PACE.
> Applied to U.S. labor supply elasticities.

**Benatia, D., Carrasco, M., & Florens, J.-P. (2017).**
Functional linear regression with functional response.
*Journal of Econometrics*, 201(2), 269--291.
DOI: [10.1016/j.jeconom.2017.08.008](https://doi.org/10.1016/j.jeconom.2017.08.008)

> IV estimation where both regressor and response are functions.
> Tikhonov regularization. Applied to immigration and native wage
> dynamics.

---

## 7. Counterfactual and Quasi-Experimental Methods

**Kurisu, D., Zhou, Y., Otsu, T., & Muller, H.-G. (2025).**
Regression discontinuity designs for functional data and random objects
in geodesic spaces.
*arXiv preprint*, 2506.18136.

> Generalizes RDD to functional outcomes and non-Euclidean data. Causal
> effects defined as geodesics between local Frechet means. Applied to
> CO concentration curves near Taipei Metro.

**Kurisu, D., Zhou, Y., Otsu, T., & Muller, H.-G. (2025).**
Geodesic synthetic control methods for random objects and functional
data.
*arXiv preprint*, 2505.00331.

> Extends synthetic control to geodesic metric spaces. Includes geodesic
> synthetic difference-in-differences with double robustness. Applied to
> Japan earthquake employment effects, German reunification fertility
> patterns, and Soviet collapse age-at-death distributions.

**Van Dijcke, D. (2025).**
Regression discontinuity design with distribution-valued outcomes.
*arXiv preprint*, 2504.03992.

> RDD for distribution-valued/functional outcomes. Introduces the local
> average quantile treatment effect (LAQTE). R package **`R3D`**.
> Applied to gubernatorial party effects on income distributions.

**Carrizosa, E., Ramirez-Ayerbe, J., & Romero Morales, D. (2024).**
A new model for counterfactual analysis for functional data.
*Advances in Data Analysis and Classification*, 18, 981--1000.
DOI: [10.1007/s11634-023-00563-5](https://doi.org/10.1007/s11634-023-00563-5)

> Counterfactual explanation methodology for functional data classifiers.
> Identifies which samples can be combined to create minimal-distance
> counterfactuals.

**Fontana, N., Ieva, F., Zuccolo, L., Di Angelantonio, E., & Secchi, P. (2025).**
Unraveling time-varying causal effects of multiple exposures: integrating
FDA with multivariable Mendelian randomization.
*arXiv preprint*, 2512.19064.

> Multivariable Functional Mendelian Randomization (MV-FMR) combining
> FPCA with cross-validation for basis selection. Handles overlapping
> instruments and mediation. Applied to UK Biobank blood pressure and
> BMI effects on coronary artery disease.

---

## 8. Nonparametric and Kernel Methods

**Kurisu, D., Otsu, T., & Xu, M. (2026).**
Nonparametric causal inference with functional covariates.
*Journal of Business & Economic Statistics*, 44(1), 53--66.
DOI: [10.1080/07350015.2025.2501563](https://doi.org/10.1080/07350015.2025.2501563)

> IPW estimator for ATE when covariates include functional variables.
> Propensity score estimated via kernel methods for functional data.
> Establishes sqrt(n)-consistency and asymptotic normality.

**Raykov, Y. P., Luo, H., Strait, J. D., & KhudaBukhsh, W. R. (2025).**
Kernel-based estimators for functional causal effects.
*arXiv preprint*, 2503.05024.

> Causal effect estimators using empirical Frechet means and
> operator-valued kernels. Uses Fisher-Rao metric for
> outcome/covariate registration. Avoids explicit propensity function
> estimation.

**Lundborg, A. R., Shah, R. D., & Peters, J. (2022).**
Conditional independence testing in Hilbert spaces with applications to
functional data analysis.
*Journal of the Royal Statistical Society Series B*, 84(5), 1821--1850.
DOI: [10.1111/rssb.12544](https://doi.org/10.1111/rssb.12544)

> Generalised Hilbertian Covariance Measure (GHCM) for testing
> conditional independence when X, Y, Z may all be functional. Critical
> building block for causal discovery with functional data. R package
> **`ghcm`** on CRAN.

---

## 9. Personalized Treatment Regimes with Functional Data

**McKeague, I. W. & Qian, M. (2014).**
Estimation of treatment policies based on functional predictors.
*Statistica Sinica*, 24(3), 1461--1485.
DOI: [10.5705/ss.2012.196](https://doi.org/10.5705/ss.2012.196)

> Model-based estimation of optimal treatment policies from RCTs using
> functional regression. Establishes asymptotic optimality and develops
> prediction intervals for policy value.

**Laber, E. B. & Staicu, A.-M. (2018).**
Functional feature construction for individualized treatment regimes.
*Journal of the American Statistical Association*, 113(523), 1219--1227.
DOI: [10.1080/01621459.2017.1321545](https://doi.org/10.1080/01621459.2017.1321545)

> Functional Q-learning for optimal individualized treatment regime
> estimation from sparse, irregularly-spaced longitudinal data. Uses
> FPCA for data-driven feature construction. Applied to the STAR*D
> depression trial.

**Ciarleglio, A., Petkova, E., Ogden, T., & Tarpey, T. (2018).**
Constructing treatment decision rules based on scalar and functional
predictors when moderators of treatment effect are unknown.
*Journal of the Royal Statistical Society Series C*, 67(5), 1331--1356.
DOI: [10.1111/rssc.12278](https://doi.org/10.1111/rssc.12278)

> Modified covariates method with augmentation (MC-A) for treatment
> decisions using both scalar and functional predictors. Group lasso for
> variable selection. Applied to EEG-based depression treatment.

**Kong, K., Guan, L., & Zhang, Z. (2025).**
Flexible regression methods for estimating optimal individualized
treatment regimes with scalar and functional covariates.
*Statistical Methods in Medical Research*, 34(7), 1459--1479.
DOI: [10.1177/09622802251340259](https://doi.org/10.1177/09622802251340259)

> Regression-based framework integrating scalar and functional
> covariates for optimal treatment regime estimation in both RCTs and
> observational studies.

---

## 10. Neuroscience Applications

**Cao, X., Sandstede, B., & Luo, X. (2019).**
A functional data method for causal dynamic network modeling of
task-related fMRI.
*Frontiers in Neuroscience*, 13, 127.
DOI: [10.3389/fnins.2019.00127](https://doi.org/10.3389/fnins.2019.00127)

> Causal dynamic network (CDN) method that estimates brain activations
> and connections simultaneously. Links observed fMRI with latent
> neuronal states modeled by ODEs using FDA basis function expansion.

**Gao, X., Wang, J., Hu, G., & Sun, J. (2024).**
Functional causal inference with time-to-event data.
*Statistics in Biosciences*.
DOI: [10.1007/s12561-024-09439-4](https://doi.org/10.1007/s12561-024-09439-4)

> Functional Accelerated Failure Time (FAFT) model with regression
> adjustment, functional IPW, and doubly robust estimators for survival
> outcomes. Applied to ADNI Alzheimer's data.

---

## Software

| Package | Language | Source | Purpose |
|---------|----------|--------|---------|
| `cfma` | R | CRAN | Causal functional mediation (fully functional) |
| `funmediation` | R | CRAN | Functional mediation with binary/continuous outcomes |
| `ghcm` | R | CRAN | Conditional independence testing in Hilbert spaces |
| `R3D` | R | GitHub | RDD for distribution-valued/functional outcomes |

---

## Key Research Groups

- **Reimherr & Kennedy** (Penn State / CMU) -- doubly robust FATE
- **Lindquist & Caffo** (Johns Hopkins) -- functional mediation, fMRI
- **Lee & Li** (UC Berkeley) -- functional SEMs, causal discovery
- **Ni & Wong** (Texas A&M) -- Bayesian functional DAGs, cyclic graphs
- **Kurisu & Otsu** (U. Tokyo / LSE) -- nonparametric causal FDA, RDD,
  synthetic control
- **Zhu** (UNC Chapel Hill) -- FLSEM, biomedical imaging
- **Fan Li** (Duke) -- causal mediation for sparse/functional data

---

## Open Directions

- **TMLE / targeted learning** for functional outcomes -- not yet developed
- **Pearl's do-calculus** has not been formally connected to FDA
- **Semiparametric efficiency bounds** for functional outcomes remain open
- **DiD with fully functional outcomes** is very new (Fang & Liebl 2025)
- No unified package covering all causal methods for functional data
