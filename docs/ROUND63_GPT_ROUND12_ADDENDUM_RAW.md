# R12 full-text addendum — definitive equation-by-equation certification

**Source reviewed:** the complete 12-page publisher PDF of F. Grönberg, M. Danielsson, and M. Sjölin, “Count statistics of nonparalyzable photon-counting detectors with nonzero pulse length,” *Medical Physics* **45**(8), 3800–3811 (2018), DOI 10.1002/mp.13063.

> **VERDICT: X.5 CLAIM CONFIRMED AFTER MANDATORY AMENDMENT (not verbatim).**
>
> **Confidence: FULL-TEXT-CERTIFIED.**

The core ROUND63 residual result survives: the paper contains neither the active-start detector-channel information surface `I_loglambda(rho,nu)`, nor a two-dimensional principal ridge, nor the cube-root law, nor the terminal-phase missing-information identity. However, the earlier R12 wording understated Grönberg et al. in one material respect. They **do** compute exact finite-window Fisher information from a scalar integrated-count PMF for a single-threshold CT density task, and their fixed-frame Fisher/SNR curve has an interior maximum followed by high-rate decline. The manuscript must acknowledge that result explicitly.

## Executive rulings

| R12 question | Full-text ruling |
|---|---|
| **1. Any finite-window Fisher-information analysis of the scalar integrated count?** | **Yes.** Equations (3), (28)–(29), and (38)–(40), together with Sec. 2.D, compute task-specific Fisher information for soft-tissue thickness from the finite-frame scalar count `Y(t)`. The exact seminonparalyzable PMF is evaluated through Eq. (3); Poisson and Gaussian alternatives use Eqs. (5) and (4). What is absent is Fisher information for `log lambda` treated as a detector-channel quantity and mapped over independent axes `(rho,nu)`. |
| **2. Any information optimum/ridge, `nu^(1/3)` scaling, or terminal-phase decomposition?** | **A one-dimensional task-specific optimum is present.** Figure 6(a) and the Discussion show that fixed-frame Fisher information/SNR rises and then decreases above a count rate. **No** two-dimensional `(rho,nu)` ridge, analytic maximizer, `nu^(1/3)` scaling, terminal residual-phase variable, or identity `I_N=E[N]-rho^2 E[Var(R_nu|N)]` appears. |
| **3. Anything narrowing/endangering X.5?** | **Yes, wording only.** X.5 must not imply that Grönberg et al. considered only moments/Gaussian surrogates, only multibin spectral-CT metrics, or no finite-window scalar-count optimum. The narrow ROUND63 formulas and asymptotic ridge remain unanticipated in the full text. |

---

# 1. Equation-by-equation certification

The categories below are: **PMF/count law**, **moments/renewal asymptotics**, **Fisher/task metric**, **Gaussian/Poisson surrogate**, **multibin extension**, and **other/simulation**.

| Eq. | What it establishes | Classification | Relevance to ROUND63 |
|---:|---|---|---|
| **(1)** | Standard ideal delta-pulse nonparalyzable formulas `E[Y(t)]=lambda t/(1+lambda tau)` and `Var[Y(t)]=lambda t/(1+lambda tau)^3`. | **Moments** | This is the same large-count moment structure underlying renewal-CLT surrogates; it is not a finite-window PMF or Fisher ridge. |
| **(2)** | Generic Fisher information definition for parameter `theta` from observation `Y(t)`. | **Fisher/task metric** | Establishes the general CRLB framework only. |
| **(3)** | Exact discrete Fisher information `sum_k (partial_theta p_Y[k;theta])^2/p_Y[k;theta]` for a finite-support count PMF. | **Fisher/task metric — exact** | This is the equation later combined with the exact finite-frame PMF (29). It is the decisive evidence that the paper performs exact scalar-count FI, albeit for material thickness rather than `log lambda`. |
| **(4)** | Fisher information of a multivariate Gaussian with parameter-dependent mean and covariance. | **Gaussian surrogate metric** | Used with asymptotic moments in the Gaussian approximation. No finite-window ridge derivation. |
| **(5)** | Fisher information of independent Poisson variables. | **Poisson surrogate metric** | Used for the Poisson approximation. |
| **(6)** | Mean vector and covariance matrix of energy-bin counts from total-count mean/variance and conditional bin probabilities. | **Multibin extension / moments** | The authors explicitly note exactness only for the delta-pulse model and possible missing correlations for nonzero pulses. It is not used for the actual Sec. 2.D task, because no spectral model is developed there. |
| **(7)** | Exponential waiting-time density. | **Other / PMF precursor** | Basic renewal ingredient. |
| **(8)** | Corresponding exponential CDF. | **Other / PMF precursor** | Basic renewal ingredient. |
| **(9)** | Heaviside function definition. | **Other** | Not an information result. |
| **(10)** | Event decomposition for the probability that the detector free time `T` is at most `Delta t`, including pre-window events in the seminonparalyzable part. | **PMF/count-law precursor** | Encodes their start/prehistory convention; it is not the ROUND63 active-start convention. |
| **(11)** | Defines `p`, `q`, and `F(Delta t)` as Poisson no-event/event probabilities. | **PMF/count-law precursor** | Leads to the mixed free-time law. |
| **(12)** | Free-time CDF `P(T<=Delta t)=pU+qF`. | **PMF/count-law precursor** | Finite-time renewal ingredient. |
| **(13)** | Mixed free-time density with a point mass at zero plus an exponential tail. | **PMF/count-law precursor** | The atom at zero reflects arrivals before `t=0`; this is a physically richer, different initial condition from active-start. |
| **(14)** | Renewal-theory asymptotic Gaussian mean and variance of the count process in terms of renewal-time moments. | **Moments / Gaussian asymptotics** | Explicitly asymptotic; no exact finite-window FI. |
| **(15)** | Defines time-normalized asymptotic mean and variance. | **Moments** | Used to fit simulation and construct approximations. |
| **(16)** | Mean and variance of the mixed free time `T`. | **Moments** | Derived in Appendix A, Proposition 1. |
| **(17)** | Closed-form asymptotic time-normalized mean and variance for the seminonparalyzable process. | **Moments / Gaussian surrogate input** | Feeds the Gaussian and Poisson performance approximations. |
| **(18)** | Defines the density of the sum of `n` free times as an `n`-fold convolution. | **PMF derivation** | Beginning of the exact finite-window count-law derivation. |
| **(19)** | Expands the convolution of the point-mass/exponential mixture into a binomial mixture of Erlang densities. | **PMF derivation** | Exact count-law ingredient. |
| **(20)** | Defines the binomial coefficient/weight `B_nk`. | **PMF derivation** | Mixture weight. |
| **(21)** | Erlang density for the sum of `k` exponential variables. | **PMF derivation** | Re-derived as Eq. (A1) in Appendix A. |
| **(22)** | Erlang CDF. | **PMF derivation** | Proved in Appendix A, Proposition 3. |
| **(23)** | Defines Poisson probability notation `P_m(x)`. | **Other / PMF notation** | Used to express the finite-time CDF compactly. |
| **(24)** | Exact CDF of the sum `S_n` of free times. | **PMF derivation** | Finite-window count-law ingredient. |
| **(25)** | Binomial normalization identity showing the mixture weights sum to one. | **PMF derivation / normalization** | No Fisher result. |
| **(26)** | Relates the event `Y(t)>=n` to `S_n+n tau<=t`. | **Finite-window count law** | Converts renewal times into count probabilities. |
| **(27)** | Expresses the finite-window count CDF through the CDF of `S_n`. | **Finite-window count law** | Exact finite-frame relation. |
| **(28)** | Gives the explicit piecewise finite-window CDF of `Y(t)` for `0<=n<=floor(t/tau)`. | **Finite-window count CDF** | This is an exact finite-support scalar-count law. The paper notes dependence on dimensionless `lambda tau`, `tau_s/tau`, and `t/tau`. |
| **(29)** | Obtains the exact finite-window scalar-count PMF by differencing adjacent CDF values. | **Finite-window count PMF — exact** | Combined with Eq. (3) in Sec. 2.D. This directly invalidates the earlier characterization that no exact finite-window scalar-count FI was present. |
| **(30)** | Discrete pulse-train simulation with photon energies and additive Gaussian noise. | **Other / Monte Carlo electronics** | Validation machinery only. |
| **(31)** | Convolution of the pulse train with the pulse-shaping kernel. | **Other / simulation** | Validation machinery. |
| **(32)** | Chosen unipolar pulse shape. | **Other / detector simulation** | Physical simulation model, not count FI. |
| **(33)** | Discretization and truncation of the pulse kernel. | **Other / simulation** | No information result. |
| **(34)** | Baseline-holder forward-difference update. | **Other / electronics simulation** | No information result. |
| **(35)** | Piecewise-clipped baseline-holder dynamics. | **Other / electronics simulation** | No information result. |
| **(36)** | Baseline subtraction and downsampling to the digital counter input. | **Other / simulation** | No information result. |
| **(37)** | Standard errors of the Monte Carlo sample mean and sample variance. | **Other / validation statistics** | Used only for stopping simulation. |
| **(38)** | Transmitted input count rate as a function of soft-tissue thickness `theta`. | **Task model** | Supplies the chain rule from the detector count law to the material-thickness Fisher information. |
| **(39)** | Defines fixed-frame `SNR^2(lambda)` as Fisher information and fixed-dose dose efficiency `F/(lambda t)`. | **Fisher/task metric** | Creates two one-dimensional paths through finite-window parameter space: fixed `t`, and fixed `lambda t`. It is not an independent `(rho,nu)` map. |
| **(40)** | Uses the CRLB to show maximal task SNR squared is proportional to Fisher information for fixed `theta`. | **Fisher/task interpretation** | Justifies plotting the FI as `SNR^2`. |
| **(A1)** | Repeats/derives the Erlang convolution density used in Eq. (21). | **PMF derivation** | No Fisher calculation. |

## Appendix derivations

- **Appendix A, Proposition 1:** derives Eq. (16) by computing `E[T]`, `E[T^2]`, and `Var[T]`. Classification: **moments**.
- **Appendix A, Proposition 2:** proves the Erlang convolution density by induction, yielding Eq. (A1)/(21). Classification: **PMF derivation**.
- **Appendix A, Proposition 3:** differentiates the proposed Erlang CDF and invokes the fundamental theorem of calculus to prove Eq. (22). Classification: **PMF derivation**.
- **Algorithm 1:** simulates the nonparalyzable digital counter. Classification: **Monte Carlo validation**, not an analytic derivation.

There is **no appendix derivation of Fisher information, no differentiation of the FI with respect to count rate, no asymptotic optimization over frame length, and no missing-information calculation**.

---

# 2. Definitive ruling on R12 Question 1

## Yes: exact finite-window scalar-count Fisher information is present

Section 2.D states that the authors compute `F(theta;lambda)` as a function of input count rate using:

1. the **exact seminonparalyzable PMF (29)** with `tau_s/tau=0.9`, scored by the exact discrete FI formula **(3)**;
2. the ideal nonparalyzable special case of the same PMF;
3. a Poisson approximation scored by **(5)**;
4. a Gaussian approximation using asymptotic mean and variance, scored by **(4)**.

The observation is one scalar single-threshold integrated count `Y(t)`. The exact PMF has finite support `0,...,floor(t/tau)`, so this is genuinely a **finite-window exact scalar-count FI calculation** within their seminonparalyzable model.

The distinction from ROUND63 is narrower:

- Grönberg parameterizes the task by soft-tissue thickness `theta`, with `lambda=lambda(theta)` from Eq. (38); by chain rule their task FI contains the detector count information multiplied by task sensitivity.
- ROUND63 isolates Fisher information for **`log lambda` itself**, independent of a material task.
- Grönberg evaluates count rate along one fixed-frame slice and one fixed-dose path; ROUND63 treats `rho=lambda tau` and `nu=T/tau` as independent axes and derives their asymptotics.
- Grönberg uses a seminonparalyzable nonzero-pulse-length model with a prehistory/start atom; ROUND63’s theorem is for the ideal active-start nonparalyzable model.

Accordingly, the sentence “Grönberg did not analyze finite-window scalar-count Fisher information” is **superseded and must not appear**.

---

# 3. Definitive ruling on R12 Question 2

## Yes: a task-specific one-dimensional information maximum is present

Figure 6(a) plots fixed-frame `SNR^2`, which by Eqs. (39)–(40) is proportional to Fisher information. The exact seminonparalyzable curve rises, reaches an interior maximum near normalized input rate `lambda tau` of order two, and then decreases. The Discussion states explicitly that the proposed model “predicts a decrease in SNR2 above a certain count rate,” whereas the ideal nonparalyzable delta-pulse model predicts a monotone increase.

Therefore the manuscript must not imply that Grönberg et al. did not observe an information optimum or high-rate information decline.

## No: the ROUND63 ridge structure is absent

The full text contains none of the following:

- an independent two-axis map in `(rho,nu)`;
- optimization of detector-channel `I_N(log lambda)` over `rho` for each `nu`;
- an analytic stationarity equation for the maximizing load;
- a ridge law of order `nu^(1/3)`;
- the expansion `rho*(nu)=(6nu)^(1/3)-2/3+O(nu^(-1/3))`;
- a terminal residual dead-time phase `R_nu`;
- a complete-data/missing-data Fisher identity;
- the formula `I_N=E[N]-rho^2 E[Var(R_nu|N)]`.

Their fixed-frame curve is one task-specific slice. Their fixed-dose curve follows `lambda t=constant`, hence a one-dimensional inverse relation between load and window length. Neither is analyzed as a principal ridge or as a large-`nu` asymptotic problem.

---

# 4. What the full text changes in the frozen claim

## Mandatory corrections

The prior R12 issue body contained two statements that are now superseded:

1. **Superseded:** “No accessible material shows Grönberg evaluating exact finite-window scalar-count FI.”  
   **Correct:** Eqs. (3), (28)–(29), and Sec. 2.D do exactly that for a material-thickness task.

2. **Superseded:** “Grönberg’s Fisher use is spectral-CT image-performance FI.”  
   **Correct:** the actual numerical Fisher example is a **single-threshold soft-tissue density/thickness task**; Eq. (6) merely discusses a possible multibin extension, and the authors explicitly say no spectral model was developed in the paper.

## Residual claim that survives

No numbered equation, appendix derivation, figure, or discussion passage anticipates the two ROUND63 formulas. The defensible residual contribution is not the existence of finite-window scalar-count FI or a finite-rate optimum. It is the combination:

1. detector-channel FI for `log lambda` rather than a material-specific task;
2. independent dimensionless axes `(rho,nu)` for the ideal active-start count-only channel;
3. a terminal-phase missing-information decomposition;
4. the analytic principal-ridge asymptotic `rho*(nu)=(6nu)^(1/3)-2/3+O(nu^(-1/3))`.

---

# 5. Amended frozen manuscript subsection X.5

The previous X.5 text does **not** survive verbatim. Replace it with the following.

```latex
\subsection{Relation to prior work}

Müller systematized renewal-theoretic analysis of extended and nonextended
dead-time counters, relating distorted interarrival laws to count statistics
and measurement-origin effects~\cite{Muller1973}. Yu and Fessler derived exact
first and second moments for several dead-time counting models and used those
moments to study intensity correction and the adequacy of Poisson
reconstruction models~\cite{YuFessler2000}. Wang et al. derived pileup count
statistics for idealized energy-discriminating nonparalyzable x-ray detectors
and evaluated task-specific Cramér--Rao performance as a function of input
rate~\cite{Wang2011}. Alvarez used renewal-process large-count moments and
multivariate-Gaussian count models to quantify pileup-induced covariance,
signal-to-noise loss, and material-estimation Cramér--Rao bounds
~\cite{Alvarez2014}. Grönberg et al. derived an exact finite-frame count
probability mass function for a seminonparalyzable detector with nonzero pulse
length, validated its moments against an electronic Monte Carlo model, and
used the exact distribution together with Poisson and Gaussian approximations
to evaluate Fisher information for a single-threshold soft-tissue-density task
as a function of count rate; their fixed-frame Fisher curve exhibits an
interior maximum followed by high-rate decline~\cite{Gronberg2018}.

A complementary line retains richer event information. Rapp et al. modeled
the sequence of dead-time-distorted detection times as a Markov chain,
enabling accurate high-flux single-photon lidar and large acquisition-time
reductions~\cite{Rapp2021}. Jorgensen and Johnson derived asymptotic
Fisher-information rates and efficient estimators for periodic dead-time
event-detection processes, showing that activation and event statistics can
retain information discarded by conventional histogram summaries
~\cite{JorgensenJohnson2026}.

Here we isolate a different detector-channel question: an ideal
nonparalyzable detector is active at the start of a finite exposure, only its
integrated scalar count is retained, and the parameter of interest is the log
arrival rate itself. Treating $\rho=\lambda\tau$ and $\nu=T/\tau$ as
independent axes, we derive the exact finite-window information surface, the
terminal-phase missing-information identity
\begin{equation}
 I_N=\mathbb{E}[N]-\rho^2\mathbb{E}\!\left[
 \operatorname{Var}(R_\nu\mid N)\right],
\end{equation}
and the principal information ridge
\begin{equation}
 \rho^*(\nu)=(6\nu)^{1/3}-\frac{2}{3}+O(\nu^{-1/3}).
\end{equation}
These statements concern the active-start ideal count-only channel and do not
replace task-specific finite-rate optima for nonzero-pulse-length detectors or
information bounds for time-resolved dead-time event records.
```

## Optional nuclear-counting sentence

The optional Bécares–Blázquez sentence from the original R12 ruling remains approved and is useful to prevent an overbroad novelty reading.

---

# 6. Final claim discipline

## Permitted

> We derive the detector-channel finite-window information surface for log arrival rate under the active-start ideal nonparalyzable scalar-count model.

> We derive a terminal-phase missing-information decomposition and the associated cube-root principal ridge.

> Prior finite-window task Fisher analyses already showed count-rate-dependent maxima and decline; our result supplies a two-axis detector-channel asymptotic law and decomposition.

## Not permitted

- “Grönberg considered only mean and variance.”
- “Grönberg used only a Gaussian approximation.”
- “Grönberg did not compute finite-window Fisher information.”
- “Prior work did not find an information optimum under dead time.”
- “The existence of a high-rate FI decline is new.”
- “Grönberg’s numerical FI study was multibin spectral CT.”
- categorical “first” or “no prior work” language.

## Final provenance status

> **FULL-TEXT CERTIFICATION COMPLETE.** The Grönberg paper narrows the novelty wording but does not contain or anticipate the ROUND63 terminal-phase identity or cube-root `(rho,nu)` ridge. The amended X.5 above is the frozen submission wording.