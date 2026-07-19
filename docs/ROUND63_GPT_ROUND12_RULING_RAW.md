# R12 — Final ruling: Grönberg check and X.5 frozen wording

**Scope:** this ruling closes the final theory/provenance blocker for the current Optics Express manuscript. It does not alter the frozen campaign results, the R9 manuscript scope lock, or the R10–R11 follow-up method line.

## Executive verdict

| Question | Ruling | Confidence |
|---|---|---|
| Does Grönberg–Danielsson–Sjölin (2018) contain a finite-window Fisher-information analysis of the scalar integrated count? | **No such analysis was found in the accessible record.** The paper derives a nonzero-pulse-length count-statistics model and uses Fisher information for spectral-CT image-performance metrics; the accessible evidence does not show an exact scalar-count Fisher calculation indexed by finite window `T` or by `(rho,nu)`. | **Moderate-high**, not equation-by-equation certified because the full article PDF was inaccessible in this environment. |
| Does it derive an information optimum/ridge, `nu^(1/3)` law, or terminal-phase missing-information identity? | **No such result was found.** Nothing accessible gives an optimization of scalar count information over load, a `nu^(1/3)` scale, or `I_N = E[N] - rho^2 E[Var(R_nu|N)]`. | **Moderate-high**, subject to the same access limitation. |
| Does another located work anticipate the exact cube-root ridge or identity? | **No located prior gives either formula.** Broad count-rate optima, pileup-induced CRLB degradation, renewal moments, and high-flux recovery are well established, but the active-start finite-window scalar-count ridge and terminal-phase decomposition were not found. | **Moderate-high best-effort literature verdict; absence is not mathematically provable from search alone.** |
| Manuscript wording | The X.5 text below is **frozen and approved**. It deliberately avoids priority adjectives and restricts the residual claim to the active-start, nonparalyzable, finite-window, scalar integrated-count observation. | High. |

---

# 1. Full-text access and bibliographic correction

## 1.1 Access limitation

I could not obtain the complete typeset article or thesis-bundle copy in this tool environment:

- Wiley’s article/full-text endpoints returned paywall/HTTP-access failures.
- The KTH DiVA thesis bundle advertises a full PDF containing Paper I, but its `FULLTEXT01.pdf` endpoint repeatedly returned a cache miss here.
- The ResearchGate entry is request-only.

Accordingly, this is **not an honest equation-by-equation certificate of the publisher PDF**. The ruling is based on the strongest accessible evidence:

1. the official Wiley article record, abstract, bibliographic data, and full reference list;
2. the KTH thesis record identifying the article as Paper I and describing the thesis contribution;
3. the companion SPIE paper’s accessible scope statement;
4. later papers from the same detector-modeling line;
5. the cited predecessor literature and a targeted search for the proposed ridge and missing-information formulas.

If an author manuscript or institutional PDF becomes available before submission, a human should still search the full article for `Fisher`, `information`, `Cramér`, `measurement time`, `count rate`, `maximum`, `optimal`, `renewal`, and `asymptotic`. No manuscript claim below depends on asserting that this search was already completed.

## 1.2 Bibliographic correction

The correct citation is:

> F. Grönberg, M. Danielsson, and M. Sjölin, “Count statistics of nonparalyzable photon-counting detectors with nonzero pulse length,” *Medical Physics* **45**(8), 3800–3811 (2018), DOI: [10.1002/mp.13063](https://doi.org/10.1002/mp.13063).

It is **Volume 45, Issue 8**, not Issue 10.

---

# 2. What Grönberg et al. established, and where their analysis ends

The official article record states that Grönberg et al.:

1. derive an analytical statistical count-distribution model for nonparalyzable photon-counting detectors with nonzero pulse length;
2. validate that model against an event-level Monte Carlo detector simulation;
3. show that a simple Gaussian approximation accurately reproduces the modeled behavior and performance;
4. compute **spectral-CT image-performance metrics using Fisher information**, comparing the proposed count model, its approximations, and the ideal nonparalyzable model.

The reference structure is also diagnostic: the paper cites spectral-CT basis-image CRLB work, energy-bin optimization, Wang’s pileup statistics, Yu–Fessler moments, Müller’s renewal/dead-time review, and Gaussian/count-covariance performance analyses. The accessible companion paper explicitly describes renewal-theoretic **asymptotic mean and variance** calculations followed by pileup correction.

## Ruling on Question 1

No accessible material shows Grönberg et al. evaluating

\[
I_{N_T}(\log\lambda)
=\sum_m p_m(\lambda,T)
\left(\partial_{\log\lambda}\log p_m(\lambda,T)\right)^2
\]

for the single scalar integrated count as an exact finite-window object, nor mapping it over

\[
\rho=\lambda\tau,
\qquad
\nu=T/\tau.
\]

Their Fisher-information use is best characterized as **task/image-performance Fisher information built from a detector count-statistics model**, including a Gaussian approximation, rather than the exact finite-window scalar-count information studied in ROUND63.

This distinction must be made narrowly: Grönberg’s detector model is physically richer than the ideal zero-pulse-length model used for the ridge theorem. The residual claim is not “our detector model is more complete”; it is “we ask a different information question about a compressed scalar observation.”

## Ruling on Question 2

No accessible evidence shows any of the following in Grönberg et al.:

- maximizing scalar count information over `rho` for each finite `nu`;
- a principal information ridge in the `(rho,nu)` plane;
- the scale `rho* = Theta(nu^(1/3))`;
- the expansion
  \[
  \rho^*(\nu)=(6\nu)^{1/3}-\frac23+O(\nu^{-1/3});
  \]
- a complete-data/missing-data decomposition involving the residual detector phase at the frame boundary;
- the identity
  \[
  I_N=\mathbb E[N]-\rho^2\mathbb E[\operatorname{Var}(R_\nu\mid N)].
  \]

Their accessible analysis ends at the count-distribution/moment model, Gaussian performance approximation, and spectral-CT Fisher/CRLB evaluation. ROUND63 continues in a different direction: exact finite-`nu` scalar-count Fisher information, its load maximizer, and a terminal-phase missing-information interpretation.

---

# 3. Best-effort search for closer prior art

## 3.1 Broad operating optima and degradation are established

The manuscript must not claim that dead-time detectors had not previously been assigned an optimal count rate or that information/performance was not known to decline at high flux.

- **Wang et al. (2011)** derived pileup statistics for idealized energy-discriminating nonparalyzable x-ray detectors and evaluated task-specific Cramér–Rao performance versus input count rate; their results show rapid loss of photon-counting advantages as pileup increases. DOI: [10.1118/1.3592932](https://doi.org/10.1118/1.3592932).
- **Alvarez (2014)** used renewal-process large-count moments, multivariate-Gaussian count surrogates, and full/constant-covariance CRLBs to quantify SNR and material-estimation degradation with pileup. DOI: [10.1118/1.4898102](https://doi.org/10.1118/1.4898102).
- **Bécares and Blázquez (2012)** derived an optimal nuclear-counting rate by balancing counting statistics, saturation, source multiplicity, and dead-time-calibration uncertainty. DOI: [10.1155/2012/240693](https://doi.org/10.1155/2012/240693).
- PET/NEMA and radiation-instrumentation literature also contains task-specific maxima such as noise-equivalent count-rate optima.

These are genuine finite-rate optima, but they optimize different task metrics and do not yield the ROUND63 scalar-count ridge law.

## 3.2 Renewal/count-statistics foundations

- **Müller (1973)** reviewed extended and nonextended dead-time models, related distorted interarrival distributions to count statistics through renewal theory, and treated the effect of choosing a random measurement origin. DOI: [10.1016/0029-554X(73)90773-8](https://doi.org/10.1016/0029-554X(73)90773-8).
- **Yu and Fessler (2000)** derived simple exact first and second moments for several dead-time counting models, proposed intensity correction, and assessed when Poisson reconstruction models remain adequate. DOI: [10.1088/0031-9155/45/7/324](https://doi.org/10.1088/0031-9155/45/7/324).
- **Vannucci and Teich (1978)** analyzed rate variation in dead-time-modified Poisson counting, another precursor on non-Poisson dead-time count statistics.

No located version of these works gives the cube-root maximizer or terminal-phase Fisher decomposition.

## 3.3 Time-resolved high-flux information

- **Rapp, Ma, Dawson, and Goyal (2021)** modeled the sequence of dead-time-distorted detection times as a Markov chain and demonstrated high-flux single-photon lidar with more than two orders of magnitude acquisition-time improvement. Their observation retains time stamps/detection sequence information rather than compressing a frame to one count. DOI: [10.1364/OPTICA.403190](https://doi.org/10.1364/OPTICA.403190).
- **Jorgensen and Johnson (2026)** develop asymptotic statistical theory for periodic dead-time event-detection processes, identify sufficient statistics beyond conventional histograms, derive a Fisher-information rate, and prove efficiency of MLE and one-step estimators. Their model is a gated/periodic event sequence and their result is an asymptotic information-rate theory, not the finite single-window integrated-count ridge. Preprint: [arXiv:2605.23210](https://arxiv.org/abs/2605.23210).

These works reinforce—not weaken—the need to state the observation model precisely: timestamps, activation statistics, histograms, and one integrated count preserve different amounts of information.

## 3.4 Search verdict on the two ROUND63 formulas

Targeted searches for combinations of `dead time`, `finite window`, `Fisher information`, `cube root`, `nu^(1/3)`, `terminal phase`, `residual phase`, and `missing information`, together with inspection of the nearest cited/citing literature, found no prior expression matching

\[
\rho^*(\nu)=(6\nu)^{1/3}-\frac23+O(\nu^{-1/3})
\]

or

\[
I_N=\mathbb E[N]-\rho^2\mathbb E[\operatorname{Var}(R_\nu\mid N)].
\]

**Internal provenance verdict:** “no prior matching formula was located” is supported at moderate-high confidence.

**Manuscript discipline:** do not state that no prior exists, do not use “first,” and do not claim that optimal dead-time operation is new in general. State only the narrow positive contribution: **“we derive…”**

---

# 4. Frozen manuscript subsection X.5

The following text is approved verbatim, subject only to replacing citation keys with the manuscript’s final BibTeX keys.

```latex
\subsection{Relation to prior work}

Müller systematized renewal-theoretic analysis of extended and nonextended
dead-time counters, relating distorted interarrival laws to the resulting count
statistics and measurement-origin effects~\cite{Muller1973}. Yu and Fessler
derived exact first and second moments for several dead-time counting models
and used those moments to study intensity correction and the adequacy of
Poisson reconstruction models~\cite{YuFessler2000}. Wang et al. derived
pileup count statistics for idealized energy-discriminating nonparalyzable
x-ray detectors and evaluated task-specific Cramér--Rao performance as a
function of input rate, establishing substantial high-rate degradation of
spectral-counting advantages~\cite{Wang2011}. Alvarez used renewal-process
large-count moments and multivariate-Gaussian count models to quantify
pileup-induced covariance, signal-to-noise loss, and material-estimation
Cramér--Rao bounds~\cite{Alvarez2014}. Grönberg et al. derived an analytical
count-statistics model for nonparalyzable detectors with nonzero pulse length,
validated it by Monte Carlo simulation, and used the model and a Gaussian
approximation to compute spectral-CT Fisher-performance metrics
~\cite{Gronberg2018}.

A complementary line retains richer event information. Rapp et al. modeled
the sequence of dead-time-distorted detection times as a Markov chain, enabling
accurate high-flux single-photon lidar and large acquisition-time reductions
~\cite{Rapp2021}. Jorgensen and Johnson derived asymptotic
Fisher-information rates and efficient estimators for periodic dead-time
event-detection processes, showing that activation and event statistics can
retain information discarded by conventional histogram summaries
~\cite{JorgensenJohnson2026}.

Here we consider a different observation: one finite active-start exposure
from which only the integrated scalar count is retained. For the ideal
nonparalyzable model, we derive the exact finite-window Fisher information of
that count as a function of $\rho=\lambda\tau$ and $\nu=T/\tau$, the
missing-information identity
\begin{equation}
 I_N=\mathbb{E}[N]-\rho^2\mathbb{E}\!\left[
 \operatorname{Var}(R_\nu\mid N)\right],
\end{equation}
and the principal information ridge
\begin{equation}
 \rho^*(\nu)=(6\nu)^{1/3}-\frac{2}{3}+O(\nu^{-1/3}).
\end{equation}
These statements are restricted to the finite-window scalar count-only
observation and do not replace task-specific spectral-CT analyses or
information bounds for time-resolved dead-time event records.
```

## Optional nuclear-counting sentence

If space permits, insert after the Alvarez sentence:

```latex
Task-specific finite-rate optima also arise in nuclear counting; for example,
Bécares and Blázquez optimized counting rate by balancing dead-time losses,
counting uncertainty, source multiplicity, and dead-time-calibration
uncertainty~\cite{BecaresBlazquez2012}.
```

This sentence is useful because it explicitly prevents an overbroad reading of the novelty claim.

---

# 5. Frozen claim discipline

Permitted:

> We derive a finite-window count-information ridge for the active-start nonparalyzable scalar-count observation.

> We derive a terminal-phase missing-information decomposition for the integrated count.

> The ridge differs from task-specific spectral-CT count-rate optima and from information-rate results for time-resolved event sequences.

Not permitted:

- “the first dead-time information optimum”;
- “dead-time information was previously assumed monotone”;
- “prior work did not study optimal count rate”;
- “Grönberg considered only mean and variance” — their official record says they derived a count distribution and computed Fisher-performance metrics;
- “Grönberg’s Fisher analysis was asymptotic only” unless a human obtains the full text and verifies that exact statement;
- “no prior work contains the cube-root law” as a categorical literature fact;
- any suggestion that the ideal zero-pulse-length model supersedes the physically richer nonzero-pulse-length detector model.

The narrow residual claim is sufficient for the paper and survives the accessible prior-art record without priority language.

---

# 6. Sources used for this ruling

- Grönberg et al. 2018 official article record: https://doi.org/10.1002/mp.13063
- Grönberg and Danielsson companion paper: https://doi.org/10.1117/12.2293095
- Müller 1973: https://doi.org/10.1016/0029-554X(73)90773-8
- Yu and Fessler 2000: https://doi.org/10.1088/0031-9155/45/7/324
- Wang et al. 2011: https://doi.org/10.1118/1.3592932
- Alvarez 2014: https://doi.org/10.1118/1.4898102
- Bécares and Blázquez 2012: https://doi.org/10.1155/2012/240693
- Rapp et al. 2021: https://doi.org/10.1364/OPTICA.403190
- Jorgensen and Johnson 2026: https://arxiv.org/abs/2605.23210

**Final status:** X.5 may be filled with the frozen text above and the theory provenance blocker may be closed, with the explicit full-text-access caveat retained in the internal audit record.