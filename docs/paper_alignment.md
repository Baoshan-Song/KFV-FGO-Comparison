# Paper-to-code audit

Reviewed versions: [v1](https://arxiv.org/html/2511.00306v1) and [v2](https://arxiv.org/html/2511.00306v2). The implementation target is v2. The archived MATLAB branch and the former Python tip `fe1a3a5` are not numerically equivalent to the corrected recurrence.

## Mapping

| Paper requirement | Implementation | Verification |
|---|---|---|
| v1 §3.1.2 single-state, predict/update structure | Re-FGO window 1; current observations follow state elimination | Window guard and stage-one tests |
| v2 Eqs. (4), (17), predictive Schur complement | `FactorGraph.marginalize`, local square-root elimination | Compare with `inv(F P F.T + Q)` |
| v2 Eqs. (18)–(20), fixed prediction during iteration | Anchored `MarginFactor`; gain-form `iterated_update` | Linear analytic posterior; nonlinear independent solvers |
| v2 Eq. (5), Table 3, factor-only QR | `mar_measurements`, retaining `R`, `d`, and anchor | Hessian, affine RHS and total cost checks |
| Historical measurements cannot relinearize after anchoring | Retire factors, keep current state and one linear prior | Factor lifecycle checks |

## What changed in this implementation

The previous Python branch deliberately froze a prior's residual as its state moved and multiplied measurement information by another robust weight on every evaluation. Its second-stage prior was constructed from the estimated state and a Hessian block, with an identity fallback. A second, stale public package could bypass graph optimization entirely. All four behaviors have been removed.

The implementation uses the sign convention `A δ = b`, followed by `x += δ`. If `a` is the linearization anchor, the stored prior is

$$\|R(x-a)-d\|^2.$$

At another state, its right-hand side is `d − R(x − a)`, not a frozen `d`. The stage-two code computes thin QR `A = Q R`, then `d = Qᵀb`. For any increment,

$$\|A\delta-b\|^2=\|R\delta-d\|^2+\|b-Qd\|^2.$$

The final term is state independent and is saved as a diagnostic. Tests verify both the information matrix and affine term, so a zero-residual prior made from an arbitrary estimated state cannot pass. Cholesky/QR factors are used with their actual orientation; no symmetric square root is assumed.

## Finite iterations and convergence

V2 describes anchoring at the converged mode. A numerical solver only approaches that mode within a tolerance. Python stores the **last linear system that produced the accepted update**, with its original anchor. At convergence this implements mode anchoring up to the stopping tolerance. With one iteration it retains the prediction-point Jacobian and therefore gives EKF. With an exhausted multi-iteration budget it implements the same truncated IEKF update, but does not assert convergence. Re-evaluating the Jacobian at the posterior only on the graph side would break this correspondence.

The original nonlinear factors are never evaluated again after stage two. The next stage-one prior absorbs the stored affine factor and the new motion factor. Ordinary sliding-window FGO is intentionally different: active nonlinear factors can be relinearized until their states are retired.

## Robust extension and scope

The Gaussian fixed-prior correspondence is the primary verification target. Robust modes use fixed nominal range variance and recompute IRLS weights from each current standardized residual. The prediction remains Gaussian and unchanged within an epoch. Both independent solvers use this same explicitly defined extension; the reported covariance is its local weighted least-squares approximation, not an exact covariance under a non-Gaussian posterior. The previous cumulative-weight MATLAB recurrence is not preserved.

The old `autoDiff` option actually used finite differences and was ignored by one caller because of a naming mismatch. It now works and is labelled correctly. Generic finite-difference factors and the displacement factor also had inconsistent residual/Jacobian signs; these are corrected and tested.

The generated circular dataset now uses `angle[k] = omega * k * dt`, avoiding the duplicated endpoint and inconsistent angular interval. The motion function retains the existing first-order discretization, so truth/model discretization error still exists. Loaded MAT files are not rewritten. The simulation check is distinct from the real-data workflow described below.

## Real GNSS/INS contribution

PR #5 supplies RINEX/Xsens readers, a ten-state pseudorange/IMU model, UTC alignment and saved MATLAB references. These inputs now run through the same canonical graph as the simulation. The real filter uses a gain update with a fixed prediction and recomputes nominal elevation/SNR variances and robust weights at each linearization. Re-FGO independently assembles these observation factors and performs both QR marginalizations. No public import path retains the old filter shortcut or frozen-residual prior.

The contribution's first-valid-sample Euler propagation, supplied AHRS attitude, discrete process covariance and observation corrections are retained. This is not a full strapdown INS or RTK implementation. Model probes and the one-update EKF still reproduce their saved MATLAB checks. Legacy tests requiring filter delegation, an unanchored frozen prior or the old sliding-window objective were replaced with independent state/covariance and causal-residual checks, rather than preserving those bugs. Saved references and provenance remain unchanged.

[Real validation](real_validation.json) reports all four matched pairs over 438 epochs. The [Colab guide](../README.md) documents model limits and parameter meanings. Real GUI accuracy uses interpolated reference ECEF positions at overlapping UTC times, with no fitted spatial alignment. Article figures use the same viewer and online histories.
