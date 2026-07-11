# MICPR: A Multi-Stage Imputation Pipeline for Suppressed Health Data

## The Problem: The "Suppression Dilemma"
In public health research, using CDC WONDER data presents a major hurdle: **Data Suppression**. To protect individual privacy, the CDC hides all death counts between 1 and 9, flagging them as "suppressed."

In my research, this created a "Missing Not at Random" (MNAR) scenario. Listwise deletion would introduce severe selection bias by systematically excluding rural counties from the analysis. Conversely, simple imputation methods—such as zero-replacement or constant value substitution—fail to capture the true underlying data variance, leading to biased and inaccurate regression estimates. My objective was to develop a rigorous, multi-stage approach to reconstruct these data points with the highest possible accuracy, ensuring the imputed values remain representative of the real-world demographic structure.

## The Solution: A Three-Tiered Approach
To solve this, I developed a hybrid pipeline that combines deterministic logic with probabilistic modeling.

### Tier 1: Rule-Based Temporal Substitution (The "Erdman Rules")
I implemented deterministic rules based on the Erdman et al. (2021) framework, which leverages the **temporal stability** of mortality rates within the same geographic unit (county). Since suppression is a binary threshold effect (any count $< 10$ is masked), the surrounding years—the "temporal neighborhood"—provide a robust context to estimate missing values.

* **Conservative Upper Bound:** If both adjacent years (preceding and following) are stable ($\ge 10$ deaths), the missing year is statistically likely to be close to this threshold. I impute this as **9 deaths**. This is a **conservative assumption**: it treats the missing value as the maximum possible count within the suppressed range, avoiding the introduction of non-existent high-frequency data while respecting the censorship threshold.
* **Temporal Interpolation:** If only one neighboring year is stable ($\ge 10$), I apply a linear temporal interpolation, estimating the missing year as 50% of the observed stable neighbor. This assumes a smooth transition in health outcomes over time rather than an abrupt spike or drop, preserving the temporal trend of the county.
* **Gap Preservation:** Where temporal stability could not be confirmed (e.g., when neighbors are also suppressed or data is missing), I retained the `NA` status. This prevents the introduction of artificial noise and delegates these more complex cases to the subsequent probabilistic Tier 3 model.

### Tier 2: The Longitudinal Anchor (6-Year Aggregation)
To bridge the gaps where Tier 1 rules failed, I created a 6-year aggregate dataset. Aggregation breaks the suppression threshold (the sum over 6 years is almost always $> 9$). I used this as a "temporal anchor" to decompose estimates back into 3-year intervals. This ensures that even in small counties, the imputed values are anchored to long-term trends rather than random noise.

### Tier 3: Probabilistic Imputation (MICE & Constrained Recalibration)
For the final missing values, I utilized **Multiple Imputation by Chained Equations (MICE)**:

* **Variance Stabilization:** After analyzing the distribution of non-censored (viable) mortality data, I identified a distinct right-skewed distribution. Consequently, I applied a logarithmic transformation (`log(deaths + 1)`) to the count data. This transformation was specifically selected because it most effectively normalizes the right-skewness observed in the viable data, providing a robust statistical basis for the MICE algorithm.
* **Conditioning:** I used population size as a predictor, which is epidemiologically essential for accurate estimates, as mortality counts are fundamentally tied to the population at risk.
* **The Breakthrough (Recalibration):** A common issue in MICE is "Clamping" (capping values at 9), which creates artificial ceilings. Instead, I introduced **Proportional Min-Max Scaling**. I mapped the imputed estimates back to the $[1, 9]$ range, preserving the relative ordinal ranking between counties while staying strictly within legal privacy boundaries.
