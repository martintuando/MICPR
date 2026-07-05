# MICPR: A Multi-Stage Imputation Pipeline for Suppressed Health Data

## The Problem: The "Suppression Dilemma"
In public health research, using CDC WONDER data presents a major hurdle: **Data Suppression**. To protect individual privacy, the CDC hides all death counts between 1 and 9, flagging them as "suppressed."

In my research, this created a "Missing Not at Random" (MNAR) scenario. If I simply deleted these rows, I would lose most rural counties (Selection Bias). If I guessed (e.g., replaced with 0 or 5), I would distort the statistical variance. I needed a way to recover this data without making up values out of thin air.

## The Solution: A Three-Tiered Approach
To solve this, I developed a hybrid pipeline that combines deterministic logic with probabilistic modeling.

### Tier 1: Rule-Based Temporal Substitution (The "Erdman Rules")
I implemented deterministic rules based on the Erdman et al. (2021) framework:
* **Stability Check:** If surrounding years show $\ge 10$ deaths, the missing year is filled with 10.
* **Proportional Filling:** If one neighbor is stable, I estimate the missing value as half of that neighbor.
* **Gap Preservation:** Where rules didn't apply, I left the values as `NA` to allow the statistical model to handle them later.

### Tier 2: The Longitudinal Anchor (6-Year Aggregation)
To bridge the gaps where yearly rules failed, I created a 6-year aggregate dataset. Aggregation breaks the suppression threshold (the sum over 6 years is almost always $> 9$). I used this as a "temporal anchor" to decompose estimates back into 3-year intervals. This ensures that even in small counties, the imputed values are anchored to long-term trends rather than random noise.

### Tier 3: Probabilistic Imputation (MICE & Constrained Recalibration)
For the final missing values, I utilized **Multiple Imputation by Chained Equations (MICE)**:
* **Variance Stabilization:** I log-transformed the counts (`log(deaths + 1)`) to handle right-skewed data.
* **Conditioning:** I used population size as a predictor, which is epidemiologically essential for accurate estimates.
* **The Breakthrough (Recalibration):** A common issue in MICE is "Clamping" (capping values at 9), which creates artificial ceilings. Instead, I introduced **Proportional Min-Max Scaling**. I mapped the imputed estimates back to the $[1, 9]$ range, preserving the relative ordinal ranking between counties while staying strictly within legal privacy boundaries.

## Scientific Value
This process transforms "Missing Data" into "Reconstructed Reality." By using this pipeline, the final dataset retains its demographic structure and statistical variance, preventing the systematic bias that would occur in a standard analysis.

## Usage
The script `cdc_imputation.R` handles the end-to-end process: loading the panels, applying the 3-tier logic, and exporting the final `Master_Imputed_Outcomes_Wide.csv` ready for regression analysis.
