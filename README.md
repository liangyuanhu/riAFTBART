# riAFTBART
R codes to implement a flexible approach for causal inferences with multiple treatments and clustered survival outcomes. We develop a flexible random-intercept accelerated failure time model, in which we use Bayesian additive regression trees to capture arbitrarily complex relationships between censored survival times and pre-treatment covariates and use the random intercepts to capture cluster-specific main effects. 
We develop an efficient Markov chain Monte Carlo algorithm to draw posterior inferences about the population survival effects of multiple treatments and examine the variability in cluster-level effects.
We further propose an interpretable sensitivity analysis approach to evaluate the sensitivity of drawn causal inferences about treatment effect to the potential magnitude of departure from the causal assumption of no unmeasured confounding.
Expansive simulation empirically validate and demonstrate good practical operating characteristics of our proposed methods.
The accompanying article is: Hu et al. A flexible approach for causal inference with multiple treatments and clustered survival outcomes arXiv preprint arXiv:2202.08318. 
