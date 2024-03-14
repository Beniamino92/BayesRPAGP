# BayesRPAGP
Bayesian Random Phase-Amplitude Gaussian Process

`BayesRPAGP` in an R software for Bayesian inference of trial-level amplitude, latency, and ERP waveforms. Motivated by the need for rigorous, flexible, and interpretable methods for trial-level analysis of ERP data, we developed the Random Phase-Amplitude Gaussian Process (RPAGP) modeling framework, which assumes that individual trials are generated as a common structural signal transformed by a trial-specific amplitude and phase shift plus ongoing brain activity. 
In the proposed RPAGP framework, the unknown signal is modeled via a Gaussian Process (GP) prior and an autoregressive process is assumed for the ongoing brain activity. We set priors on the trial-specific model parameters and design an efficient algorithm for posterior inference.  Further details are provided in Pluta, D., Hadj-Amar, B., Li., M., Zhao, Y., Versace, F., Vannucci., M,  (2024) ["Improved Data Quality and Statistical Power of Trial-Level Event-Related Potentials with Bayesian Random-Shift Gaussian Processes"](https://www.provanonesiste2231123.com), published in Scientific Reports. 



## Example 

We provide a snapshot of `tutorial.Rmd`, which contains a tutorial for using our software in R

* Run RAPGP sampler
  ```R
  results <- fit_rpagp(y = dat$y, n_iter = n_iter,
                         theta0 = theta0, hyperparam = hyperparam,
                         pinned_point = pinned_point,
                         pinned_value = pinned_value)
  ```
<p align="center">
<img src="https://github.com/BayesRPAGP/blob/main/plots/example.png" width="400" heigth="140"/> 
</p>

