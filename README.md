# BayesRPAGP
Bayesian Random Phase-Amplitude Gaussian Process

`BayesRPAGP` in an R software for Bayesian inference of trial-level amplitude, la-
tency, and ERP waveforms. Motivated by the need for rigorous, flexible, and interpretable methods for trial-level analysis of ERP data, we developed the Random Phase-Amplitude Gaussian Process (RPAGP) modeling framework, which assumes that individual trials are generated as a common structural signal transformed by a trial-specific amplitude and phase shift plus ongoing brain activity. 
In the proposed RPAGP framework, the unknown signal is modeled via a Gaussian Process (GP) prior and an autoregressive process is assumed for the ongoing brain activity. We set priors on the trial-specific model parameters and design an efficient algorithm for posterior inference.  Further details are provided in Pluta, D., Hadj-Amar, B., Li., M., Zhao, Y., Versace, F., Vannucci., M,  (2024) ["_Improved Data Quality and Statistical Power of Trial-Level Event-Related Potentials with Bayesian Random-Shift Gaussian Processes"](https://www.tandfonline.com/doi/full/10.1080/01621459.2019.1623043), published in Scientific Reports. 



## Example 

We provide a snapshot of `illustrative_example.Rout`, which contains a tutorial for using our software in Julia

* Run RJMCMC sampler
  ```Julia
  MCMC_simul = AutoNOM(data, hyperparms; s_start = [40])
  ```
<p align="center">
<img src="https://github.com/Beniamino92/AutoNOM/blob/master/figures/posterior_data.png" width="400" heigth="140"/> 
</p>

