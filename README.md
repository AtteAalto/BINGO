# BINGO

This is the code for BINGO (Bayesian Inference of Networks with Gaussian prOcess dynamical models), a method for gene 
regulatory network inference from time series data. It is introduced in the article

A. Aalto, L. Viitasaari, P. Ilmonen, L. Mombaerts and J. Goncalves. "Gene regulatory network inference from sparsely sampled noisy data", [_Nature Communications_ 11: 3493 (2020)](https://www.nature.com/articles/s41467-020-17217-1).

Please cite the article if you use the method in your project. Method is tested with Matlab 2017a.


## General information

BINGO assumes that the discrete time series data is sampled from a continuous time trajectory _x_, satisfying _dx/dt = f(x) + u_. The function _f_ is modelled as a Gaussian process (GP) and the network structure is encoded in the hyperparameters of the covariance of the GP. The mean of the _i_<sup> th</sup> component of _f_ is _m_<sub>_i_</sub>(_x_) = _b_<sub>_i_</sub> - _a_<sub>_i_</sub>_x_<sub>_i_</sub>. The method is based on MCMC sampling, meaning that it collects a large number of network samples. If a particular link appears in 60% of the samples, that is the posterior probability for the existence of that link, given the data and the underlying assumptions. It is recommended that at least 5000 samples is collected. Nine out of ten samples are discarded (thinning), so that 5000 samples requires 50000 iterations. It is possible to run a smaller amount of iterations in the beginning to check that the sampler works, and then run more samples. Sampling is fairly slow even for small systems, but the computation time scaling is rather good as the dimension grows.

## Basic use

The code consists of two main functions, _BINGO_init_, that initialises the MCMC sampler and does some data processing. The sampling is done in the function _BINGO_. The inputs to this function are the data, given as a Matlab structure as described below, the state of the sampler (initialised by _BINGO_init_), and the sampling parameters (initialised by the _BINGO_init_, more information below). The output of the function consists of a matrix _Plink_, which is an _n_ x _n_ matrix where _n_ is the number of genes. An element (_i_,_j_) of the matrix _Plink_ tells the number of the collected network samples that contain the link _j_ —> _i_. The function also returns the number of the collected samples, which is called _chain_. The posterior probabilities of links are obtained by dividing _Plink_/_chain_. The output variable _xstore_ contains the posterior mean of the continuous trajectory. The output variable state contains the current state of the sampler (needed if the sampling is continued, for example). The output variable _stats_ contains statistics of the MCMC sampling, which are also displayed in the command window after the sampling is finished.

The main file (_BINGO_main_) used for running the simulations contains first a burn-in, that is, running the sampler for a while and discarding all the samples. Then it contains a sampling part, and finally a part continuing the sampling process. This makes it possible to run some iterations first to get some idea, and then continue the sampling to improve the results. Notice that the method is also easily parallelisable. Just run the method on different processors, and in the end, sum up the _Plink_ matrices. A simple example dataset is provided together with an example file (_BINGO_example_) where the various functionalities of the method (described below in this document) are presented.

>### WARNING!
>
> In the data pre-processing (in the file _BINGO_init_), the time series data is scaled and moved so that the minimum value of each variable is 0 and the maximum is 1. This scaling is equivalent to giving scale-dependent priors to the hyperparameters _beta_. However, if the data contains genes that are not really expressed and contain only noise, this scaling leads to noise amplification. Therefore, the data should be checked and selected beforehand, and genes that are not at all expressed should be discarded. If, however, a gene is expressed in one experiment, and not in another, then it should be kept. Reducing the number of genes will also speed up the algorithm.

## Data structure

The data is given in a Matlab structure called _data_. The time series data is included as a cell field _data.ts_ containing different experiments. Each cell corresponds to one experiment, and should consist of a matrix whose number of rows is the number of genes, and the number of columns is the number of time points in that experiment.

>__Example:__
>
>Number of genes is _n_. Two experiments have been conducted with _m1_ time points in the first experiment, and _m2_ time points in the second experiment. The data is in matrices _y1_ (_n_ x _m1_) and _y2_ (_n_ x _m2_). Then the _data.ts_ field is generated as follows:
>
>`data.ts = {y1 , y2};`

__If your data is uniformly sampled, with no external input, no missing measurements, no prior knowledge of the network, and none of the time series correspond to knockout/knockdown experiments, you are good to go. If this is not the case, read further.__

### Sampling times

If the sampling times are not explicitly given by the user, uniform sampling (with sampling time 1) is assumed. If the sampling is not uniform (or if different sampling frequency is used in the different experiments), this can be specified in the cell field _data.Tsam_. 

>__Example:__
>
>Possible inputs for the cell field _Tsam_ are:  
>  - Constant sampling time intervals (_dt_) for all experiments: `data.Tsam = {dt}`  
>  - Different constant sampling time intervals (_dt1_ and _dt2_) for the two experiments: `data.Tsam = {dt1 , dt2};`  
>  - Completely non-uniform sampling times: `data.Tsam = {[t1_1, t1_2 , … , t1_m1] , [t2_1, t2_2 , … , t2_m2]}`. In this case, the number of given time points must match the size of the corresponding measurement matrix.  
>  - A combination of the two above: `data.Tsam = {[t1 , t2 , … , tm1] , [dt2]};`

__Additional information:__

The main function _BINGO_ wants the field _Tsam_ in the format of a cell field where each cell has all the measurement times (corresponding to the case of completely non-uniform sampling). If the format is something else, then the initialisation file _BINGO_init_ changes the format. However, if the user gives the sampling times in some other form, an error message will be displayed.

### Explicit input (perturbations etc.)

Explicit inputs can be given in the field data.input. Its use is best explained through examples.

>__Example:__
>
>Say _m1_ = _m2_ = 8. Consider the following example cases:  
>  1. The second experiment is a perturbation experiment, where some unknown perturbation is applied half way through the experiment. This kind of input is fed as follows:  
>		`data.input = {zeros(1,8) , [0 0 0 0 1 1 1 1]};`  
>  2. Both experiments are perturbation experiments with two different perturbations. This is given as follows:  
>	  `data.input = {[0 0 0 0 1 1 1 1; zeros(1,8)] , [zeros(1,8); 0 0 0 0 1 1 1 1]};`  
>  3. Both experiments have the same perturbation. This is given as follows:  
>	   `data.input = {[0 0 0 0 1 1 1 1] , [0 0 0 0 1 1 1 1]};`  
>  4. In the first experiment, a perturbation is applied at time zero, and in the other experiment, there is no perturbation:  
>    `data.input = {ones(1,8) , zeros(1,8)};`  
>
>   However, if there is only one experiment with a perturbation applied at time zero, then an input of the form {ones(1,8)} does not have any effect, and it should not be included. In such case, assuming the system was in steady state before the perturbation, it is possible to duplicate the first measurement writing
>
>    `data.ts = {[y(:,1) , y]};`
>
>   and then giving input
>
>    `data.input = {[0 1 1 . . . 1]};`
>
>  5. The input can also be a continuous variable, such as temperature:
>
>  	  `data.input = {[T1_1 , T1_2, …, T1_8] , [T2_1 , T2_2, …, T2_8]};`
>
>    and so on.

When external input is given, the size of the output matrix _Plink_ is _n_ x (_n_ + _n_in_) where _n_in_ is the dimension of the input. The method also tells which genes does the perturbation affect directly. If the perturbation target is known, read the section “Prior knowledge of the network” below.

__Additional information:__

The above warning concerns also external inputs. For example, if the temperature variation, and its effect on the process is negligible, then it should not be included as a possible input. 

### Missing measurements

If the measurement of some genes fail at some time points, this can be included in the information given to the method in a field _data.missing_. Each cell in this field should be a matrix with size _N_miss_ x 2, where each row corresponds to a missing measurement of a gene. The first element on the row identifies the gene in question, and the second element tells the index of the time point where the measurement is missing. The values will be ignored in the method, so it does not matter what is the actual value of the element in the original data matrix.

>__Example:__
>
>Say the measurement of gene 2 failed at the 2nd and 3rd measurement times, and in the second experiment, the measurement of gene 3 failed in the 7th time point and the measurement of gene 5 failed in the first and fourth time points. This is given as follows:
>
> `data.missing = {[2 2; 2 3] , [3 7; 5 1; 5 4]};`

__Additional information:__

Missing measurements are treated so that the Crank—Nicolson sampler draws Brownian bridge samples between the previous and next non-missing sample (corresponding precisely to the measurement model where the missing measurements are omitted). This is achieved by first replacing missing measurements by their linear interpolates (in the file _BINGO_init_), and then increasing the variance of the corresponding measurement noise in the _BINGO_ file.

### Prior knowledge of the network

If some links are certainly known to exist or not to exist, this information can be encoded in the field _data.sure_. This should be an _n_ x _n_ matrix (or _n_ x (_n_ + _n_in_) where _n_in_ is the input dimension if one also has prior knowledge on the targets of potential inputs) with ones in the positions indicating links that are known to exist, minus ones on the places known not to exist, and zeros elsewhere. Such prior information improves the results and speeds up the sampling as the number of candidate networks decreases.

>__Example:__
>
>In a 3 x 3 network, it is known for sure that gene 3 is a regulator of gene 1. In addition, there is one perturbation (external input) which is known to affect only gene 2. This is written as
>
> `data.sure = [0 0 1 -1; 0 0 0 1; 0 0 0 -1];`

Links can also be given different Bayesian prior probabilities of existence. The prior proabilities are given in a field _parameters.link_pr_. This should be a matrix of size _n_ x _n_ matrix (or _n_ x (_n_ + _n_in_) where _n_in_ is the input dimension. If a certain link exists with (prior) probability _p_, then the corresponding element in _link_pr_ matrix is  _p_ / (1-_p_). If a link is known to exist or to not exist for sure, it is advised to use the _data.sure_ field, since that will speed up the MCMC sampling.






### Knockout time series

If some time series correspond to gene knockout/knockdown experiments, then that time series should not be taken into account in the inference of links pointing to the knocked-out/down gene. Such information is included in the field _data.ko_. For example, if there are four time series, and the third corresponds to an experiment where the second gene is knocked out and the fourth time series corresponds to an experiment where the fourth and fifth genes are knocked out, it is given as follows:

`data.ko = {[ ] , [ ] , [2] , [4 5]};`

## Method parameters

Method parameters are given to the _BINGO_-function as a structure called parameters. The values are set in the file _BINGO_init_, and it is not necessary to change these parameters. The parameter fields are:  

  - _link_pr_: The prior for networks is that the existence of each link is independent of others. If a link exists with (prior) probability _p_, then _link_pr_ = _p_ / (1-_p_). This parameter controls the level of sparsity in the network samples. This can be specified by the user, but it is not necessary. For example, if there are 50 nodes, and it is expected that each node has about 2-3 (incoming/outgoing) links, then one should set, for example, _p_ = 2.5 / 50. However, if the user does not provide the parameter, the default value is _link_pr_ = 1/_n_. Note that this parameter can also be a matrix, if different priors are given for different links.      
  - _its_: number of iterations.   
  - _nstep_: number of steps in the “continuous” trajectory between the measurement times. Default value is 4. It can be reduced to 3 or even 2 if the sampling takes a very long time. If it is one, the method reduces to a kind of a discrete-time GPDM.  
  - _Theur_: a temperature variable used in a heuristic scheme to speed up the topology sampling. If _Theur_ is one, it corresponds to exact sampling. Higher _Theur_ improves the acceptance probability of the topology sampling.
    
A number of different __step size parameters__ can be set by the user. Shorter step sizes increase the acceptance rates of the MCMC sampling, but it has a negative effect on the efficiency of the exploration of the parameter space. Note that hyperparameters _gamma_, _beta_, _a_, and _b_ are sampled together. If their acceptance probability drops, either one or more of these can have too long step sizes. Similarly, the trajectory is sampled together with the noise variance _q_.  

  - _etraj_: step size for the Crank—Nicolson sampler for the continuous trajectory. Default value is 0.06, and it is recommended that this is at least 0.05. It can be adjusted if the acceptation probability for the trajectory becomes low, which may happen if there are several experiments in the data.  
  - _ea_: step size for sampling the _a_-parameters in the mean function.  
  - _eb_: step size for sampling the _b_-parameters in the mean function.  
  - _egamma_: step size for sampling the _gamma_-parameter in the GP covariance.  
  - _ebeta_: step size for sampling the _beta_-parameter in the GP covariance.  
  - _er_: step size for sampling the measurement error variance _r_.  
  - _eq_: step size for sampling the process noise covariance _q_.  

Not part of the parameters-structure is _nr_pi_, which is the number of pseudoinputs used in the GP regression in the code. Higher number improves accuracy, but slows down sampling.

## Troubleshooting

__Problem:__ The acceptance probabilities for the trajectory and the topology are (almost) zero.

__Possible solutions:__ The step size of the trajectory (_parameters.etraj_) can be reduced, but it should not be below 0.05. It is possible that the initial value for the process noise covariance is too small. Try to increase it by changing it inside the file _BINGO_init_ (_state.q_). There is also a heuristic temperature variable (_parameters.Theur_) which is one by default. It can be increased to about 1.5 to speed-up topology sampling. However, this method is heuristic, and then the sampling is not exactly what is claimed.

__Problem:__ The method is very slow.

__Possible solutions:__ 
  1. Consider parallelisation (see below).
  2. The gene expression trajectory is sampled on a finer grid where each measurement interval is divided into four parts (by default). This can be made coarser by changing _nstep_ = 3 inside the function _BINGO_init_.
  3. Number of pseudoinputs can be reduced inside _BINGO_init_ by setting _nr_pi_ = 40 (for example, default is 50).
  4. Reduce dimension by throwing out genes whose expression data contains merely noise (which actually should be done anyway). 
  5. If some genes are only interesting as potential regulators, they can be given to the method as inputs rather than as part of the time series data.

## Parallelisation

  1. Generate the _data_-structure;
  
`parfor`-loop  

  2. Run _BINGO_init_;  
  3. The burn-in can be shortened by setting _parameters.its_ = 1000 (for example, default is 3000);  
  4. Run _BINGO_ with a suitable number of iterations;  
  5. Store results;  
`end`

  6. Sum up the different _Plink_ matrices and divide by the total number of samples.

  ## Updates on the method since publication of the article
  
  __March 22, 2021:__ Hyperparameter sampling is changed from random walk to random walk in log-domain to make it scale-independent.

  __April 21, 2021:__ Prior probabilities of links can be different.



