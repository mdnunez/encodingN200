<img src="./extra/small_hnl_logo.png" height="128"> <img src="./extra/small_cidlab_logo.png" height="128">

# encodingN200 0.1.2
EEG measures of neural processing reflect human visaual encoding time (encodingN200).

**Authors: Michael D. Nunez, Aishwarya Gosai, Joachim Vandekerckhove, and Ramesh Srinivasan from the Cognitive Sciences Department at the University of California, Irvine**

### Significance

Encoding of a sensory stimulus is the first step in perceptual decision making. From this analysis we report evidence that human evoked potentials (EPs) around 200 milliseconds after stimulus presentation track the time of extrastriate cortical processing before visual evidence accumulation during decision making. These EP latencies vary across individuals and depend upon external visual noise. We present linear regressions, framed in the context of cognitive theory, to test 1-to-1 relationships between these EP latencies and reaction time distributions. Hierarchical Bayesian model-fitting was also used to relate these EPs to a specific cognitive parameter that tracks time related to visual encoding. This work adds to the growing literature that suggests that EEG records can track the onset of evidence accumulation.

### Hypothesis

Non-decision times for each session of data, as estimated by drift-diffusion models of accuracy and reaction time distributions, will have a 1-to-1 linear relationship with N200 negative peak latency time (an event-related potential possibly reflecting visual encoding time)

### Prerequisites

[MATLAB](https://www.mathworks.com/)

[MATLAB Repository: artscreenEEG](https://github.com/mdnunez/artscreenEEG)

[MCMC Sampling Program: JAGS](http://mcmc-jags.sourceforge.net/)

[Program: JAGS Wiener module](https://sourceforge.net/projects/jags-wiener/)

[Scientific Python libraries](https://www.continuum.io/downloads)

[Python Repository: pyjags](https://github.com/tmiasko/pyjags)

[R](https://www.r-project.org/) (for figures)

### Downloading

The repository can be cloned with `git clone https://github.com/mdnunez/encodingN200.git`

The repository can also be may download via the _Download zip_ button above.

### Installation

After downloading/unzipping the repository, users will need to add these functions to the MATLAB path. In MATLAB, add the repository to the PATH with

```matlab
%Set 'artloc' to full directory path
emloc = 'C:\Users\MATLAB\encodingN200';
addpath(genpath(emloc));
```

### License

encodingN200 is licensed under the GNU General Public License v3.0 and written by Michael D. Nunez, Aishwarya Gosai, Joachim Vandekerckhove, and Ramesh Srinivasan from the Cognitive Sciences Department at the University of California, Irvine.

### Further Reading

Nunez, M. D., Vandekerckhove, J., & Srinivasan, R. (2017).
[How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters.](https://www.researchgate.net/publication/298275031_How_attention_influences_perceptual_decision_making_Single-trial_EEG_correlates_of_drift-diffusion_model_parameters)
Journal of Mathematical Psychology, 76, 117-130.

Nunez, M. D., Srinivasan, R., & Vandekerckhove, J. (2015). 
[Individual differences in attention influence perceptual decision making.](https://www.researchgate.net/publication/273466831_Individual_differences_in_attention_influence_perceptual_decision_making) 
Frontiers in Psychology, 8.

