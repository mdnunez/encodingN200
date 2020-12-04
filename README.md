<img src="./extra/small_hnl_logo.png" height="128"> <img src="./extra/small_cidlab_logo.png" height="128">

### Citation

Nunez, M. D., Gosai, A., Vandekerckhove, J., & Srinivasan, R. (2019).
[The latency of a visual evoked potential tracks the onset of decision making.](https://sci-hub.st/https://www.sciencedirect.com/science/article/pii/S1053811919303386) NeuroImage. doi: 10.1016/j.neuroimage.2019.04.052

[Preprint](https://www.researchgate.net/publication/332630466_The_latency_of_a_visual_evoked_potential_tracks_the_onset_of_decision_making)

[Elsevier source](https://www.sciencedirect.com/science/article/pii/S1053811919303386)


# encodingN200  
#### (Repository version 0.5.6)
EEG measures of neural processing reflect visual encoding time (encodingN200) before evidence accumulation during human decision making

**Authors: Michael D. Nunez, Aishwarya Gosai, Joachim Vandekerckhove, and Ramesh Srinivasan from the Cognitive Sciences Department at the University of California, Irvine**

### Significance

Encoding of a sensory stimulus is believed to be the first step in perceptual decision making. In this study, we report evidence that visual evoked potentials (EPs) around 200 ms after stimulus presentation track the time of visual figure-ground segregation before the onset of evidence accumulation during decision making. These EP latencies vary across individuals, are modulated by external visual noise, and increase response time by x milliseconds when they increase by x milliseconds. Hierarchical Bayesian model-fitting was also used to relate these EPs to a specific cognitive parameter that tracks time related to visual encoding in a decision-making model of response time. This work adds to the growing literature that suggests that EEG signals can track the component cognitive processes of decision making.

### Hypothesis

In this study, we directly tested the hypothesis that N200 latencies (EEG stimulus-locked negative peaks occurring between 150 and 275 milliseconds in human subjects) reflect visual processing time that occurs before evidence accumulation during quick decision making.

### Prerequisites

[MATLAB](https://www.mathworks.com/)

[MATLAB Repository: artscreenEEG](https://github.com/mdnunez/artscreenEEG)

[R](https://www.r-project.org/) (for figures)

[Python 2 and Scientific Python libraries](https://www.anaconda.com/products/individual)

For these next 3 install steps in Ubuntu, see [directions here](https://github.com/mdnunez/pyhddmjags/blob/master/jags_wiener_ubuntu.md).

[MCMC Sampling Program: JAGS](http://mcmc-jags.sourceforge.net/)

[Program: JAGS Wiener module](https://sourceforge.net/projects/jags-wiener/)

[Python Repository: pyjags](https://github.com/michaelnowotny/pyjags), can use pip:
```bash
pip install pyjags
```

### Downloading

The repository can be cloned with `git clone https://github.com/mdnunez/encodingN200.git`

The repository can also be may download via the _Download zip_ button above.

### License

encodingN200 is licensed under the GNU General Public License v3.0 and written by Michael D. Nunez, Aishwarya Gosai, Joachim Vandekerckhove, and Ramesh Srinivasan from the Cognitive Sciences Department at the University of California, Irvine.

### Further Reading

Lui, K. K., Nunez, M. D., Cassidy, J. M., Vandekerckhove, J., Cramer, S. C., & Srinivasan, R. (2020).
[Timing of readiness potentials reflect a decision-making process in the human brain.](https://sci-hub.st/https://link.springer.com/article/10.1007/s42113-020-00097-5) Computational Brain & Behavior.

Nunez, M. D., Vandekerckhove, J., & Srinivasan, R. (2017).
[How attention influences perceptual decision making: Single-trial EEG correlates of drift-diffusion model parameters.](https://sci-hub.st/https://www.sciencedirect.com/science/article/abs/pii/S0022249616000316)
Journal of Mathematical Psychology, 76, 117-130.

Nunez, M. D., Srinivasan, R., & Vandekerckhove, J. (2015). 
[Individual differences in attention influence perceptual decision making.](https://www.frontiersin.org/articles/10.3389/fpsyg.2015.00018/full) 
Frontiers in Psychology, 8.

### Posterior samples

Posterior samples for Models 1 and 2 from the paper can be found on [here.](https://figshare.com/articles/Posterior_samples_from_hierarchical_Models_1_and_2_for_the_study_associated_with_the_paper_The_latency_of_a_visual_evoked_potential_tracks_the_onset_of_decision_making_/8244746)

