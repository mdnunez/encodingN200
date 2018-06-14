# pdm5b_papermodels.py - Contains hierarchical models used in manuscript
#
# Copyright (C) 2018 Michael D. Nunez, <mdnunez1@uci.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Record of Revisions
#
# Date            Programmers                         Descriptions of Change
# ====         ================                       ======================
# 12/04/17      Michael Nunez                         Cleanup for release
# 12/12/17      Michael Nunez                           Add basic model
# 12/28/17      Michael Nunez                           Add linear models
# 01/09/18      Michael Nunez                          Clear up language
# 01/26/18      Michael Nunez                       Adding 'all_n200lat_random'
# 06/13/18      Michael Nunez     Add model with properly accounted-for lapse trials

# ### JAGS Models
jagsmodels = dict()

# Linear regression of trial-averaged or single-trial N200 latencies on reaction time percentiles
jagsmodels['n1lat_regression'] = '''
model {
	##########
    #Prior of Regression intercept (remaining reaction time)
    tercond ~ dnorm(.3, pow(.25,-2))
    #Prior of effect of single-trial N200 latency on reaction time
    n1gammault ~ dnorm(1,pow(3,-2)) #'Informative' prior for Bayes Factor calculation
    #Prior of Variability in reaction time
    tersubsd ~ dgamma(.2,1)

    ##########
    #Linear regression (normal likelihood)
    for (i in 1:N) {
        rtpercentile[i] ~ dnorm(tercond + n200lat[i]*n1gammault, pow(tersubsd,-2) )
    }
}
'''

# 3 Parameter Model with random effects of true N1 on all parameters across EEGsessions, split by experiment
jagsmodels['all_n1lat_random'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-EEGsession variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in N1 latency
	n1subsd ~ dgamma(.2,1)

	for (f in 1:3) {
		n1gammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		n1gammasd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {
		for (e in 1:nexps) {
			#Condition-level N1 latency
			n1cond[e,k] ~ dnorm(.2, pow(.1,-2))

	        #Condition-level non-decision time
			tercond[e,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level drift rate
			deltacond[e,k] ~ dnorm(1, pow(2, -2))

			#Condition-level boundary separation
			alphacond[e,k] ~ dnorm(1, pow(.5,-2))

			#Condition-level effects of N200 latency on non-decision time, drift rate, and boundary separation
			for (f in 1:3) {
				n1gammacond[f,e,k] ~ dnorm(n1gammault[1,f],pow(n1gammasd[1,f],-2))
			}
		}
		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[experiment[ses]+1,k]
				+ n1gammacond[1,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[1,2,k]*experiment[ses]*n1sub[k,ses],
				pow(tersubsd, -2))T(0,1)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[experiment[ses]+1,k]
				+ n1gammacond[2,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[2,2,k]*experiment[ses]*n1sub[k,ses],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[experiment[ses]+1,k]
				+ n1gammacond[3,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[3,2,k]*experiment[ses]*n1sub[k,ses],
				pow(alphasubsd, -2))T(.1,3)

		    #EEGsession-level N1 latency
			n1sub[k,ses] ~ dnorm(n1cond[experiment[ses]+1,k],
				pow(n1subsd, -2))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphasub[condition[i],EEGsession[i]],
		tersub[condition[i],EEGsession[i]],
		beta,
		deltasub[condition[i],EEGsession[i]])
	}
}
'''

# 3 Parameter Model with random effects of non-decision time on true N1 across EEGsessions, split by experiment, additive effect of noise
jagsmodels['all_n1lat_request2'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-session variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-session variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-session variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-session variability in N1 latency
	n1subsd ~ dgamma(.2,1)

	##########
	#Block-level parameters
	##########

	#Condition-level effects of N1 latency on non-decision time
	for (f in 1:3) {
		n1gammacond[f,1,3] ~ dnorm(1,pow(3,-2))
		n1gammacond[f,2,3] ~ dnorm(0,pow(1,-2))
		for (e in 1:nexps) {
			for (k in 1:2) {
				n1gammacond[f,e,k] ~ dnorm(0,pow(1,-2))
			}
		}
	}

	for (k in 1:nconds) {
		for (e in 1:nexps) {
			#Condition-level N1 latency
			n1cond[e,k] ~ dnorm(.2, pow(.1,-2))

	        #Condition-level non-decision time
			tercond[e,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level drift rate
			deltacond[e,k] ~ dnorm(1, pow(2, -2))

			#Condition-level boundary separation
			alphacond[e,k] ~ dnorm(1, pow(.5,-2))

		}
		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[experiment[ses]+1,k]
				+ n1gammacond[1,1,3]*n1sub[k,ses]
				+ n1gammacond[1,1,2]*(k<3)*n1sub[k,ses]
				+ n1gammacond[1,1,1]*(k<2)*n1sub[k,ses]
				+ n1gammacond[1,2,3]*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[1,2,2]*(k<3)*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[1,2,1]*(k<2)*experiment[ses]*n1sub[k,ses],
				pow(tersubsd, -2))T(0,1)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[experiment[ses]+1,k]
				+ n1gammacond[2,1,3]*n1sub[k,ses]
				+ n1gammacond[2,1,2]*(k<3)*n1sub[k,ses]
				+ n1gammacond[2,1,1]*(k<2)*n1sub[k,ses]
				+ n1gammacond[2,2,3]*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[2,2,2]*(k<3)*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[2,2,1]*(k<2)*experiment[ses]*n1sub[k,ses],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[experiment[ses]+1,k]
				+ n1gammacond[3,1,3]*n1sub[k,ses]
				+ n1gammacond[3,1,2]*(k<3)*n1sub[k,ses]
				+ n1gammacond[3,1,1]*(k<2)*n1sub[k,ses]
				+ n1gammacond[3,2,1]*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[3,2,2]*(k<3)*experiment[ses]*n1sub[k,ses]
				+ n1gammacond[3,2,3]*(k<2)*experiment[ses]*n1sub[k,ses],
				pow(alphasubsd, -2))T(.1,3)

		    #EEGsession-level N1 latency
			n1sub[k,ses] ~ dnorm(n1cond[experiment[ses]+1,k],pow(n1subsd, -2))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphasub[condition[i],EEGsession[i]],
		tersub[condition[i],EEGsession[i]],
		beta,
		deltasub[condition[i],EEGsession[i]])
	}
}
'''

# Random effects of single-trial N200 latency on all parameters, split by experiment
jagsmodels['all_n200lat_random'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-subject variability in residual non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-subject variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-subject variability in log() boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-subject variability in N200 latency
	n200subsd ~ dgamma(.2,1)

	#Between-trial variability in N200 latency
	n200trialsd ~ dgamma(.2,1)

	for (f in 1:3) {

		#Hierarchical-level effect of N200 latency
		n1gammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior

		#Between-condition variability in effect of N200 latency
		n1gammasd[1,f] ~ dgamma(1,1)

		#Between-subject variability in N200 effect
		n1gammasubsd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {
		for (e in 1:nexps) {
			#Condition-level N200 latency
			n200cond[e,k] ~ dnorm(.2, pow(.1,-2))

	        #Condition-level residual non-decision time
			tercond[e,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level residual drift rate
			deltacond[e,k] ~ dnorm(1, pow(2, -2))

			#Condition-level residual log() boundary separation
			alphacond[e,k] ~ dnorm(0, pow(.7,-2))

			#Condition-level effects of N200 latency on non-decision time, drift rate, and boundary separation
			for (f in 1:3) {
				n1gammacond[f,e,k] ~ dnorm(n1gammault[1,f],pow(n1gammasd[1,f],-2))
			}
		}
		
		#Subject-level parameters
		for (ses in 1:nses) {
			#Subject-level residual non-decision time
			tersub[k,ses] ~ dnorm(tercond[experiment[ses]+1,k], 
				pow(tersubsd, -2))

			#Subject-level residual drift rate
			deltasub[k,ses] ~ dnorm(deltacond[experiment[ses]+1,k], 
				pow(deltasubsd, -2))

			#Subject-level diffusion coefficient
			alphasub[k,ses] ~ dnorm(alphacond[experiment[ses]+1,k], 
				pow(alphasubsd, -2))

		    #Subject-level N200 latency
			n200sub[k,ses] ~ dnorm(n200cond[experiment[ses]+1,k], 
				pow(n200subsd, -2))

			#Subject-level effects of N200 latency on non-decision time, drift rate, and boundary separation
			for (f in 1:3) {
				n1gammasub[f,k,ses] ~ dnorm(n1gammacond[f,experiment[ses]+1,k],
					pow(n1gammasubsd[1,f],-2))
			}
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		#Note that N200 latencies are censored between 150 and 274 ms
		n200lat[i] ~ dnorm(n200sub[condition[i], EEGsession[i]],
			pow(n200trialsd,-2))T(.151,.274)

		y[i] ~ dwiener(exp(alphasub[condition[i],EEGsession[i]]
		+ n1gammasub[3, condition[i], EEGsession[i]]*n200lat[i]),
		tersub[condition[i],EEGsession[i]]
		+ n1gammasub[1, condition[i],EEGsession[i]]*n200lat[i], 
		beta, 
		deltasub[condition[i],EEGsession[i]]
		+ n1gammasub[2, condition[i],EEGsession[i]]*n200lat[i])
	}
}
'''

# 3 Parameter Model with random effects of true N1 on all parameters across EEGsessions, split by experiment
#  Properly accounting for lapse trials
jagsmodels['all_n1lat_random_lapse'] = '''
model {
	##########
	#Fixed Parameters
	beta <- .5
	##########
	#Between-EEGsession variability in non-decision time
	tersubsd ~ dgamma(.2,1)

	#Between-EEGsession variability in drift
	deltasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in boundary separation
	alphasubsd ~ dgamma(1,1)

	#Between-EEGsession variability in N1 latency
	n1subsd ~ dgamma(.2,1)

	for (f in 1:3) {
		n1gammault[1,f] ~ dnorm(1,pow(3,-2)) #'Informative' prior
		n1gammasd[1,f] ~ dgamma(1,1)
	}

	##########
	#Block-level parameters
	##########
	for (k in 1:nconds) {
		for (e in 1:nexps) {
			#Condition-level N1 latency
			n1cond[e,k] ~ dnorm(.2, pow(.1,-2))

	        #Condition-level non-decision time
			tercond[e,k] ~ dnorm(.3, pow(.25,-2))

			#Condition-level drift rate
			deltacond[e,k] ~ dnorm(1, pow(2, -2))

			#Condition-level boundary separation
			alphacond[e,k] ~ dnorm(1, pow(.5,-2))

			#Condition-level effects of N200 latency on non-decision time, drift rate, and boundary separation
			for (f in 1:3) {
				n1gammacond[f,e,k] ~ dnorm(n1gammault[1,f],pow(n1gammasd[1,f],-2))
			}
		}
		#EEGsession-level parameters
		for (ses in 1:nses) {
			#EEGsession-level non-decision time
			tersub[k,ses] ~ dnorm(tercond[experiment[ses]+1,k]
				+ n1gammacond[1,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[1,2,k]*experiment[ses]*n1sub[k,ses],
				pow(tersubsd, -2))T(0,1)

			#EEGsession-level drift rate
			deltasub[k,ses] ~ dnorm(deltacond[experiment[ses]+1,k]
				+ n1gammacond[2,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[2,2,k]*experiment[ses]*n1sub[k,ses],
				pow(deltasubsd, -2))T(-9, 9)

			#EEGsession-level boundary separation
			alphasub[k,ses] ~ dnorm(alphacond[experiment[ses]+1,k]
				+ n1gammacond[3,1,k]*(1-experiment[ses])*n1sub[k,ses]
				+ n1gammacond[3,2,k]*experiment[ses]*n1sub[k,ses],
				pow(alphasubsd, -2))T(.1,3)

		    #EEGsession-level N1 latency
			n1sub[k,ses] ~ dnorm(n1cond[experiment[ses]+1,k],
				pow(n1subsd, -2))

			#EEGsession-level lapse trials
        	probsub[k, ses, 1:2] ~ ddirch(c(1,1))
		}
	}
	##########
	#Wiener likelihoods
	for (i in 1:N) {
		y[i] ~ dwiener(alphasub[condition[i],EEGsession[i]],
		tersub[condition[i],EEGsession[i]],
		beta,
		deltasub[condition[i],EEGsession[i]])
	}
	##########
	# Wiener likelihoods
	for (i in 1:N) {
        # Log density for DDM process
        ld_comp[i, 1] <- dlogwiener(y[i],alphasub[condition[i],EEGsession[i]],
        tersub[condition[i],EEGsession[i]],
        beta,
        deltasub[condition[i],EEGsession[i]])

        # Log density for lapse trials (negative max RT to positive max RT)
        ld_comp[i, 2] <- logdensity.unif(y[i], -maxrt[condition[i],EEGsession[i]], maxrt[condition[i],EEGsession[i]])

        # Select one of these two densities (Mixture of nonlapse and lapse trials)
        density[i] <- exp(ld_comp[i, component_chosen[i]] - Constant)
		
        # Generate a likelihood for the MCMC sampler using a trick to maximize density value
        Ones[i] ~ dbern(density[i])

        # Probability of mind wandering trials (lapse trials)
        component_chosen[i] ~ dcat(probsub[condition[i],EEGsession[i],1:2])
	}
}
'''