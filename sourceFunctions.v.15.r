# load libraries
# To track the fitness landscape it has become practical to have SQL lite databases to hold this object
# It is only truly necessary with genotypes longer than ~ 30, and thus considering this programs intent...
library("RSQLite")
library("DBI")
library("sn") # This allows the skewwed normal distribution
library("evd") # This allows the extreme value distributions to be used.
library("VGAM") # This includes the Fretchet distribution to be called
library("compiler") # This allows me to pre-compile byte code functions.

###################################################################################################
####################################### BEGIN OF FUNCTIONS ########################################
###################################################################################################

# This allows us to simulate the birth process using a deterministic form, then to round to integer values if requested
birthFunction <- function(func_inSize, func_inProb = 1, func_roundValues = TRUE){
	# Calcualte the deterministic value and then update if we are rounding.
	func_tmpReturn <- func_inSize * func_inProb
	if(func_roundValues && any(func_tmpReturn %% 1 != 0)){
		func_tmpUpdates <- which(func_tmpReturn %% 1 != 0)
		func_tmpReturn[func_tmpUpdates] <- unlist(lapply(func_tmpReturn[func_tmpUpdates],function(thisPop){ 
																floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1) 
														}))
	}
	return( func_tmpReturn )
}

# This allows us to simulate the death process as a deterministic value, and may be density dependent.  If 
# the death rate is density dependent then a value must be passed to func_depDensity or we simply use the full func_inProb value
deathFunction <- function(func_inSize, func_inProb = 0, func_roundValues = TRUE,
							func_depDensity = FALSE, func_densityMax = NULL, func_densityPower = 4){
	# The number of death events defined deterministically using the probability of death and the popualion size, with rounding being 
	# the only stochastic portion.  However, we can ignore any zero value populations, they won't have further deaths

	# Here we check if the calculation is density dependent or not
	func_tmpReturn <- if(func_depDensity && !is.null(func_densityMax)){
					func_inSize * func_inProb * (sum(func_inSize)/func_densityMax)^ func_densityPower
				} else {
					func_inSize * func_inProb
				}
							
	# this is where the stochastic rounding is performed.
	if(func_roundValues && any(func_tmpReturn %% 1 != 0)){
		func_tmpUpdates <- which(func_tmpReturn %% 1 != 0)
		func_tmpReturn[func_tmpUpdates] <- unlist(lapply(func_tmpReturn[func_tmpUpdates],function(thisPop){ 
																floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1) 
														}))
	}
	return( func_tmpReturn )
	#return( unlist(lapply(tmpReturn,function(thisPop){ floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1) })) )
}

# This allows us to simulate the mutation process as a deterministic value, at present I'm not permitting anything but integer results
mutationFunction <- function(func_inSize, func_inProb = 0){
	# Calcualte the deterministic value and then update if we are rounding.
	func_tmpReturn <- func_inSize * func_inProb
	# Here is where I enforce updating of value to be integer.
	func_tmpUpdates <- which(func_tmpReturn %% 1 != 0)
	func_tmpReturn[func_tmpUpdates] <- unlist(lapply(func_tmpReturn[func_tmpUpdates],function(thisPop){ 
																floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1) 
														}))

	return( func_tmpReturn )
}

# This is a simple little function used to represent drift (through poisson distribution calls) to a vector
# It is possible to add drift via a call to the gamma function, and by using the shape value to redefine the 
# rate and shape parameters 
addDrift <- function(func_inVector, func_integerValues = TRUE){
	# At present I don't have a proper continuous Poisson like distribution from which to sample
	# Thus I suppress the ability to not round values.
	func_integerValues = TRUE ### THIS EXISTS TO FORCE RPOIS AS A COMPARISON DUE TO LACK OF CONTINUOUS POISSON LIKE DISTRIBUTION TO SAMPLE
	
	# If we're rounding to integer values then we use the poisson distribution, else we use gamma with default rate parameter of 1
	sign(func_inVector) * unlist(lapply(abs(func_inVector),function(thisValue){ 
									if(func_integerValues){
										rpois(1,thisValue)
									} else {
										### TEST HAVE SHOWN THIS TO NOT WORK WELL
										func_tmpShape <- thisValue
										func_tmpRate <- 1 + (1.6/func_tmpShape)
										rgamma(1, func_tmpShape * func_tmpRate, func_tmpRate)
									}
								}))
}

# This is the logistic equation which I use in order to adjust my stochastic rounding to up and down weight the probability of 
# rounding my fractions as a value of 1.  I have presets based on experience working with my own simulation system.  I've not implemented
# a means of adjusting these within a run, but that could be done by simply updating the arguments passed through function calls.
### NOTE: this is not currently in use.... It's been found to generate worse fit to theoretical results.
logisticRounding <- function(func_decimal,func_midPoint = 0.5, func_slopeAngle = 12, func_maxReturn = 1){
	return( func_maxReturn/(1 + exp(-func_slopeAngle * (func_decimal - func_midPoint))) )	
}


# this is just a copy of the logistic growth general form equation but expanded to allow users to set the midpoint
#f(x) = K / (1 + ((K - N_0)/N_0) *exp-k(x-x_0))  ; Where x_0 is an adjustment to the position of the midpoint of the curve's maximum value
# K = the curves maximum value, k = the steepness of the curve (growth rate), and N_0 is the starting population
# This can be used to generate a curve of the growth between time points and we can use the normalised
# proportion of that curve in each step to "allocate" the amount of birth between points
# This script is only setup to give a final population based on the start and carrying capacity, again we can change the basal exponent at users risk.
# If we don't enforce a particular midpoint, nor different basal exponent, then this is the logistic growth equation
### NOTE: The midAdjust value allows the user to shift the point around which exponential growth in centered by setting the generation value at which
###			the population should reach the starting population value
logisticGrowth <- function(func_rate, func_step, func_startPop = NULL, func_maxPop = NULL, func_midAdjust = 0, func_basalExponent = exp(1)){ 
	return( func_maxPop / (1+ (((func_maxPop - func_startPop)/(func_startPop)) * func_basalExponent^(-func_rate * (func_step - func_midAdjust))) ) )
}

# This function is the logistic map, a discrete form of the logistic equation that will allow me to calculate growth
# of the population through steps, but there is no obvious means of changing the mid-point of inflection as per the 
# differential of logisticGrowth.  This map is given by:  N_t+1 = N_t + r * (N_t (K - N_t)/K)
logisticMap <- function(func_rate, func_startPop, func_maxPop){
	func_startPop + func_rate * (func_startPop * ((func_maxPop-func_startPop)/func_maxPop))
}

# This function will allow the starting or ending population size to be determined 
# based on a given rate per step, step number, either  the start or stop population size,
# and with the option to set the basal value to exponentiate (it defaults to e making this an exponential growth function) 
expGrowth <- function(func_rate, func_step,func_startPop = NULL, func_endPop = NULL, func_basalExponent = exp(1)){
	# If the user has defined a starting population value we assume we're calculating the final value
	if(is.null(func_endPop) && !is.null(func_startPop)){
		func_tmpRate <- func_basalExponent ^ (func_rate * func_step)
		if(any(func_tmpRate < 1)){ func_tmpRate[which(func_tmpRate < 1)] <- 1 }
		return( func_startPop * func_tmpRate )
	
	# If the user has entered an end value we assume we're calculating the starting value
	} else if(!is.null(func_endPop) && is.null(func_startPop)){
		return( log(func_endPop,base = func_basalExponent)/(func_rate * func_step) )
	} else {
	# otherwise we complain that this function is not meant to perform calculation of other terms.
		stop(paste("There was a problem when calculating the growth as starting pop was ",func_startPop," and final pop was ",func_endPop," please review",sep=""))
	}
}

# This function takes the simulation run's parameters for type of growth, as well as the population's data, and returns the amount of growth
# for each lineage.  As values of growth calculated may include decimals, and this simulation tool deals in whole individuals, we round
# the values in a stochastic fashion with probability == size of the remainder.  Rounding is done in this way as to not prejudice small populations.
# This mehtod is prepared to handle values such as exponential, logistic, or constant where logistic is the default.
# the second line of variables are optinal and depend on which growth type is being called.  The func_deathScale argument is to control if 
# we scale the number of births by the number of deaths, if this is done then the growth between time steps should match up to the analytical solutions
# of the respective growth forms; NOTE: That constant is not affected since it already considers deaths, nor is Poisson since it's growth form
# is simply a statement that the number of offspring is poisson distributed, and is not in fact explicitly a function about growth.
# The func_drift argument is a logical toggle for selecting if births are deterministic or simply an expectation that will be passed through a Poisson call.
# Last part is the control of rounding values, this affects whether or not values will be integer or not.
growthFunction <- function(func_inSize, func_inFitness, func_bProb, func_sizeStep, func_growthForm = c("logistic","exponential","constant","poisson"),  
						func_deaths = NULL, func_carryingCapacity = NULL, func_basalRate = NULL, func_deathScale = FALSE, func_drift = TRUE, 
						func_roundValues = TRUE){
	# We start by triming down to the first growth form values passed, or this means that it will default to logistic
	# then we check that what is being used is something that can be recognised, else we complain.
	func_growthForm <- func_growthForm[1]
	if(!is.element(func_growthForm[1],c("logistic","exponential","constant","poisson"))){
		stop(paste("Incorrect growth form was passed as being ", func_growthForm," please review",sep=""))
	}
	# Other sanity check, if the vector of populations is all zero, we throw this as well
	if(all(func_inSize == 0)){
		# No individual means no births.... end of story.
		return( func_inSize )
	} else if( any(is.element(func_inSize,c("NaN","-Inf","Inf"))) ){
		# We complain about this, though it shouldn't be a problem...
		stop("One of the population elements had NaN, or Inf size, please review")
	}
	
	# Great now we handle growth differently based on which form was called
	func_numBirths <- if(func_growthForm == "exponential"){
							# this is simply exponential growth for our lineage(s) which can be computed by the increase in individuals
							# And we remove the original amount of each population since this value will be added to the existing numbers.
							tmpBirths <- expGrowth(func_rate= func_inFitness * func_bProb *
													log(if(!is.null(func_basalRate)){func_basalRate}else{2}),
													func_step= func_sizeStep,
													func_startPop= func_inSize) - 
										(func_inSize - if(func_deathScale){func_deaths}else{0})
							if(func_drift){
								tmpUpdate <- which(tmpBirths %%1 != 0)
								tmpBirths[tmpUpdate] <- addDrift(tmpBirths[tmpUpdate], func_integerValues = func_roundValues) 
							}
							tmpBirths
						} else if (func_growthForm == "logistic"){
							## Now we calculate the growth rate, per lineage, which is modified by the density dependent logisticMap term
							func_densityTerm <- (logisticMap(func_rate= func_inFitness * func_bProb * func_sizeStep * 
																			if(!is.null(func_basalRate)){func_basalRate}else{2},
																func_startPop= sum(func_inSize),
																func_maxPop= func_carryingCapacity)/sum(func_inSize)) - 1	
							# Now we computed the actual number of newborn, if asked to compensate for the dead then we add back in a number of newborn
							# given the number of dead as determined by a call to growthFunction with constant growth form
							tmpBirths <- (func_inSize * func_densityTerm) + 
											if(func_deathScale && !all(func_deaths == 0)){
												growthFunction(func_inSize = func_inSize,
																func_inFitness = func_inFitness, 
																func_bProb = func_bProb, 
																func_sizeStep = func_sizeStep, 
																func_growthForm = "constant", 
																func_deaths = func_deaths, 
																func_carryingCapacity = func_carryingCapacity, 
																func_basalRate = func_basalRate, 
																func_deathScale = FALSE, 
																func_drift = FALSE,
																func_roundValues = func_roundValues)
											}else{
												0
											}
							if(func_drift){
								tmpUpdate <- which(tmpBirths %%1 != 0)
								tmpBirths[tmpUpdate] <- addDrift(tmpBirths[tmpUpdate], func_integerValues = func_roundValues) 
							}
							tmpBirths										
						} else if (func_growthForm == "poisson"){
							# This is Poisson growth, where an individual has a number of offspring drawn from the poisson distribution
							# this is meant to keep a population approximately equal so long as prob_b == prob_d. 
							unlist(lapply(func_inSize * func_inFitness * func_bProb * func_sizeStep,function(thisLineage){ rpois(1,thisLineage) }))
						} else if(func_growthForm == "constant"){
							# To consider constant growth, I calculate the deterministic growth, then if there is drift this value gets
							# passed to the drift function with basis on rounding.  These values are then scaled back to the number of deaths
							tmpBirths <- unlist(lapply(func_inSize * calc_relativeFitness(func_fitVector = func_inFitness) * func_bProb * func_sizeStep,
													function(thisLineage){ 
														if(func_drift){
															addDrift(thisLineage, func_integerValues = func_roundValues)
														} else {
															thisLineage
														} }))
							# If there were no births we can simply return that vector of zeroes...
							if(sum(tmpBirths) != 0){
								tmpBirths <- tmpBirths/sum(tmpBirths) * sum(func_deaths)
							}
							tmpBirths

						}
						
	# We round our population values stochastically, if we're told to round values to integers, and if the value is not already an integer 
	# - found to save time, and becuase above when we add drift, the values may already be integers if we're roundign values (due to Poisson stochasticity)
	if(func_roundValues){
		tmpUpdate <- which(func_numBirths %%1 != 0)
		if(length(tmpUpdate) > 0){
			func_numBirths[tmpUpdate] <- unlist(lapply(func_numBirths[tmpUpdate], function(thisValue){ 
						# We handle that the number of births may be a negative value 
						if(sign(thisValue) == 1){
							floor(thisValue) + rbinom(1,1, thisValue %%1)
						} else {
							ceiling(thisValue) - rbinom(1,1, thisValue %%1)
						} 
				}))
		}
	}
	# Now we adjust the number of births calculated, this is likely not actually required except that there is some stochastic rounding
	# which occurs meaning that we may overshoot population size.  
	### NOTE: this obviously only matters for constant growth  
	if(func_growthForm == "logistic"){
		# We only need to adjust the size if the proposed growth is not ba;anced AND there is at least some proposed growth 
		if((sum(func_inSize, func_numBirths) - if(func_deathScale){sum(func_deaths)}else{0}) > func_carryingCapacity && !all(func_numBirths == 0)){
			func_numBirths <- adjustBirths(func_adjVector= func_numBirths,
											func_sumTotal= func_carryingCapacity - sum(func_inSize) + if(func_deathScale){sum(func_deaths)}else{0},
											func_roundValues = func_roundValues)
		}
	} else if(func_growthForm == "constant"){
		if(sum(func_numBirths) != sum(func_deaths) && !all(func_numBirths == 0)){
			func_numBirths <- adjustBirths(func_adjVector= func_numBirths,
											func_sumTotal= sum(func_deaths),
											func_roundValues = func_roundValues)
		}
	}
	# We now return the new lineage sizes given the births
	return( func_numBirths )
}

# This function is used to take a vector of values, and a desired total sum, and ensure the vector sums to this.
# This is basically a means to rescale our proposed number of births to the amount permissible.
adjustBirths <- function(func_adjVector, func_sumTotal, func_roundValues = func_roundValues){
	# This is a likely redundant validation as we ought not to be in this function unless this is true, but....
	if(sum(func_adjVector) != func_sumTotal){
		func_adjVector <- func_adjVector/sum(func_adjVector) * sum(func_sumTotal)
	}
	# If we're meant to round our values then we do so
	if(func_roundValues && any(func_adjVector %% 1 != 0)){
		func_tmpUpdates <- which(func_adjVector %% 1 != 0)
		func_adjVector[func_tmpUpdates] <- unlist(lapply(func_adjVector[func_tmpUpdates],function(thisPop){ 
																floor(thisPop) + rbinom(1, size = 1, prob = thisPop %% 1) 
														}))
	}
	# Due to rounding the value may anew not sum properly to the permissible value, but the problem ought not to be possible unless rounding is considered!
	# Hence our resolution to this using integer adjustments, but as a safety we won't bother if the difference is less than 1
	if(abs(diff(c(sum(func_adjVector), func_sumTotal))) >= 1){
		# The adjustment will be the difference between our permissible and proposed values
		func_adjustments <- func_sumTotal - sum(func_adjVector)
		# we distribute our adjustment randomly among the vector indexes with probability proportional to their relative to the total
		# Also, because this is an adjustment to the number of births, and we're doing integer adjustments, we won't adjust any values that are less than 1 
		func_tmpAdjustable <- which(abs(func_adjVector) >= 1)
		func_proposedAdjustments <- NULL
		func_tmpCounter <- 0
		# Now, for as long as the porposed number of adjustments does not equal the change required, we recalculate, unless we've tried this 100 times in
		# which case we complain and break.
		while(sum(func_proposedAdjustments) != abs(func_adjustments) && func_adjustments >= 1 && func_tmpCounter < 100){
			# We propose adjustments 	
			func_proposedAdjustments <- table(sample(func_tmpAdjustable, 
														abs(func_adjustments),
														replace=TRUE, 
														prob = abs(func_adjVector[func_tmpAdjustable]/sum(func_adjVector[func_tmpAdjustable]))))
			# We reduce adjustments to be be no greater than the value passed 
			func_proposedAdjustments <- sapply(names(func_proposedAdjustments) ,function(thisNamed){
												return( min(func_adjVector[as.numeric(thisNamed)], func_proposedAdjustments[thisNamed]) )			
											})
			# This is an object format work-around for times when there are no proposed adjustments
			if(length(func_proposedAdjustments) == 0){ func_proposedAdjustments <- NULL }
				
			# We update our counter
			func_tmpCounter <- func_tmpCounter + 1
		}
		# We report if the counter was reached to break this chain
		if(func_tmpCounter >= 100){
			stop("While adjustingBirths and trying to aportion the adjustments got stuck in a loop for 100 replicates, please review...")	
		}
		if(!is.null(func_proposedAdjustments)){
			# Otherwise we apply this proposed adjustment by either adding or subtracting the sampled numbers based on the sign 
			func_adjVector[as.numeric(names(func_proposedAdjustments))] <- func_adjVector[as.numeric(names(func_proposedAdjustments))] + 
																			if(func_adjustments > 0){ 
																				as.vector(func_proposedAdjustments)
																			} else {
																				- as.vector(func_proposedAdjustments)
																			}
		}
	}
	# We now return our vector which should now be adjusted.	
	return( func_adjVector )
}

# Now we need a function to generate the fitness values for a vector of genotypes that are passed along.
# We also provide a means to send forward the fitness of focal genotype(s) which may be the immedtae ancestor (Additive requires), 
# an optimal genotype (RMF required), the independent fitness of all sites (NK required)
# (some models care about this, e.g. Fisher's Geometric....)
### NOTE: As we do not define all genotypes in the space (due to the genotype space of large genomes being impractical)
###			we can only reasonably implement landscape models which are not based on the fitness of anything other than the local space
fitnessLandscape <- function(tmpGenotypes, tmp_focalFitness, landscapeModel = "HoC", tmp_ancestralFitness = const_ancestFitness,
									tmp_weightsRMF = const_RMF_indWeight, tmp_optimaRMF = const_RMF_globalOptima, 
									tmp_correlationsNK = const_NK_interactionMat, tmp_NK_ancestDep = const_DepbySite_ancestFitness, 
									relativeFitness = TRUE, func_genomeLength = genomeLength, 
									func_distribution = constDist, func_distParameters = const_distParameters,
									func_distAsS = const_distAsS, func_sepString = sepString){
											#startTime <- proc.time() # This is for reporting on run times but is not needed....
	# We create a list of the mutations in each of the genotypes having been passed, only if not the "Fixed" model...
	if(landscapeModel != "Fixed"){
		tmp_genotypeList <- lapply(strsplit(tmpGenotypes,func_sepString),as.numeric)
	}
	# This is a function of small convenience to convert the state of the tmp_optimaRMF object
	if(!is.na(tmp_optimaRMF) && is.character(tmp_optimaRMF)){ tmp_optimaRMF <- as.numeric(strsplit(tmp_optimaRMF,func_sepString)[[1]]) }
	
	# This creates the overall fitness associated to mutations 
	fitnessVec <- if(landscapeModel == "HoC"){
						# For the HoC since each site is genotype is uncorrelated we will simply generate a random number for each.
						tmpReturn <- fitnessDist(length(tmp_genotypeList), tmpDistribution = func_distribution, tmpParameters = func_distParameters)
						# If the distribution draws relative fitnesses but needs to be used as selection coefficients...
						if(func_distAsS){
							tmpReturn <- tmpReturn - 1
						}
						tmpReturn
					# Now we will deal with the instance of a simple additive model
					} else if(landscapeModel == "Additive") {
						# In the additive model there is no interaction, so we find the number and position of mutations in the genotype(s) - in tmp_genotypeList - 
						# then review the constant independent fitness per site matrix which should be passed similar to the NK model.
						as.vector(sapply(tmp_genotypeList,function(thisGenotype){
									# So each genotype is already split into it's constituent mutations, so we simply add the values using the 
									# siteBystate matrix which should be represented in tmp_focalFitness
									sum(c(tmp_ancestralFitness,  tmp_focalFitness[setdiff(thisGenotype,1:nrow(tmp_focalFitness)),"0"], tmp_focalFitness[thisGenotype,"1"]))
								}))					
					} else if(landscapeModel == "NK") {
						###### This is as per the model information in Kauffman 1989
						# In this case we simply need to know what is the independent fitness of all sites, in each state, thus tmp_focalFitness should be a
						# matrix with as many columns as there are states each site can assume, as many rows as genome length, and all states possible are included.
						if(!is.matrix(tmp_focalFitness) || ncol(tmp_focalFitness) != length(const_siteStates) || nrow(tmp_focalFitness) != func_genomeLength || !all(is.element(const_siteStates,colnames(tmp_focalFitness))) ){ 
							stop("The number of states, and/or the shape of the tmp_focalFitness passed to a call for fitnessLandscape function, is not correct please review")
						}
						# AND we need to know which others sites each focal site depends upon.  BUT if the interactions == 0 then we can't ahve been passed a matrix
						# and thus tmp_correlationsNK should remain as a NULL value which can be handled downstream
						if(const_numInteractions != 0){
							if(!is.matrix(tmp_correlationsNK) || ncol(tmp_correlationsNK) != const_numInteractions || nrow(tmp_correlationsNK) != func_genomeLength){
								stop("The shape of the tmp_correlationsNK passed to a call for fitnessLandscape function, is not a matrix and or does not have as many columns as suggested by const_numInteractions, please review")
							}
						} else if(const_numInteractions == 0 && !is.null(tmp_correlationsNK)){
							# If there are no interactions then the object should be null
							stop("The tmp_correlationsNK object should have been passed as null since there are no reported NK correlations, please review")
						}
						
						# Ok, then we can calculate what the fitness is for a genotype by working out mean the fitness effect of each site (w_i).
						# The fitness value of each site is a function of the state of the site and all other sites upon which it depends.
						# We're defining this function as simply the mean of the independent fitness values
						## So for each genotype, we pass the states of the sites, which in the initial implementation was recorded by the ID string being the mutant sites
						#### NOTE: As a computational shortcut, I've built the "const_DepbySite_ancestFitness" object so that we know the bySite dependent fitness
						####	for sites which do not hiold mutations, and thus we update the per site dependent fitness only for those sites which have dependence
						####	on a mutant site.	 
						as.vector(sapply(tmp_genotypeList,function(thisGenotype){
									# This finds the independent siteWise fitness value of each site, based on their states as per thisGenotype
									tmp_indFitness <- tmp_focalFitness[,"0"]
									tmp_indFitness[thisGenotype] <- tmp_focalFitness[thisGenotype,"1"]
									# We simply return the mean of the ancestral dependent fitness values, having updated those sites which
									# are mutants.
									tmp_depFitness <- tmp_NK_ancestDep
									tmp_mutCorrelates <- unique(which(apply(tmp_correlationsNK,MARGIN=1,function(theseSites){ any(is.element(thisGenotype,theseSites)) })))
									for(thisSite in union(thisGenotype,tmp_mutCorrelates)){
										tmp_depFitness[thisSite] <- mean(tmp_indFitness[c(thisSite, tmp_correlationsNK[thisSite,])])
									}
									return( mean(tmp_depFitness) )
						}))
					} else if(landscapeModel == "RMF"){
						# in this case we need to know the optimal type (the string passed by tmp_focalFitness) and the weightings of independent
						# and interaction components of the genotype (assumed to be the first and second component of the tmp_wightsRMF  length 2 vector)
						if(length(tmp_weightsRMF) != 1){ stop("The vector of weights passed to RMF fitnessLandscape model is not length 1, please review.") }
						# From Neidhart 2014, we find that the fitness value for a genotype, in their modified RMF landscape model, is given by:
						# F(genotype) = -cD(genotype,optima_genotype) + random_component ; where c is the scaling of the independent fitness contribution
						# and the D() function is a measure of the Hamming distance, between a genotype with the independent fitness optima and the focal genotype.
						# The random component is drawn from a distribution and each genotype has it's own value, meaninging in a binary state genotype there are 2^(func_genomeLength) elements
						# hence we'll be drawing only when we create new genotype's and assign them fitness.
						tmp_randomTerm <- fitnessDist(length(tmp_genotypeList), tmpDistribution = func_distribution, tmpParameters = func_distParameters)
						# If the distribution draws relative fitnesses but needs to be used as selection coefficients...
						if(func_distAsS){
							tmp_randomTerm <- tmp_randomTerm - 1
						}
						## The math is just the independent term which needs to know the distance to optima, and then a random component
						## The random component should be of the order of selection coefficients 
						as.vector( -tmp_weightsRMF * 
									unlist(lapply(tmp_genotypeList,function(thisGenotype){ 
										# The distance between two genotypes can be calculated as the length of one - the number of sites 
										# the second shares in common, + the number of sites the second has in excess
										tmpOverlap <-is.element(thisGenotype, tmp_optimaRMF)
										length(tmp_optimaRMF) - length(which(tmpOverlap)) + length(which(! tmpOverlap))
									})) + 
									tmp_randomTerm )
					
					}	else if (landscapeModel == "Fixed"){
						# This assumes the user has passed a named vector as the tmp_focalFitness with a complete description 
						# of the fitness value for each genotype that could arise in simulations.  The names should be "binaryString" and values the "fitness".
						# The binary string values can be used to return the fitness values
						if(all(is.element(tmpGenotypes,names(tmp_focalFitness)))){
							tmp_focalFitness[tmpGenotypes]
						} else {
							stop("There was a problem trying to use the Fixed fitness landscape model, review input make certain all possible genotypes are declared.")
							Sys.sleep(20)
							q(save="no")
						}	
					}# This is the fitnessVec creation closing
	
	# Now, if we're using relativeFitness, then we calculate the relative fitness value using our function, which already handles instances
	# of the ancestral fitness being passed as a negative value or of zero or NULL.
	if(relativeFitness){ fitnessVec <- calc_relativeFitness(fitnessVec, func_ancestFit = tmp_ancestralFitness) }	
	# We round the fitness vector off to the 4th decimal place ought of simple convenience and aesthetics.	
	return( round(fitnessVec,4) )
}

# This is a function to draw fitness values from a distribution,
### NOTE: The user must pass at least the minimum correct number of parameters for the distribution chosen
###			and they must be in the same order as the R function expects them 
fitnessDist <- function(tmpDraws, tmpDistribution = constDist, tmpParameters = const_distParameters ){
	if(tmpDistribution == "Fixed"){
		# It has been found that if the number of values from which to be sampled is == 1 (meaning we passed two values), then it samples integers...
		# so we force a value to be replciated twice if this is the case
		return( sample(if(length(tmpParameters) == 2){ rep(tmpParameters[2],2) } else { tmpParameters[2: length(tmpParameters)]}, 
						tmpDraws, replace = as.logical(tmpParameters[1])) )
	}else if(tmpDistribution == "Gamma"){
		return( rgamma(tmpDraws, shape = tmpParameters[1], rate = tmpParameters[2]) )
	} else if(tmpDistribution == "Uniform") {
		return( runif(tmpDraws, min = tmpParameters[1], max = tmpParameters[2]) )
	} else if(tmpDistribution == "Normal") {
		return( rnorm(tmpDraws, mean = tmpParameters[1], sd = tmpParameters[2]) )
	} else if(tmpDistribution == "Chi2") {
		return( rchisq(tmpDraws, df = tmpParameters[1], ncp = tmpParameters[2]) )
	} else if(tmpDistribution == "beta") {
		return( rbeta(tmpDraws, shape1 = tmpParameters[1] , shape2 = tmpParameters[2]) )
	}  else if(tmpDistribution == "exp") {
		return( rexp(tmpDraws, rate = tmpParameters[1]) )
	}  else if(tmpDistribution == "evd") {
		return( rgev(tmpDraws, loc = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
	}  else if(tmpDistribution == "rweibull") {
		return( rrweibull(tmpDraws, loc = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
	}  else if(tmpDistribution == "frechet") {
		return( rfrechet(tmpDraws, location = tmpParameters[1], scale = tmpParameters[2], shape = tmpParameters[3]) )
	} else if(tmpDistribution == "skewNorm") {
		# We set the location parameter to 1.05 so the mean is near 1 defined by the omega and alpha
		tmpReturn <- rsn(tmpDraws, xi = tmpParameters[1], omega = tmpParameters[2], alpha = tmpParameters[3], tau = tmpParameters[4])
		# This distribution allows for negative value space, thus we adjust....
		tmpReturn[which(tmpReturn < 0)] <- 0
		return( as.vector(tmpReturn) )
	}
}


# This function will take a focal genotype and then create all the unique genotypes which have one mutation more or less
# which are already not within the reference list.
# NOTE: I'm having this function inherit some pre-definitions, these must be respected!
createGenotypes <- function(tmp_focalGenotype, tmp_focalFitness, maxHamming, tmp_landModel = "HoC", tmp_sepString = sepString,
								tmpDirection = allow_backMutations,  tmp_relativeFitness = const_relativeFitness, tmp_currNeighbours = NULL, 
								tmp_genCon = connections_dataBase$genotypeSpace, tmp_tableSplit = db_splitTables, tmp_genomeLength = genomeLength,
								tmp_distAsS = const_distAsS, ...){
	# If we haven't previously defined what are the possible neighbours for tmp_focalGenotype, we do so here
	if(is.null(tmp_currNeighbours)){ tmp_currNeighbours <- comp_defineNeighbours(func_tmpGenotype = tmp_focalGenotype, func_tmpDirection = tmpDirection)}
	# We can define aspects of our focalGenotype from the binaryString which ought to have been passed 
	tmp_numMuts <- length(strsplit(tmp_focalGenotype, tmp_sepString)[[1]])
	# We find the range of mutants that could be reached in the next step, and bound this by the no-mutant type and the genomeLength
	tmp_mutsRange <-  max(0,(tmp_numMuts - if(tmpDirection){ maxHamming }else{ 0 })): min(tmp_genomeLength ,(tmp_numMuts + (maxHamming)))
	# Unless we allow back mutations AND the maxHamming is > 1, we don't consider the current numMuts as part of the neighbourhood.
	if(!tmpDirection && maxHamming <= 1){
		tmp_mutsRange <- tmp_mutsRange[-which(tmp_mutsRange == tmp_numMuts)]
	}
	
	# We'll need to query what tables are within the database, both within this function and the tmp_found_neededNeighbours call
	tmp_dbTables <- dbListTables(tmp_genCon)
														#startTime <- proc.time() # This is for reporting on run times but is not needed....
	# We now search for which neighbours already exist within our SQL database
	tmp_found_neededNeighbours <- comp_find_neededNeighbours(tmp_possibleNeighbours = tmp_currNeighbours, 
															tmp_focal_numMuts = tmp_numMuts,
															tmpRange_numMuts = tmp_mutsRange, 
															tmp_refTables = tmp_dbTables,
															tmp_genCon = tmp_genCon)
	# Now so long as there are neighbours that need to be created, we do so.
	if(length(tmp_found_neededNeighbours) != 0){
		# Now we assign each of these neighbours to a numMuts category for table updating reasons
		# This is done by splitting the strings of all the new neighbours to be assigned, and finding how many mutations based 
		# on the length of mutation positions recorded in the compressed nomenclature
		neighbourMuts <- lapply(tmp_mutsRange,function(x){ 
							return( tmp_found_neededNeighbours[which(sapply(strsplit(tmp_found_neededNeighbours, sepString),function(y){ length(y) }) == x)] ) })
		names(neighbourMuts) <- sapply(tmp_mutsRange,nameTable, "func_subNaming"=tmp_tableSplit)
		# Now the names of neighbourMuts will be the stringConstant table name for names in our database which reference to a particular numMuts value
		# So we look if there are any tables which contain these names, if not we create a table, if so see the  <else>  section.
		# Also, if there any neighbourMuts list elements without length... meaning there are not mutants needed to be created, we don't pass that instance. 
		for(thisTable in names(neighbourMuts)[which(lapply(neighbourMuts,length) > 0)]){	
			# We'll be assigning new fitnessLandscape space defined by thisTable's neighbours into some table..
														#startTime <- proc.time() # This is for reporting on run times but is not needed....
			### NOTE: This call would have a problem generating our ancestral genotype while it's binary string is considered as "",
			###			I consider this a trivial problem as I expect to always have started runs by defining the no-mutant genotype.
			assign(tmp_genotypeObject,create_genotypeFrame(tmpID = tmp_newID: (tmp_newID + length(neighbourMuts[[thisTable]]) - 1),
														tmpStrings = neighbourMuts[[thisTable]],
														tmpFitnesses = comp_fitnessLandscape(tmpGenotypes = neighbourMuts[[thisTable]], 
																							tmp_focalFitness = tmp_focalFitness, 
																							landscapeModel = tmp_landModel, 
																							relativeFitness = tmp_relativeFitness,
																							func_distAsS = tmp_distAsS)), 
					pos=".GlobalEnv")
			# We now update the value of the tmp_newID value for the next iteration.	
			tmp_newID <<- tmp_newID + length(neighbourMuts[[thisTable]])
			# This checks if we need to make a new table or are adding to existing tables, it's controlled by an external logical 
			# Of whether or not we have a max size of DB table sizes
			if(!any(grepl(thisTable, tmp_dbTables)) || !tmp_tableSplit){
				# In this case we are simply taking the genotypeFrame for all neighbours in thisTable and copying it to the database
				dbWriteTable(tmp_genCon,
							name = ifelse(tmp_tableSplit, paste(thisTable,1,sep=""), thisTable), 
							value = eval(as.name(tmp_genotypeObject)), 
							append=TRUE)
			} else {
				# This is trickier, this means we need to identify the number of tables which share this string, and for the last created one
				# we querry if it has the maximum number of rows, if so we create a new one, if not we add to it with these
				tmp_similarTables <- tmp_dbTables[which(grepl(thisTable, tmp_dbTables))]
				# Now to find which is the oldest we look for the highest value in the 3rd separated position
				tmpOldest <- tmp_similarTables[which.max(sapply(tmp_similarTables,function(x){ as.numeric(strsplit(x,tmp_sepString)[[1]][3]) }))]
				# We querry the number of rows in this table and if there are fewer than our maxRows we'll add, otherwise we create a new one
				if(nrow(dbGetQuery(tmp_genCon,paste("SELECT genotypeID FROM ", tmpOldest,sep=""))) < maxRows){
					# This means we can simply add our new table to the existing one
					dbWriteTable(tmp_genCon,
									name = tmpOldest,
									value = eval(as.name(tmp_genotypeObject)), 
									append=TRUE)
				} else {
					# We write a new table which simply takes the index (third indexed sepString position) and advance to the next one....
					dbWriteTable(tmp_genCon,
									name = paste(paste(strsplit(tmpOldest,tmp_sepString)[[1]][-3],collapse=tmp_sepString),as.numeric(strsplit(tmpOldest,tmp_sepString)[[1]][3])+1,sep=tmp_sepString),
									value = eval(as.name(tmp_genotypeObject)))
				}
				
			}# This closes out  if we needed to make a new table or add to an existing one	
		} # This closes out going through all the table types that might need to be created as a result of numMuts in genotypes
		
		# Now we go to the ancestral genotype's table and we update that it has been explored
		# This requires we find it's table, so we search for all tables in it's mut range
		tmp_neededTables <- tmp_dbTables[which(grepl(nameTable(tmp_numMuts, func_subNaming=tmp_tableSplit), tmp_dbTables))]
		# If the focalGenotype is anything other than our absolute ancestor we need to find the tables and use it's binaryString
		# However, for the ancestor we don't do this since it has no binaryString value.... we'd have needed it's genotypeID
		tmp_ancestFind <- 0
		tmp_ancestTable <- tmp_dbTables[which(grepl(nameTable(0, func_subNaming=tmp_tableSplit), tmp_dbTables))]
		names(tmp_ancestFind) <- tmp_ancestTable
		if(tmp_focalGenotype != ""){
			# We find the proper table and then the genotypeID of our ancestor 
			tmp_ancestFind <- sapply(tmp_neededTables, function(thisTable){ 
										as.vector(unlist(dbGetQuery(tmp_genCon, paste("SELECT genotypeID FROM ",thisTable,' WHERE binaryString = "', tmp_focalGenotype,'"',sep=""))))
								},simplify=FALSE)
			# There should be a single returned non NULL value
			tmp_ancestTable <- names(tmp_ancestFind)[which(sapply(tmp_ancestFind,length) == 1)]
		}
		dbExecute(tmp_genCon, 
					paste("UPDATE ", tmp_ancestTable ,' SET isExplored = 1 WHERE genotypeID = ', tmp_ancestFind[tmp_ancestTable],sep=""),
					synchronous = NULL)
		
	} # This closes out if we had any neighbours that needed to be created.
	
	# We don't return anything, we've performed work....
	return( NULL )
}


# This function will help us investigate the local neighbourhood of a focal genotype and see if we need to create new genotype information
# NOTE: It assumes the tmp_focalGenotype being passed to it is either a character string, or a numeric vector
find_neededNeighbours <- function(tmp_possibleNeighbours, tmp_focal_numMuts, tmp_refTables, maxHamming = max_numMutations, 
							tmp_tableSplit = db_splitTables, tmp_genomeLength = genomeLength,
							tmpDirection = allow_backMutations, tmpRange_numMuts = NULL, tmp_genCon = connections_dataBase$genotypeSpace){
													#startTime <- proc.time() # This is for reporting on run times but is not needed....
	# If we haven't been passed the vector of all tables in our reference database, then we define it here
	if(is.null(tmp_refTables)){ tmp_refTables <- dbListTables(tmp_genCon) } 
	# Now we find which of these are not elements of binaryStrings already within our database, we look for possible neighbours that are not 
	# found within the database (called using the select function of "dplyr"
	# This is started by querrying for binaryStrings from among meaningfull tables in our database , we define which meaningful tables exist:
	# If we've passed a range for the number of mutants then we skip this, otherwise we redefine the NULL parameter
	tmp_neededTables <- NULL
	if(is.null(tmpRange_numMuts)){
		tmpRange_numMuts <- max(0,(tmp_focal_numMuts - if(tmpDirection){ maxHamming }else{ 0 })): min(tmp_genomeLength,(tmp_focal_numMuts + maxHamming))
		tmp_neededTables <- tmp_refTables[unlist(lapply(tmpRange_numMuts[-which(tmpRange_numMuts == tmp_focal_numMuts)],function(x){ which(grepl(nameTable(x, func_subNaming=tmp_tableSplit), tmp_refTables)) }))]
	} else {
		tmp_neededTables <- tmp_refTables[unlist(lapply(tmpRange_numMuts,function(x){ which(grepl(nameTable(x, func_subNaming=tmp_tableSplit), tmp_refTables)) }))]
	}
	# We now ask if there are any table from which we can compare for Neighbours, if there aren't (length == 0) we needn't do anything as all neighbours need to be made.
	if(length(tmp_neededTables) == 0){
		# Well in this case we can return all possible neighbours as ones that need to be created.
		return( unlist(tmp_possibleNeighbours) )
	} else {
		# Then we querry our database, among the meaningfull tables, and take the possible neighbours that are not in the database yet.
		# We grab as a vector all the binary strings which are within those meaningfull table
		tmp_binaryStrings <- gsub("[[:space:]]","",paste("\'", tmp_possibleNeighbours,"\'",collapse=','))

		# We now look for which of these possible neighbours exist within our database
		return ( setdiff(tmp_possibleNeighbours, 
						unlist(dbGetQuery(tmp_genCon, 
											paste('SELECT binaryString FROM ',
													tmp_neededTables, 
													' WHERE binaryString IN (',
													tmp_binaryStrings,')',
													collapse=" UNION ")) 
						)) )
	}
}

# This is a function which informs us of which binary strings are the nearest neighbours of a focal genotype
# We expect that at least a binaryString single length character vector is passed, from this we can work out number of mutations,
# Otherwise if a vector is passed and/or number of mutations therein we use that directly, 
# lastly we need to be informed if we're allowing back mutations, technically the user can define the maxHamming distance, but we've only coded it for 1...
defineNeighbours <- function(func_tmpGenotype, func_tmpDirection = allow_backMutations, func_maxHamming = max_numMutations, 
								func_sepString = sepString, func_genomeLength = genomeLength){
	# We expect to be passed the func_tmpGenotype as the binaryString which is the position of each mutation with a separator
	func_tmpGenotype <- as.numeric(strsplit(func_tmpGenotype, func_sepString)[[1]])
	# Now, if there are no mutations in this genome, we can only add mutations, so we define each position in the genome as a possible mutant type
	return( if(length(func_tmpGenotype) == 0){
				as.character(seq(1, func_genomeLength))
			################################################################################################################################
			#### THESE TWO "else if" CALLS are hard coded for a maxHamming of 1, these will require meaningful updates if this changes. ####
			################################################################################################################################
			# Otherwise we need to consider adding mutations to each other position in the forward direction if one direction only 	
			} else if(length(func_tmpGenotype) > 0 && !func_tmpDirection){
				unlist(lapply(seq(1, func_genomeLength)[-func_tmpGenotype],function(thisPos){
						paste(c(func_tmpGenotype,thisPos)[order(c(func_tmpGenotype,thisPos))],collapse= func_sepString)
				}))
			# This is the case where we allow back mutations, meaning we permit the addition of mutations and their removal
			} else if(length(func_tmpGenotype) > 0 && func_tmpDirection){
				# So we'll return the concatenation of all those values where we add a mutation and remove one.
				c(unlist(lapply(seq(1, func_genomeLength)[-func_tmpGenotype],function(thisPos){
						paste(c(func_tmpGenotype,thisPos)[order(c(func_tmpGenotype,thisPos))],collapse= func_sepString)
					})),
					unlist(lapply(1:length(func_tmpGenotype),function(thisPos){
						paste(c(func_tmpGenotype[-thisPos])[order(c(func_tmpGenotype[-thisPos]))],collapse= func_sepString)
					})) )
			}) 
}

# This is a function for creating genotype_refDatabse index data.frames
create_genotypeFrame <- function(tmpID, tmpStrings,tmpFitnesses){
	if(length(tmpStrings) != length(tmpFitnesses)){
		stop("There was a problem trying to create a genotype data frame, as a result of the number of genotypes and fitnesses passed do not match")
	}
	return( data.frame("genotypeID"=tmpID,"binaryString"= tmpStrings, "fitness"= tmpFitnesses, "isExplored"= 0, stringsAsFactors=FALSE) )
}


# This is a function for the naming of tables to be copied to the SQL database.  NOTE: The func_tmpIndex value will be reduced to a single string element
# We also incldue a logical to suppress the use of sub-inxeding when creating names
nameTable <- function(func_tmpMutations, func_tmpIndex = NULL, func_baseString = string_tableNames, func_sepString = sepString, 
						func_splitName = FALSE, func_subNaming = db_splitTables){
	# As a afety I've included an option to split rather than build table names to recover the numMuts stored therein
	if(func_splitName){
			return( as.vector(unlist(lapply(func_tmpMutations, function(x) { strsplit(x, func_sepString)[[1]][2] }))) )
	}
	# Now if we're building names we need to know if we are indexing or not....
	if(func_subNaming){
		if(!is.null(func_tmpIndex)){
			return( paste(func_baseString,func_tmpMutations,paste(func_tmpIndex,collapse=func_sepString),sep=func_sepString) )
		} else {
			# We've added this trailing func_sepString element so that numeric indexces can later be better tracked within names using grepl calls with fixed = TRUE
			return( paste(func_baseString,func_tmpMutations,"",sep=func_sepString) )
		}
	# This ignores the presence of the func_tmpIndex value
	} else {
		return( paste(func_baseString,func_tmpMutations,sep=func_sepString) )
	}
}

# This is to be used for the names of our step report object tables
nameTable_step <- function(func_Index, funcSplit = FALSE){
	# This checks if we're asking for the name to be split or not
	if(funcSplit){
		# We return only the second piece of information given the setup of this function's naming practice.
		return( unlist(lapply(strsplit(func_Index,"_"),function(x){ x[2] })) )
	} else {
		return( paste("Step",func_Index,sep="_") )
	}
}

# This is to be used for the names of our neighbourhood object tables
nameTable_neighbourhood <- function(func_Index, funcSplit = FALSE){
	# This checks if we're asking for the name to be split or not
	if(funcSplit){
		# We return only the second piece of information given the setup of this function's naming practice.
		return( unlist(lapply(strsplit(func_Index,"_"),function(x){ x[2] })) )
	} else {
		return( paste("genotypeID",func_Index,sep="_") )
	}
}

# This is a little cheater function to close and re-open our database connection after each write
# We pass a filename for the connection we want.
resetDatabase <- function(func_conName, func_existingCon = NULL, func_type = "connect"){
	# Now we either open or close the connection for the conName passed
	if(func_type == "connect"){
		# If there is an existing connection we'll close it
		if(!is.null(func_existingCon)){ dbDisconnect(func_existingCon) }
		return( dbConnect(SQLite(), func_conName) )
	} else {
		return( dbDisconnect(SQLite(), func_conName) )
	}
}


# This is a function which will allow us to generate matrices easily for reporting on the population(s) in a step
reportPopulations <- function(func_numMuts, func_genotypeID, func_popSizes, func_fitnesses, func_births, func_deaths, func_mutants, func_progenitor){
	# We require that the user send us each piece of information being the same length, otherwise we complain.  Similarly if there are less
	# elements than epected (as per reportMat_colnames) we complain as well
	length_allInputs <- c(length(func_numMuts), length(func_genotypeID), length(func_popSizes), length(func_fitnesses),
						 length(func_births), length(func_deaths), length(func_mutants), length(func_progenitor))
	if( all(sapply(length_allInputs[-1],function(tmpLengths){ length_allInputs[1] == tmpLengths  })) ){
		if(length(length_allInputs) != length(reportMat_colnames)){
			stop("Trying to create population matrix with insufficient number of elements passed, please review")
			return( NULL )
		} else {
			# If so then we'll use these to build a matrix, otherwise we'll complain about the information passed along
			tmpReturn <- data.frame(func_numMuts, func_genotypeID, func_popSizes, func_fitnesses, func_births, func_deaths, func_mutants, func_progenitor,stringsAsFactors=FALSE)
			dimnames(tmpReturn) <- list(rep(func_genotypeID,max(length_allInputs)/length(func_genotypeID)), 
										reportMat_colnames)
			return( tmpReturn )
		}
	} else {
		stop("Trying to create population matrix with vectors of information that are not all the same size, please review")
		return( NULL )
	}
}

# This is a function to search our database and then find the binary string of the genotypeID passed
# We ask, for convenience, that the user help this search by defining the numMuts, this reduces
# the number of tables that we'll have to search, the two items passed ought to be single length vectors.
retrieve_binaryString <- function(func_genotypeID, func_numMuts = NULL, func_subNaming = db_splitTables,
									func_landscapeCon = connections_dataBase$genotypeSpace){
	# We find the needed tables by using the func_numMuts value passed
	tmp_refTables <- dbListTables(func_landscapeCon)
	# If we've been informed about the number of mutations this ID has, then we subset our tables
	if(!is.null(func_numMuts)){
		tmp_refTables <- tmp_refTables[which(grepl(nameTable(func_numMuts, func_subNaming = func_subNaming), tmp_refTables))]
	}
	# We now query all the table(s) for the information
	func_tmpReturn <- dbGetQuery(func_landscapeCon, paste("SELECT binaryString,isExplored FROM ", 
												tmp_refTables,
												' WHERE genotypeID IN(', 
												paste(func_genotypeID,collapse=","),
												')',
											collapse=" UNION "))
	# If nothing was found then we'll have a nrow == 0 object, in which case we simply search all tables in the DB
	# I don't use recursion to avoid being caught in an infinite loop, and instead put a crash out marker
	if(nrow(func_tmpReturn) == 0){
		func_tmpReturn <- dbGetQuery(func_landscapeCon, paste("SELECT binaryString,isExplored FROM ", 
												dbListTables(func_landscapeCon),
												' WHERE genotypeID IN(', 
												paste(func_genotypeID,collapse=","),
												')',
											collapse=" UNION "))
	}
	# Last check
	if(nrow(func_tmpReturn) == 0){
		stopError(paste("Problem finding the binaryString information for: ",func_genotypeID," please review",sep=""))
	} else {
		return( func_tmpReturn )
	}
}


# This function will compute the size of disturbance which occurs in a population, and thus the number of steps until its recovery
# The argument for maxSize exists but I am uncertain that the code's postAnalysis script will handle if this value changes.... it may, haven't thought about it very hard.
# The cappedGrowth logical is to define if the maxSize value is meant to be the carrying capacity or some initial value of the population.
compute_distGrowth <- function(func_distFactor, func_growthType, func_distType = const_distType, func_growthRate = const_growthRate,
								func_popSize = NULL, func_focalSize = NULL, func_manualGenerations = NULL, func_stepDivs = size_timeStep){
	# If the growth form is constant growth them we simply return that there will never be a disturbance
	# If the user has not set a number of generations prior to disturbance we set this as a large number.
	if(is.element(func_growthType, c("poisson","constant"))){
		return( c("factor"=0,
					"popLost"= 0,
					"stepReq"= if(is.null(func_manualGenerations )){
									1e12
								} else {
									ceiling(func_manualGenerations/size_timeStep)
								}) )
	}
	
	# We calculate the factor of population loss based on the growth type, for exponential or poisson this is simply a resseting of the population
	# to it's initial size
	func_tmpFactor <- as.vector(if(func_growthType == "logistic"){
									if(func_distType == "bottleneck"){
										func_distFactor["factor"]
									} else if(func_distType == "random") {
										rnorm(1,mean= func_distFactor["factor"],sd= func_distFactor["random"])
									} else {
										stop(paste('There was a problem with the definition of <func_distType> = ', func_distType, ' please review',sep=""))
									}
								} else if (func_growthType == "exponential") {
									(sum(func_popSize)/func_focalSize)
								})
	# Since, in theory, the dilution factor could be a value smaller than 1 (due to distribution size or growth form), this would mean that 
	# a proposed dilution would be of a form that increases the population size.  Thus if func_tmpFactor < 1 we set it to 1
	if(func_tmpFactor < 1){ func_tmpFactor <- 1 }
	
	# Now knowing the factor of loss, we calculate the number of steps in which growth will occur based on the division type
	# We also adjust this by the number of divisions with each step of growth.  We take the ceiling since we cannot run partial growth steps....
	func_growthSteps <- if(is.null(func_manualGenerations)){
							# If there aren't a mannual number of generations defined, then we take the log, base of growth rate,
							# of the disturbance factor.  However, if there was no disturbance, because perhaps this is an initialisation
							# step then we will set the this value to a default of 100 fold dilution to replicate serial passaging
							log(if(func_tmpFactor <= 0){ 100 } else { func_tmpFactor }, base = func_growthRate)
						} else {
							func_manualGenerations
						}
	func_growthSteps <- ceiling(func_growthSteps/func_stepDivs)
	# We return what we;ve calculated
	return( c("factor" = func_tmpFactor,
				"popLost" = 0,
				"stepReq" = as.vector(func_growthSteps)) )
}


# This is a function which calculates the amount of loss suffered by a population based on a dilution factor suffered by the population
# We return the lesser of the initial size OR the drawn value - as we cannot lose more than the initial number 
lossSampling <- function(func_inPopulation, func_dilutionFactor){
	# Here, the loss experienced is quite simply based on the the size of each lineage	
	return( apply(cbind(func_inPopulation, func_dilutionFactor),MARGIN=1,function(thisPop){ 
					min(thisPop[1],rpois(1, lambda = thisPop[1] * thisPop[2])) 
				}) )
}

# This function will allow me to dynamically build unique names for jobs performed with SHAPE, sepString is presumed to be in the evironment
# When splitting the jobID, setID, and repID values should be passed as logicals to control if we expect those to exist in what was passed.
name_batchString <- function(funcBase, func_setID = NULL, func_jobID = NULL, func_repID = NULL, 
								funcSplit = FALSE, func_sepString = sepString){
	# Only the base need be defined, if the others are not supplied then the base is assumed to be the name
	# The split option allows a user to extract information from the name if passed as the base
	if(!funcSplit){
		# I return a vector of name(s) built from the information passed, warning recycling may happen I don't control it.
		return( unname(apply(cbind(funcBase, func_setID, func_jobID, func_repID),MARGIN=1,function(thisInfo){ 
									paste(paste(thisInfo,collapse=func_sepString,sep=""),sep="")
							})) )
	} else {
		# We check that all the c(func_jobID, func_setID, func_repID) are logicals
		if(!all(is.logical(c(func_setID, func_jobID, func_repID)))){ 
			stop(print("Problem splitting the batchString of ",paste(funcBase,collapse=" ")," missing logicals.",sep=""))
		}
		# I return a matrix, possibly with only a single row, of the split information.
		func_tmpReturn <-  matrix(sapply(strsplit(funcBase,func_sepString),function(thisSplit){ 
							# The tmpReturn is built by building a vector based on the presence of ID logicals
							tmp_returnString <- thisSplit[1]
							thisSplit <- thisSplit[-1]
							for(func_tmpIndex in c(func_setID, func_jobID, func_repID)){
								if(func_tmpIndex){
									tmp_returnString <- c(tmp_returnString,thisSplit[1])
									thisSplit <- thisSplit[-1]
								} else {
									tmp_returnString <- c(tmp_returnString,NA)
								}
							}
							return( tmp_returnString )
						}),
						ncol = length(funcBase))
		dimnames(func_tmpReturn) <- list(c("base",
											if(!is.null(func_setID)){if(as.logical(func_setID)){"setID"}else{NULL}}else{NULL},
											if(!is.null(func_jobID)){if(as.logical(func_jobID)){"jobID"}else{NULL}}else{NULL},
											if(!is.null(func_repID)){if(as.logical(func_repID)){"repID"}else{NULL}}else{NULL}),
										funcBase)
		return( func_tmpReturn )
	}
}

# This is a function to calculate the relative fitness for a vector of fitnesses, it can use either an ancestral fitness value
# or it will center the vector of fitness values around 1.  Negative fitness values will be treated as zero.  If the ancestral fitness
# value is zero, then we must assume the fitness values are selection coefficients and we'll center this around 1.
# If the user has passed some value to the weights, we'll calculate a weighted mean fitness.  
# The deault operation, when landscpae is RMF, will be to ignore magnitude of the ancestor and instead we simply calculate the absolute distance from that value
calc_relativeFitness <- function(func_fitVector, func_ancestFit = NULL, func_weights = NULL, func_absDistance = (simModel == "RMF")){
	# If we've passed an ancestral fitness and it is zero ,then we're calculating fitness values as selection coefficients
	# thus we want to use a distance measure for our relative fitness.
	if(!is.null(func_ancestFit)){
		if(func_ancestFit == 0){
			# This is irrespective of the simModel passed as it's based on what the ancestral fitness value was defined to be
			func_absDistance <- TRUE
		}
	}
	# There are some special circumstances to be handled.  Simplest case, is that we've been passed along
	# and ancestral fitness value, that is not zero.  If the ancestral fitness value is zero, we simply ignore it
	# and use our other mehtods of centering values and calculating relative fitness.  Also if we pass an ancestral fitness
	# but it is zero, then we're calculating distributions to represent (s), so we use distance.
	if( (!is.null(func_ancestFit) && func_absDistance) ){
		# We return the distance from the ancestor to a minimum of zero
		return( unlist(lapply(1 + (func_fitVector - func_ancestFit),function(func_thisFitness){ max(0, func_thisFitness) })) )
	} else if(!is.null(func_ancestFit) && func_ancestFit != 0){
		# Simplest case is the ancestral fitness is not zero. This is simple, we just divide the fitness vector by the ancestral fitness value
		# But multiply by it's sign so that the "direction" of fitness increase is preserved in the event of a negative ancestor.  We also consider
		# that negative fitness values
		return( unlist(lapply(sign(func_ancestFit) *func_fitVector/func_ancestFit, function(func_thisFitness){ max(0, func_thisFitness) })) )
		# If the ancestral fitness is zero, then we won't be using it, thus we pass to our other methods. 
	} else {
		# This means we use the fitness vector that was passed will be centered around 1, and negative numbers are treated as zero
					# While this tep would account for the magnitude of values, it's not appropriate in so much as the user ought to be 
					# defining DFE parameters such that the values are reasonably selection coefficients or relatiev fitness values.
					#func_fitVector <- func_fitVector/(max(abs(func_fitVector)))  
		return( unlist(lapply(1 + (func_fitVector - if(is.null(func_weights)){mean(func_fitVector)}else{weighted.mean(func_fitVector,func_weights)}), function(func_thisFitness){ max(0, func_thisFitness) })) )	
	}
}

# This function is used to find which elements of a population matrix are deemed as established, it returns the rows of
# the inMatrix which are established.  If the func_estProp is numeric, then we calculate that a line is established
# if it's popSize > sum(popSize) * func_estProp.  Otherwise we try to evaluate the expression passed.  I have no real
# handlers for the validity of this expression. USER BEWARE!
# A common expression uses the suggestions of Desai 2007 where a lineage establishes when it has 1/s individuals
querryEstablished <- function(func_inMatrix, func_sizeCol = "popSize", func_fitCol = "fitness", func_estProp = 0.01){
	if(is.numeric(func_estProp)){
		# If the value is less than or equal to 1 we take it to mean a proportion of the population
		if(func_estProp <= 1){
			return( func_inMatrix[which(func_inMatrix[, func_sizeCol] >= (func_estProp * sum(func_inMatrix[, func_sizeCol]))),] )
		# Otherwise we take it to mean an exact value, we use the floor call here so that a value of ex: 1.1
		# could be used in order to force return that any lineage which arises is established.
		} else {
			return( func_inMatrix[which(func_inMatrix[, func_sizeCol] >= floor(func_estProp)) ,] )	
		}
	# NOTE: This will only work as long as all fitness values passed are greater than 1.
	} else if(func_estProp == "Desai") {
		return( func_inMatrix[which(apply(func_inMatrix[, c(func_sizeCol, func_fitCol)],MARGIN=1,function(thisRow){
					thisRow[func_sizeCol] >= 1/log(thisRow[func_fitCol])
				})),] )
	} else {
		return( func_inMatrix[eval(parse(text= func_estProp)),] )	
	}
}


# This is a standard critical error response
stopError <- function(func_message){
		stop(func_message)
		traceback()
		Sys.sleep(20)
		q(save="no")
}


# This is a function to trim a string by removing the first and last character, it's used to trim quotation marks
# used in the parameter input
trimQuotes <- function(funcIn){ substr(funcIn,2,nchar(funcIn)-1) }
addQuotes <- function(funcIn){ paste('"',funcIn,'"',sep="") }

# This is a function to control my naming of submission and template scripts
name_subScript <- function(inVar){
	paste("SHAPE_submit_",inVar,".sh",sep="")
}
name_batchSubmit <- function(inVar){
	paste("submit_jobBatch_",inVar,".sh",sep="")
}
name_bodyScript <- function(inVar){
	paste("SHAPE_body_",inVar,".r",sep="")
}
name_parameterScript <- function(inVar){
	paste("SHAPE_parameters_",inVar,".r",sep="")
}

# This is a function for writting out a job's submission script, it allows a submission script to be built for
# either remote or local job submits.
write_subScript <- function(func_subScipt, func_outDir, func_inCombos, func_inParms, func_maxJobs,
							func_appLocation = appLocation, func_commonArgs = commonArgs, func_submitArgs = submitArgs,
							func_remoteLocation = remoteLocation, func_passedArgs = passedArgs, func_sepString = sepString,
							func_externalStopper = external_stopFile){
	func_submitLines <- readLines(func_subScipt,warn=FALSE)
	# If we're doing external selfing we add in some lines to the start and end of the script
	if(as.logical(func_inParms$externalSelfing)){
		# I find the bash line and add a first line after it but before the rest
		func_tmpLine <- which(grepl("#!/bin/bash",func_submitLines,fixed=TRUE))
		func_submitLines <- c(func_submitLines[1:func_tmpLine],
								paste('echo "" > fake_workDir/',func_externalStopper,sep=""),
								func_submitLines[(func_tmpLine+1):length(func_submitLines)],
								paste('if ! [ -e "fake_workDir/',func_externalStopper,'" ]',sep=""),
								'then',
								'	fake_workDir/fakeSubmit',
								'fi')
	}
	# We build the general outline of the server submission script file, but not all elements are included if we're not on the server
	func_submitLines <- gsub("fake_appPath",appLocation,func_submitLines)
	func_submitLines <- gsub("fake_commandArgs",commonArgs,func_submitLines)
	func_submitLines <- gsub("fake_passedArgs",passedArgs,func_submitLines)
	if(as.logical(func_inParms$serverFarm)){
		# We now add in the submission commands to this templated call
		func_tmpLine <- which(grepl("fake_subMit_command",func_submitLines))
		func_submitLines <- c(func_submitLines[1:(func_tmpLine-1)],
								paste("#SBATCH ",func_submitArgs,sep=""),
								func_submitLines[(func_tmpLine+1):length(func_submitLines)])
		func_submitLines <- gsub("fake_serverPath",remoteLocation,func_submitLines)
	}
	
	# now since we're building a submission sccript for each of the job, we loop through them, but tracking how many
	# are placed into any one file so that each submission script has no more than maxJobs worth
	func_tmpCounter <- 0
	func_jobBatch <- 1
	func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
	for(thisCombo in 1:nrow(func_inCombos)){
		for(thisRep in 1:as.numeric(func_inParms$uniqueReplicates)){
			tmp_jobString <- name_batchString(funcBase = trimQuotes(func_inParms$save_batchBase),
												func_setID = thisRep,
												func_jobID = func_inCombos[thisCombo,"save_batchJob"],
												func_sepString = func_sepString)
			# I define this job's working directory
			tmp_jobDir <- paste(func_outDir,
								tmp_jobString,
								"/",sep="")
			tmpSubmit_fileName <- name_subScript(tmp_jobString)
			# I now also build the submission script file for this job, some elements are for the server only
			tmp_submit <- gsub("fake_workDir",tmp_jobDir,func_submitLines)
			tmp_submit <- gsub("fake_tmpScript.r",name_bodyScript(tmp_jobString),tmp_submit)
			
			if(as.logical(func_inParms$serverFarm)){
				tmp_submit <- gsub("fakeOut",paste(tmp_jobString,".o",sep=""),tmp_submit)
				tmp_submit <- gsub("fakeJob",tmp_jobString,tmp_submit)
				tmp_submit <- gsub("fakeDir",tmp_jobString,tmp_submit)
			} else {
				tmp_submit <- gsub("fake_serverPath/fakeDir/",tmp_jobDir,tmp_submit)
			}
			# If we're selfing externally then I update the submission script name
			if(as.logical(func_inParms$externalSelfing)){
				tmp_submit <- gsub("fakeSubmit",tmpSubmit_fileName,tmp_submit)
			}
			
			# I now write out the job's submission script
			writeLines(tmp_submit, con = paste(tmp_jobDir,tmpSubmit_fileName,sep=""))
			# Now I write to the jobBatch submission script
			func_tmpCounter <- func_tmpCounter + 1
			# This updates our counters if need be and differ based on local vs server farming.
			if(as.logical(func_inParms$serverFarm)){
				if(func_tmpCounter > func_maxJobs){
					func_tmpCounter <- func_tmpCounter - (func_maxJobs + 1)
					func_jobBatch <- func_jobBatch + 1
					# We close the old batch script and open a new one
					close(func_masterCon)
					func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
				}
			} else {
				if(func_tmpCounter > ceiling((nrow(func_inCombos) * as.numeric(func_inParms$uniqueReplicates))/func_maxJobs)){
					func_tmpCounter <- 0
					func_jobBatch <- func_jobBatch + 1
					close(func_masterCon)
					func_masterCon <- file(description = paste(func_outDir,name_batchSubmit(func_jobBatch),sep=""), open = "w")
				}
			}
			tmp_masterSubmit <- paste(if(as.logical(func_inParms$serverFarm)){"sbatch "}else{NULL},
										tmp_jobDir,tmpSubmit_fileName,sep="")
			writeLines(tmp_masterSubmit, con = func_masterCon, sep="\n")
		}
	}
	# We close the connection and report completing the tasks.
	close(func_masterCon)
	return( "Done writting out submission scripts" )
}
	
	
# This is a file for updating the post analysis plotting script and creating an updated copy in the experiment's folder
writePlotting <- function(func_infile, func_inParms, func_outFile, func_refSearch = inputReference,
							func_sepString = sepString){
	# Step 1 we read in the lines of the script we're looking to update
	func_tmpLines <- readLines(func_infile, warn=FALSE)
	# I now define the parameters that I want updated in the plotting script
	func_updateParms <- c("serverFarm","baseDir_local","baseDir_server","workDir","save_batchBase","fileName_functions")
	# I update one of the search terms for the baseDir information
	for(thisUpdate in c("serverFarm","baseDir_local","baseDir_server")){
		func_refSearch[[thisUpdate]] <- sub("serverFarm","plot_onServer",func_refSearch[[thisUpdate]])
	}
	func_tmpLines <- updateLines(func_inLines = func_tmpLines, func_searchPattern = c(list("sepString"="sepString <-"),func_refSearch[func_updateParms]),
									func_values = c(list("sepString"=addQuotes(func_sepString)),func_inParms[func_updateParms]))
	# we write out the plotting file
	writeLines(func_tmpLines, con = func_outFile)
	return( "Wrote out plotting file" )
}

# This will write out a parameter file for each job, these will be inherited by the script as it runs
writeParameters <- function(func_infile, func_inParms, func_inCombos, func_outDir, func_bodyScript,
							func_refSearch = inputReference, func_comboRef = comboReference, 
							func_indepRef = notExpanded_reference, func_condRef = conditionReference,
							func_sepString = sepString, func_ExternalStopper = external_stopFile){
	# Step 1 we read in the lines of the script we're looking to update
	func_parmLines <- readLines(func_infile, warn=FALSE)
	func_bodyLines <- readLines(func_bodyScript, warn=FALSE)
	# I now define the parameters that I want updated that are constant throughout jobs
	func_constParms <- c("serverFarm","results_removeSteps","externalSelfing","toggle_forceCompletion",
						"baseDir_local","baseDir_server","save_batchBase","fileName_functions",
						"fileName_preProcessing","maxReplicates","workDir")
	func_parmLines <- updateLines(func_inLines = func_parmLines, 
									func_searchPattern = c(func_refSearch[func_constParms],
															list("externStopper" = "external_stopFile <-")),
									func_values = c(func_inParms[func_constParms],
													list("externStopper"=addQuotes(func_ExternalStopper))))
	# Right, now we need to update and write out a parameter file and body file for each parameter combination.
	for(thisCombo in 1:nrow(func_inCombos)){
		# We now update this parameter file for everything except the unique Set information
		tmp_jobParms <- c("disturbanceType","disturbanceSize","disturbanceSpread",
							"generations_betweenDisturbances","numGenerations","genomeLength","targetNumber",
							"probabilityBirth","growthRate","probabilityDeath","death_byDensity",
							"death_densityFactor","deathRate_sizeCap","scale_births_by_deaths","mutationRate",
							"mutations_only_at_birth","allow_backMutations","growthModel","simModel",
							"const_ancestFitness","constDist","const_distParameters","const_distAsS",
							"const_RMF_initiDistance","const_RMF_theta","const_numInteractions","const_fixedFrame",
							"stochasticBirths","track_only_integer_individuals","size_ofEstablishment","trackingThreshhold")
		tmp_parmLines <- updateLines(func_inLines = func_parmLines, func_searchPattern = func_refSearch[tmp_jobParms],
									func_values = func_inCombos[thisCombo,tmp_jobParms])
		# Now we cycle through each independent statistical replicate and write out files
		for(thisRep in 1:as.numeric(func_inParms$uniqueReplicates)){
			tmp_jobString <- name_batchString(funcBase = trimQuotes(func_inParms$save_batchBase),
												func_setID = thisRep,
												func_jobID = func_inCombos[thisCombo,"save_batchJob"],
												func_sepString = func_sepString)
			# I define this job's working directory
			tmp_jobDir <- paste(func_outDir,
								tmp_jobString,
								"/",sep="")
			tmpParm_fileName <- name_parameterScript(tmp_jobString)
			tmpBody_fileName <- name_bodyScript(tmp_jobString)
			# I now update the body and parm lines with respect to the working directory and jobSet
			tmp_bodyLines <- updateLines(func_inLines =func_bodyLines,
										func_searchPattern = list("source"="sourceParms <-"),
										func_values = list("source"= addQuotes(paste(tmp_jobDir,tmpParm_fileName,sep=""))))
										
			tmp_jobParms <- updateLines(func_inLines =tmp_parmLines,
										func_searchPattern = list("thisSet"="save_batchSet <-",
																	"thisJob"="save_batchJob <-",
																	"selfScript"=c("tmp_selfScript <-","paste(outDir","\"[^\"]*\"")),
										#func_values = list("workDir"= addQuotes(paste(sub(func_outDir,trimQuotes(func_inParms$workDir),tmp_jobDir),sep="")),
										func_values = list("thisSet"= thisRep,
															"thisJob"=func_inCombos[thisCombo,"save_batchJob"],
															"selfScript"=addQuotes(name_subScript(tmp_jobString))))
			# Now I write out both files into their respective job directory
			writeLines(tmp_bodyLines, con = paste(tmp_jobDir,tmpBody_fileName,sep=""))
			writeLines(tmp_jobParms, con = paste(tmp_jobDir,tmpParm_fileName,sep=""))
		}
	}
	return( "Have written out all job's parameters and body scripts" )
}


# This is a function which is used to update lines that are searched and replace in a manner conditional to this script's circumstances
# The input lines can be a vector of any length, and the search patterns can be a list of any length where each list vector is used together.
# The values should be a list of information used as replacement info.
updateLines <- function(func_inLines,func_searchPattern, func_values){
	# There is an asssumption that the searchPattern and values are equally named
	if(!all(is.element(names(func_searchPattern),names(func_values)))){
		stopError("There were not similar searchPattern and values passed to updateLines, please review")
	}
	# For each of the names search patterns we'll try to update the first instance of the line in the code
	func_updateLines <- NULL
	for(thisUpdate in names(func_searchPattern)){
		func_updateLines <- func_inLines
		# The way lines are update depends on the length of the searchPattern information
		func_tmpLine <- which(grepl(func_searchPattern[[thisUpdate]][1],func_updateLines))[1]
		if(length(func_searchPattern[[thisUpdate]]) == 1){
			# This means we simply search for the line and replace it as the pattern and the value pasted
			func_updateLines[func_tmpLine] <- paste(func_searchPattern[[thisUpdate]],func_values[[thisUpdate]][1],sep=" ")
		} else if(length(func_searchPattern[[thisUpdate]]) == 2){
			# If there are two search patterns, the second is used to inform what is replaced in the line
			func_updateLines[func_tmpLine] <- sub(func_searchPattern[[thisUpdate]][2],func_values[[thisUpdate]][1],func_updateLines[func_tmpLine])
		} else if(length(func_searchPattern[[thisUpdate]]) == 3){
			# This means we'll be using the first element to find the line, and the other two to know how to position our change.
			# We use the first and second to identify positions in the string, and then use the region earliest in the string
			# but after the first search term.  I build these as matrices for convenience
			func_tmpPositions <- list(regexpr(func_searchPattern[[thisUpdate]][2],func_updateLines[func_tmpLine],fixed=TRUE),
										gregexpr(func_searchPattern[[thisUpdate]][3],func_updateLines[func_tmpLine]))
			func_tmpPositions <- lapply(func_tmpPositions,function(x){
									matrix(c(unlist(x),attr(if(is.list(x)){x[[1]]}else{x},"match.length")),nrow=2,
											byrow = TRUE,dimnames=list(c("location","length"),NULL))
									})
			# We check that the two positions were found
			if(!all(sapply(func_tmpPositions,function(x){ any(x["location",] > 0) }))){
				stopError(paste("Could not find the proper replacement positions for ",thisUpdate," please review",sep=""))
			}
			# The positions to be used will be the first of the second searches which is larger than the first
			func_usePosition <- which(func_tmpPositions[[2]]["location",] > func_tmpPositions[[1]]["location",1])[1]
			func_tmpPositions <- func_tmpPositions[[2]][,func_usePosition]
			
			# This means the line will become what was previously there, cut around the positions, with the values inserted
			func_updateLines[func_tmpLine] <- paste(substr(func_updateLines[func_tmpLine],1,func_tmpPositions[1]-1),
												func_values[[thisUpdate]],
												substr(func_updateLines[func_tmpLine],sum(func_tmpPositions),nchar(func_updateLines[func_tmpLine])),
												sep="")
		}
		func_inLines <- func_updateLines
	}
	return( func_inLines )
}


# This is a function to take the input parameters and build the parameter combinations
buildCombos <- function(func_inLines, func_comboRef = comboReference, func_indepRef = notExpanded_reference, 
						func_condRef = conditionReference){
	# The func_inLines contain the variety of parameters and the func_refSearch explains how to find and build replacement reference
	
	# Step one, we'll build all possible combinations of parameters that are going to be expanded
	func_parmGrid <- expand.grid(func_inLines[setdiff(names(func_inLines),c(func_comboRef,func_indepRef,unname(unlist(func_condRef))))],
								stringsAsFactors=FALSE, KEEP.OUT.ATTRS = FALSE)
	# Now we build all the combinatorial references but start by ensuring the vectors are of the same length, they should be, 
	# but if the user has not done this we'll recycle the smaller vectors until they are of length for the longest one 
	func_tmpLengths <- sapply(func_inLines[func_comboRef],length)
	if(!all(max(func_tmpLengths) == func_tmpLengths)){
		func_tmpUpdate <- which(func_tmpLengths != max(func_tmpLengths))
		for(thisUpdate in func_tmpUpdate){
			func_inLines[func_comboRef[thisUpdate]] <- list(unname(unlist(rep(func_inLines[func_comboRef[thisUpdate]],max(func_tmpLengths))))[1:max(func_tmpLengths)])
		}
	}
	# We now build on the conditional parameter combinations into a matrix with the conditional combinations
	func_comboAdd <- matrix(unlist(func_inLines[func_comboRef]),nrow=max(func_tmpLengths),dimnames=list(NULL,func_comboRef))
	for(thisCond in names(func_condRef)){
		# We record the columnames prior to updating this is ued in the later steps of this loops
		func_startCols <- colnames(func_comboAdd)
		# We find which rows in the combinations hold the simulation model related to our conditional
		func_tmpLines <- which(grepl(thisCond,func_comboAdd[,"simModel"]))
		# We build the combinations of parameters for this conditional
		func_condAdd <- as.matrix(expand.grid(func_inLines[func_condRef[[thisCond]]],
									stringsAsFactors=FALSE, KEEP.OUT.ATTRS = FALSE))
		# we now add the first condAdd elements to the func_comboAdd matrix
		for(thisAddition in 1:ncol(func_condAdd)){
			# We add the first conditional combination, then any additional required
			func_comboAdd <- cbind(func_comboAdd,func_condAdd[1,thisAddition])
		}
		colnames(func_comboAdd)[(ncol(func_comboAdd)-(ncol(func_condAdd)-1)):ncol(func_comboAdd)] <- colnames(func_condAdd)
		if(length(func_tmpLines) > 0){
			if(nrow(func_condAdd) > 1){
				for(thisRow in 2:nrow(func_condAdd)){
					for(thisUnique in func_tmpLines){
						# We need to add a number of rows equal to the number of func_tmpLines
						func_comboAdd <- rbind(func_comboAdd, c(unlist(func_comboAdd[thisUnique,func_startCols]),unlist(func_condAdd[thisRow,])))
					}
				}
			}
		}
	}
	# Now for each row in func_comboAdd, we need to have nrow(func_parmGrid) replicates of that row
	func_tmpNew1 <- NULL
	for(thisRow in 1:nrow(func_comboAdd)){
		func_tmpNew1 <- rbind(func_tmpNew1,matrix(rep(func_comboAdd[thisRow,],nrow(func_parmGrid)),
													nrow=nrow(func_parmGrid), byrow = TRUE, 
													dimnames=list(NULL,colnames(func_comboAdd))))
	}
	# We now replicate the func_parmGrid a number of times equal to the rows of the other combinations
	func_tmpNew2 <- NULL
	for(thisAddition in 1:nrow(func_comboAdd)){
		func_tmpNew2 <- rbind(func_tmpNew2,func_parmGrid)
	}
	# We now stick the pieces together
	func_parmGrid <- cbind(func_tmpNew2,func_tmpNew1,stringsAsFactors=FALSE)
	# We only keep unique combinations
	func_parmGrid <- unique(func_parmGrid)
	# I now build the starting_parameterCombination vector
	func_parmGrid <- cbind(func_parmGrid,as.numeric(func_inLines$starting_parameterCombination):
											(nrow(func_parmGrid)+(as.numeric(func_inLines$starting_parameterCombination)-1)),
							stringsAsFactors=FALSE)
	colnames(func_parmGrid)[ncol(func_parmGrid)] <- "save_batchJob"
	
	# I now return this
	return( func_parmGrid )
}



###################################################################################################
######################################## END OF FUNCTIONS #########################################
###################################################################################################





###################################################################################################
###################################################################################################
##################################### FUNCTION PRE-COMPILING ######################################
###################################################################################################
###################################################################################################
######### Here we create pre-compiled versions of the most costly functions - Done here so that all constants and pre-definitions are done.
######### This uses the compiler package in R
######### For further icnrease to speed there is the option to use Rcpp, though the data shapes passed through these may make it impractical
comp_fitnessLandscape <- cmpfun(fitnessLandscape)
comp_defineNeighbours <- cmpfun(defineNeighbours)
comp_find_neededNeighbours <- cmpfun(find_neededNeighbours)
comp_createGenotypes <- cmpfun(createGenotypes)
comp_reportPopulations <- cmpfun(reportPopulations)
###################################################################################################
###################################################################################################
###################################################################################################
