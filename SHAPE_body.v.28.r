# This is the SHAPE body file but should not itself be called directly.  Instead use the farming file to actually create and submit jobs for running.
# The program will simulate the population demographics and genotype evolution
# A broad range of parameters are currently supported to allow comparison of dynamics across various models and assumptions

############ DEPENDENCIES ############
# It uses the parameter input file named by the sourceParms object
# Running this script should pass along elements which can be captured by the commandArgs() function 
# which should return at least a 'currReplicate' argument, and possibly an 'outDir'  
# then just the libraries within the source files called below

#### NOTE: This is only established to work for when mutational steps are single mutants, this could be fixed if we address how neighbours are found (see function)

# We tidy up our space and then load libraries
rm(list=ls())

# This is the file path to the parameters that will be used for this run 
sourceParms <- "E:/Research_Data/SimulationRuns/multiMutants/localRun/SHAPEtest/SHAPE_parameters.v.1.r"
if(!file.exists(sourceParms)){
	stop(paste("Could not find the parameters file located at:",sourceParms,sep=" "))
	q(save="no")
} else {
	source(file=sourceParms)
}

# We create a vector that controls which replicates this job should be cycling across,
# if we replicate internally this is between the currReplicate and maxReplicates value
# otherwise it's simply the currReplicate
workingReplicates <- if(externalSelfing){
						currReplicate
					} else {
						seq(currReplicate, maxReplicates,by=1)
					}

# This controls how many, and the number for which, replicates are to be run in this instance.
for(thisRep in workingReplicates){
	####### THESE ARE FILENAME AND OBJECT NAME CONSTANTS USED IN THE SCRIPT #######
	# This defines a naming system that provides unique names for each experiment and parameter combination as controlled by the 
	# batchString and replicate number information.
	### CAUTION: NEVER include digits or underscores (ie: "_" ) in this string literal it may break post analysis.
	save_batchIndex <- sapply(c("Landscape", "Steps", "Parameters", "Neighbourhood"),function(thisIndex){
							paste(save_batchString, thisIndex,if(run_isRecycling[thisIndex]){ recycle_repStart } else { thisRep },sep="_")
						})
	# This is the name of the file output for saving a run
	save_fileName <- paste(paste(save_batchIndex["Parameters"],sep="_"),".RData",sep="")
	# This is a name used for our database files, they'll be unique to this run.
	fileName_dataBase <- list("genotypeSpace"=paste(save_batchIndex["Landscape"],".sqlite",sep=""),
								"timeStep_States"=paste(save_batchIndex["Steps"],".sqlite",sep=""),
								"nearestNeighbours"=paste(save_batchIndex["Neighbourhood"],".sqlite",sep=""))
	
	#################################################################################################
	############################### LANDSCAPE MODEL SPECIFIC INFORMATION ############################
	#################################################################################################
	# If we're recycling our landscape, we must recapture the NK and RMF parameters in order for
	# the runs to be properly interpretable, else we're redefining this.
	if(run_isRecycling["Landscape"] && c(thisRep > recycle_repStart || file.exists(paste(outDir, paste(save_batchString, "Parameters", recycle_repStart,sep="_"),".RData",sep=""))) ){
		# These store information that would be stored in the loaded file but must be maintained as is for this run
		# these objects will be used downstream to restore the values.
		tmp_realReplciate <- thisRep
		tmp_realDist <- track_distSize
		tmp_real_batchStrings <- save_batchIndex
		tmp_realNames <- fileName_dataBase
		
		# We create an environment , into which we can load the existing information so the needed parts can be obtained
		assign("recycleParms", new.env(), pos = ".GlobalEnv")
		if(file.exists(paste(outDir, paste(save_batchString, "Parameters", recycle_repStart,sep="_"),".RData",sep=""))){
			load(paste(outDir, paste(save_batchString, "Parameters", recycle_repStart,sep="_"),".RData",sep=""), envir= recycleParms )
		} else {
			stop("There was a problem trying to load the seed runs' parameters during call to recycle Landscape.  Please review")
		}
		# Now dependeing on wheter or not we're recycling all parameters or just the fitness landscapes
		# which happens to be the DFE, FitnessLandscape, NK_modelElements, RMF_modelElements, sets....
		# we load different sets of parameters
		loadSet <- if(!run_isRecycling["Parameters"]){
						c("FitnessLandscape","DFE","NK_modelElements", "RMF_modelElements")
					} else {
						names(saveParameters)
					}
		
		# we now simply go through the sets to load an place them into the global environment
		for(thisSet in loadSet){
			for(thisNamed in names(recycleParms$runParameters[[thisSet]])){
				assign(thisNamed, recycleParms$runParameters[[thisSet]][[thisNamed]], pos=".GlobalEnv")
			}
		}
		
		# Here is where we restore the necessary information back to its initialised values.
		thisRep <- tmp_realReplciate
		track_distSize  <- tmp_realDist
		save_batchIndex <- tmp_real_batchStrings
		fileName_dataBase <- tmp_realNames
		# We now remove our loaded parameters to save workspace.
		rm(recycleParms)
		
	} else {
		### This section builds fitness landscape model objects for the NK, RMF, Additive and Fixed models ###
		
		# This is the number of interactions a site has in an NK model, this value is not used unless simModel == "NK", but should be 
		# some real positive integer.  The maximum value can be 1 less than the length of the genome.
		const_numInteractions <- min(const_numInteractions,genomeLength - 1)  
		
		# This is the matrix of interactions between sites for the NK model, if we have no interactions among sites, or are not using the NK model
		# we return this as NULL, otherwise we want a matrix with the number of columns equal to the number of interactions
		const_NK_interactionMat <- if(simModel != "NK" || const_numInteractions == 0){
										NULL
									} else if(simModel == "NK" && const_numInteractions == 1){
										as.matrix(sapply(1:genomeLength,function(thisSite){ sample((1:genomeLength)[-thisSite], const_numInteractions,replace=FALSE) }),ncol=1)
									} else {
										t(sapply(1:genomeLength,function(thisSite){ sample((1:genomeLength)[-thisSite], const_numInteractions,replace=FALSE) }))
									}
		# This builds the indepednent effect for each site, based on state, where the wild-type values are zeroes
		const_siteBystate_fitnessMat <- if(is.element(simModel,c("NK","Additive"))){
											sapply(const_siteStates,function(thisState){ 
												# If we're using the NK model then the ancestral state is the ancestral fitness value, for additive
												# this is the selection coefficient.
												return(  if(thisState == "0" && const_relativeFitness){
															 rep(if(simModel == "NK"){
																	const_ancestFitness
																} else if(simModel == "Additive"){
																	0
																}, genomeLength)
															# If we're considering the derived state, then we use the distribution information
														} else {
															tmpReturn <- fitnessDist(genomeLength, tmpDistribution = constDist, tmpParameters = const_distParameters)
															# If the distribution draws relative fitnesses but this needs to be used as selection coefficients...
															if(const_distAsS){
																tmpReturn <- tmpReturn - 1
															}
															tmpReturn
														}
													) })
											} else if(simModel == "Fixed") {
												# This assumes the const_fixedFrame object is a table created by the user and  with the colnames binaryString and fitness
												if(file.exists(paste(workDir,const_fixedFrame,sep=""))){
													const_fixedFrame <- read.csv(paste(workDir,const_fixedFrame,sep=""),header=TRUE,stringsAsFactors = FALSE)
												} else {
													stopError(paste("Could not find the fixed fitness landscape file: ",workDir,const_fixedFrame," please review",sep=""))
												}
												if(all(is.element(colnames(const_fixedFrame),c("binaryString","fitness")))){
													# Ok this means we need to convert what was passed into a named vector.  It is assumed the information passed 
													# for binaryString is a sequence of 0 and 1:  ex:   000,   010,   100,  etc.....
													tmpReturn <- const_fixedFrame[,"fitness"]
													names(tmpReturn) <- unlist(lapply(strsplit(const_fixedFrame[,"binaryString"],""),function(x){
																					if(any(x == "1")){
																						paste(which(x == "1"),sep="",collapse=sepString)
																					} else {
																						""
																					}
																				}))
													tmpReturn
												} else {
													stopError("Fixed fitness landscape information does not contain the two column names binaryString and fitness")
												}
											} else {
												NULL
											}
		# If a siteBystate object was created then we need to give it some column names, also we ensure it is a matrix
		# It may not be a matrix if our genome length = 1.... which is used for trivial purposes but can happen.
		if(is.element(simModel,c("NK","Additive"))){ 
			if(!is.matrix(const_siteBystate_fitnessMat)){
				if(length(const_siteBystate_fitnessMat) != 2){
					print("There was a problem trying to create the <const_siteBystate_fitnessMat> matrix, it's length was ",length(const_siteBystate_fitnessMat),"please review",sep="")
					q(save="no")
				} else {
					const_siteBystate_fitnessMat <- matrix(const_siteBystate_fitnessMat,nrow=1,ncol= length(const_siteStates))
				}
			}
			# Now the const_siteBystate_fitnessMat should be a matrix with a number of  columns equal to the number of states
			colnames(const_siteBystate_fitnessMat) <- const_siteStates 
		}
	
		# Lastly, as a computational convenience for calculating the per genotype fitenss values.  We calculated the per-site, dependent
		# fitness values of the ancestral genotype.  If this way, for each novel genotype we need only update the value of sites which depend
		# upon a mutant sites
		const_DepbySite_ancestFitness <- if(simModel == "NK"){
											sapply(1:nrow(const_siteBystate_fitnessMat),function(thisSite){
													mean(const_siteBystate_fitnessMat[c(thisSite, const_NK_interactionMat[thisSite,]),"0"])
													})
										} else {
											NA
										}				
										
		### These following values are the are only meaningful when the RMF model is implemented
		
		# This value is the weighting of the independent fitness component > 0, but since in the RMF fitness model equation the distance to the 
		# global optima impacts the dynamics so strongly (through c - the independent contribution), and we want to have our dynamics to have some meaningful 
		# relationship, we'll base this value on some idealised theta value (which is a measure of the contribution of independent part and the random component).  
		# From Neidhart 2014 theta is measured as:  theta = c/(sqrt(var(random_component)))
		const_RMF_indWeight <-  if(simModel == "RMF"){
									# We make a large call to our fitness function, find the variance, and then use this to help set our value
									const_RMF_theta * sqrt(var(fitnessDist(5e7, tmpDistribution = constDist, tmpParameters = const_distParameters)))
								} else { 
									NA
								}
		# For RMF some optimal genotype exists and here we'll define what that genotype is.  If the model is not RMF, we return NA
		const_RMF_globalOptima <- sample(1: genomeLength, const_RMF_initiDistance,replace=FALSE)
		const_RMF_globalOptima <- if(simModel == "RMF"){
										paste(const_RMF_globalOptima[order(const_RMF_globalOptima,decreasing=FALSE)],collapse= sepString)
									} else {
										NA
									}
	}
	
	###################################################################################################
	##################################### END OF PREDEFINITIONS #######################################
	###################################################################################################
	
	
	###################################################################################################
	######################################### BEGIN OF BODY ###########################################
	###################################################################################################	
	# We want to report on the run time so initialise that reporter here where the simulations are near to starting.
	startTime <- proc.time()
	
	######## HERE WE ESTABLISH DATABASE CONNECTIONS, THE FITNESS LANDSCAPES, AND REPORTING OBJECTS #######
	# This creates our SQL Lite database which will hold our growing landscape and reporting objects, this is done since our method
	# may (and is intended to) create large objects which are best not held in R directly.  If the connections exist, this resets them.
	connections_dataBase <- sapply(names(fileName_dataBase),function(x){ resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
																						func_type = "connect") })

	##################################################################################################################
	########################## NOW WE ESTABLISH THE INITIAL POPULATION AND OTHER OBJECTS #############################
	##################################################################################################################
	# We initialise the fitness landscape and then the step-wise population tracker database with the initial population.  
	# Here is where if we allow a user to redefine the starting population we'd make changes to SHAPE.
	
	# This is the fitness landscape database
	if(!run_isRecycling["Landscape"]  || c(thisRep == recycle_repStart  && length(dbListTables(connections_dataBase$genotypeSpace)) == 0) ){  
		assign(tmp_genotypeObject,
				create_genotypeFrame(tmp_newID, "", if(const_relativeFitness){ calc_relativeFitness(const_ancestFitness) } else { const_ancestFitness }), 
				pos=".GlobalEnv")
		dbWriteTable(connections_dataBase$genotypeSpace,
					name = if(db_splitTables){nameTable(0,1)}else{nameTable(0)}, 
					value = eval(as.name(tmp_genotypeObject)))
	}
	
	# This is the step-wise population tracking database
	if(!run_isRecycling["Steps"]  || thisRep == recycle_repStart){
		# We name our steps using the step naming function, and setup our population matrix with the population reporting function.
		dbWriteTable(connections_dataBase$timeStep_States,
						name = paste("Step",0,sep="_"), 
						value = comp_reportPopulations(func_numMuts = 0, 
														func_genotypeID = 0, 
														func_popSizes = const_focal_popValue/if(const_growthForm == "logistic"){init_distSteps["factor"]}else{1}, 
														func_fitnesses = if(const_relativeFitness){
																				calc_relativeFitness(const_ancestFitness)
																			} else {
																				const_ancestFitness
																			},
														func_births = const_focal_popValue/if(const_growthForm == "logistic"){init_distSteps["factor"]}else{1}, 
														func_deaths = 0, 
														func_mutants = 0, 
														func_progenitor = "WT"),
						overwrite = TRUE)
	}
	# We initialise the object for tracking disturbance and timing.
	track_distSize <- rbind(c("Step"=0,init_distSteps))
	
	##################################################################################################################
	############################ RESSETTING ANCESTRAL FITNESS VALUE FOR RMF MODEL ####################################
	##################################################################################################################							
	# Now if using the RMF model, we reset the ancestral fitness using the fitness calculation
	# But we have to check if we're recycling our information here, otherwise we'll have loaded this
	if(!run_isRecycling["Landscape"]  || c(thisRep == recycle_repStart  && !file.exists(paste(outDir, save_fileName,sep=""))) ){
		# For the RMF model we'll want the ancestral fitness to be the weighted distance component plus
		# what has been set as the const_ancestFitness value.  We get the distance compoenent only by fudging the distribution
		# as being a null zero distribution of values.  The assumption here is that "const_ancestFitness" was originally set to 
		# reflect the upper end of the distribution passed for future use.
		if(simModel == "RMF"){
			# I assign in here so that the ancestral value is based solely on the distance and not the distribution
			# this makes our calculations of relative fitness, in RMF landscape models, feasible 
			# -- as magnitude was found to at times (when near 0) to be a problem to work with
			# We also suppres the distAsS so that the uniform zero isn't changed to a -1 value.
			assign("const_ancestFitness",comp_fitnessLandscape("", 
																const_ancestFitness, 
																landscapeModel = simModel,
																relativeFitness = FALSE,
																func_distribution = "Uniform",
																func_distParameters = c(0,0),
																func_distAsS = FALSE),
					pos=".GlobalEnv")
		}
	}
	
	# This saves the parameters object file 
	runTime <- 0
	if(!run_isRecycling["Parameters"] || c(thisRep == recycle_repStart && !file.exists(paste(outDir,save_batchIndex["Parameters"],".RData",sep=""))) ){
		# This builds a list of the parameter objects we want saved so the values can easily be recovered later.
		runParameters <- sapply(saveParameters, function(x){ sapply(x,function(y){ eval(as.name(y)) },simplify=FALSE) },simplify="list")
		save(runParameters, runTime, connections_dataBase, fileName_dataBase,  
				file=paste(outDir, save_fileName,sep=""))
	# If we are recycling parameters and we've already saved the repStart parameter file, we'll just load all that information into the global space.
	} else if(run_isRecycling["Parameters"] && file.exists(paste(outDir,save_batchIndex["Parameters"],".RData",sep="")) ){
		# Similar to previous instance of loading, this saves some objects to be restored after the load.
		tmp_realReplciate <- thisRep
		tmp_realDist <- track_distSize
		tmp_real_batchStrings <- save_batchIndex
		tmp_realNames <- fileName_dataBase
		
		# We create an environment, load the values, assign what is needed
		assign("recycleParms", new.env(), pos = ".GlobalEnv")
		if(file.exists(paste(outDir, paste(save_batchString, "Parameters", recycle_repStart,sep="_"),".RData",sep=""))){
			load(paste(outDir, paste(save_batchString, "Parameters", recycle_repStart,sep="_"),".RData",sep=""), envir= recycleParms )
		} else {
			stop("There was a problem trying to load the seed runs' parameters during call to recycle Landscape.  Please review")
		}
		# Even if we're not re-cycling the specific parameters with a logical call, we may need to recover the fitness landscape values.
		loadSet <- if(!run_isRecycling["Parameters"]){
						c("FitnessLandscape","DFE","NK_modelElements", "RMF_modelElements")
					} else {
						names(saveParameters)
					}
		# We now simply go through the sets to load an place them into the global environment
		for(thisSet in loadSet){
			for(thisNamed in names(recycleParms$runParameters[[thisSet]])){
				assign(thisNamed, recycleParms$runParameters[[thisSet]][[thisNamed]], pos=".GlobalEnv")
			}
		}
		# We re-update our real replicate value - this would only matter if run_isRecycling["Parameters"] == FALSE
		thisRep <- tmp_realReplciate
		track_distSize  <- tmp_realDist
		save_batchIndex <- tmp_real_batchStrings
		fileName_dataBase <- tmp_realNames
		# We now remove our loaded parameters to save workspace.
		rm(recycleParms)
	}
	# At this point we check what was the max previously recorded genotypeID and pass this to the run 
	tmp_newID <<- max(unlist(lapply(dbListTables(connections_dataBase$genotypeSpace), function(thisTable){ 
									dbGetQuery(connections_dataBase$genotypeSpace,paste("SELECT MAX(genotypeID) FROM ", thisTable,sep="")) 
							}))) + 1
	
	# We now close the connections to avoid database malformations, recall it was opened to initialise the population.
	for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
	
	# We initialise our disturbance counter and a starting step object that is not currently usefull outside of troubleshooting.
	counter_logSteps <- startingStep <- 1
	
	################################################################################################################################
	################################################### ACTUAL ITERATIVE RUN BODY ##################################################
	################################################################################################################################
	# Now we run through the steps of evolution, this loop is wherein we find the real magic.
	for(thisStep in startingStep:(numGenerations / size_timeStep)){
		# This is simply for reporting reasons, will be found in the stdout
		if(thisStep %% 100 == 0){ 
			print(paste("We're now on thisStep: ",thisStep,sep=""))
			print( proc.time() - startTime )
		}
		# We open the connections, load the population states of the last step, then close the connection.
		connections_dataBase <- sapply(names(fileName_dataBase),function(x){ resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
																						func_type = "connect") })
		workMatrix <- dbReadTable(connections_dataBase$timeStep_States, nameTable_step(thisStep - 1))
		for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
		
		# This loop is where the disturbance tracking occurs, this could be made into a function given it writes to global.
		if(counter_logSteps > track_distSize[nrow(track_distSize),"stepReq"]){
			# We check which populations have at least a half individual otherwise they aren't considered alive.
			currLineages <- as.vector(which(workMatrix[,"popSize"] >= 0.5, useNames=FALSE))
			if(length(currLineages) == 0){
				print(paste("All lineages have died as of timeStep ", thisStep,sep=""))
				break
			}
			# This is an unfortunate requirement to make certain that we've got a workMatrix which is in the proper format
			workMatrix <- comp_reportPopulations(func_numMuts= workMatrix[currLineages,"numMuts"], 
												func_genotypeID= workMatrix[currLineages,"genotypeID"], 
												func_popSizes= workMatrix[currLineages,"popSize"], 
												func_fitnesses= workMatrix[currLineages,"fitness"],
												func_births=workMatrix[currLineages,"births"], 
												func_deaths=workMatrix[currLineages,"deaths"], 
												func_mutants=workMatrix[currLineages,"mutants"], 
												func_progenitor=workMatrix[currLineages,"progenitor"])
			# Now that we've trimmed our workMatrix to the living lineage we can perform lossSampling.
			tmp_lossFactors <- compute_distGrowth(func_distFactor = if(const_growthForm == "logistic"){init_distPars}else{0},
													func_growthType = const_growthForm, 
													func_distType = const_distType, 
													func_growthRate = const_growthRate,
													func_popSize = workMatrix[,"popSize"],
													func_focalSize = const_focal_popValue,
													func_manualGenerations = if(const_growthForm == "constant"){ 
																					numGenerations 
																			} else {
																				const_growthGenerations
																			},
													func_stepDivs = size_timeStep)
			# If the dilution factor is not greater than 1 we won't be performing any dilution.
			if(tmp_lossFactors["factor"] > 1){
				# Now that we've calculated the dilution factor that will be used we call our lossSampling function to see which individuals will remain
				tmp_sampledPop <- lossSampling(func_inPopulation = workMatrix[,"popSize"],
												func_dilutionFactor = 1/if(tmp_lossFactors["factor"] == 0){1}else{tmp_lossFactors["factor"]})
				# We update our loss tracker firstby the amount of loss realised, then by the actual factor of loss this represented
				tmp_lossFactors["popLost"] <- sum(workMatrix[,"popSize"] - tmp_sampledPop)
				tmp_lossFactors["factor"] <- sum(workMatrix[,"popSize"])/sum(tmp_sampledPop)
				workMatrix[,"popSize"] <- tmp_sampledPop
				# We track the population change which occurs as a result of this change
				track_distSize <- rbind(track_distSize, c("Step"= thisStep, tmp_lossFactors))
			}
			# We update our counter_logSteps object back to the initial value.
			counter_logSteps <- 1
		}
		# We check which populations have at least a half individual otherwise they aren't considered alive.
		currLineages <- as.vector(which(workMatrix[,"popSize"] >= 0.5, useNames=FALSE))
		if(length(currLineages) == 0){
			print(paste("All lineages have died as of timeStep ", thisStep,sep=""))
			break
		}
		# This is an unfortunate requirement to make certain that we've got a workMatrix which is in the proper format
		workMatrix <- comp_reportPopulations(func_numMuts= workMatrix[currLineages,"numMuts"], 
											func_genotypeID= workMatrix[currLineages,"genotypeID"], 
											func_popSizes= workMatrix[currLineages,"popSize"], 
											func_fitnesses= workMatrix[currLineages,"fitness"],
											func_births=workMatrix[currLineages,"births"], 
											func_deaths=workMatrix[currLineages,"deaths"], 
											func_mutants=workMatrix[currLineages,"mutants"], 
											func_progenitor=workMatrix[currLineages,"progenitor"])
		# We calculate the number of deaths which occur in our lineage(s)
		num_birthDeaths <- matrix(c(deathFunction(func_inSize = workMatrix[,"popSize"],
													func_inProb = const_deathProb * size_timeStep,
													func_depDensity = death_byDensity,
													func_densityMax = death_densityCap,
													func_densityPower = death_densityCorrelation,
													func_roundValues = track_asWhole)),
									ncol=1)
		# We calculate the number of births for each of our lineage(s)
		num_birthDeaths <- cbind(growthFunction(func_inSize = workMatrix[,"popSize"],
												func_inFitness = workMatrix[,"fitness"],
												func_bProb = const_birthProb,
												func_sizeStep = size_timeStep,
												func_growthForm = const_growthForm,
												func_deaths = num_birthDeaths[,1],
												func_carryingCapacity = const_focal_popValue,
												func_basalRate = const_growthRate,
												func_deathScale = scaleGrowth_byDeaths,
												func_drift = includeDrift,
												func_roundValues = track_asWhole), 
								num_birthDeaths)
		dimnames(num_birthDeaths) <- list(rownames(workMatrix),c("births","deaths"))																																										#print("starting reportPopulations") # This is for reporting on run times but is not needed....
		# Now, because the number of births may be negative (for reasons of constant growth or growth function with carrying capacities adjusting sizes - its rare...)
		# we adjust the num_birthDeaths so that negative births are shuffled to deaths.
		tmpUpdates <- which(num_birthDeaths[,"births"] < 0)
		if(length(tmpUpdates) > 0){
			# We subtract the negative births from (basically add to...) the amount of deaths
			num_birthDeaths[tmpUpdates,"deaths"] <- num_birthDeaths[tmpUpdates,"deaths"] - num_birthDeaths[tmpUpdates,"births"]
			num_birthDeaths[tmpUpdates,"births"] <- 0
		}
		
		# Now is the section where we verify the number of mutants which may arise in each lineage
		# We check if there are any lineages which could be selected
		tmp_mutableLineages <- if(!allow_backMutations && is.element(genomeLength,workMatrix[,"numMuts"])){ 
									# If we don't allow back mutations and there is a mutant which is fully derived we don't consider it
									rownames(num_birthDeaths)[which(workMatrix[rownames(num_birthDeaths),"numMuts"] != genomeLength)]
								} else {
									# Otherwise all lineages can be considered
									rownames(num_birthDeaths)
								}
		# A lineage cannot be mutable if, after births and deaths, the population size is not at least 1
		tmp_mutableLineages <- tmp_mutableLineages[which(workMatrix[tmp_mutableLineages,"popSize"] + 
														num_birthDeaths[tmp_mutableLineages,"births"] - 
														num_birthDeaths[tmp_mutableLineages,"deaths"] >= 1)]
		
		# So for all lineages which could generate a mutant, we call the mutation process to define the number which arise
		# Since a birth events include not just the offpsring but also the parental type, we multiply the number which
		# may have gained a mutation during replication by the basic growth rate value.  If mutations occur not only at birth
		# then we adjust the pop-Size by the (births * (basic growth rate - 1))
		proposedMutants  <- mutationFunction(func_inSize = unlist(lapply(num_birthDeaths[tmp_mutableLineages,"births"] * const_growthRate/(const_growthRate - 1),function(x){ max(0,x) })),
												func_inProb = const_mutProb) + 
								if(muts_onlyBirths){
									0
								} else {
									mutationFunction(func_inSize = unlist(lapply(workMatrix[tmp_mutableLineages,"popSize"] - 
																				num_birthDeaths[tmp_mutableLineages,"deaths"] -
																				num_birthDeaths[tmp_mutableLineages,"births"] * 1/(const_growthRate - 1), function(x){ max(0,x) })),
													func_inProb = const_mutProb * size_timeStep)
								}						
		names(proposedMutants) <- tmp_mutableLineages
		# We now udpate our num_birthDeaths object with the mutants generated but noting that a lineage 
		# cannot generate more mutants than it has living memebers.
		num_birthDeaths <- cbind(num_birthDeaths,"mutants"= as.vector(sapply(rownames(workMatrix),function(x){ 
																	if(is.element(x,names(proposedMutants))){
																		return( min(proposedMutants[x],
																					workMatrix[x,"popSize"] + 
																					num_birthDeaths[x,"births"] - 
																					num_birthDeaths[x,"deaths"]) )
																	} else {
																		0
																	} })))
		# We initialise the object for reporting on which mutants are created.
		all_newMutants <- NULL
		if(sum(num_birthDeaths[,"mutants"]) > 0){
			for(thisLineage in rownames(num_birthDeaths)[which(num_birthDeaths[,"mutants"] > 0)] ){
				# This re-opens the connections so we can access the fitness landscape object.
				connections_dataBase <- sapply(names(fileName_dataBase),function(x){ resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
																						func_type = "connect") })
				# We define the focal genotype we're working with
				this_focalGenotype <- retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
															func_numMuts = workMatrix[thisLineage,"numMuts"])
				# This is an unfortunate work-around for an instance where I've found that a lineage had been miss-classified by numMuts
				# which caused there to be no value returned, I similarly check for multiple instances....
				# I believe to have fixed all situations that would have led to this
				if(length(this_focalGenotype) == 0 || nrow(this_focalGenotype) != 1){
					# I'm interested in how often this sanity check is used...
					print(paste("Safety net for mis-classified genotypeID has been used for ", this_focalGenotype," on step ",thisStep,sep=""))
					# We look through all of the tables and return the one who has a search result for our genotypeID
					tmp_mutSearch <- sapply(nameTable(dbListTables(connections_dataBase$genotypeSpace), func_splitName = TRUE), function(this_numMuts){
													return( nrow(retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
																								func_numMuts = as.numeric(this_numMuts))) )
										})
					# Now if we've found only one table which held our value, we update the workMatrix's numMuts and re-run our call.
					if(length(which(tmp_mutSearch == 1)) == 1){
						workMatrix[thisLineage,"numMuts"] <- as.numeric(names(tmp_mutSearch)[which(tmp_mutSearch == 1)])
						# We now update the value for this_focalGenotype, as we should have been able to define a proper table.
						this_focalGenotype <- retrieve_binaryString(func_genotypeID = as.numeric(thisLineage),
																		func_numMuts = workMatrix[thisLineage,"numMuts"])
					} else { 
						stop(paste("There was a problem while trying to define a unique table containing ",thisLineage," during ",thisStep," please review",sep=""))
					}
				}
				# We check to see if thisLineage is already within the neighbourhood reference database, this is a time saving database.
				tmp_allNeighbours <- NULL
				if(is.element(nameTable_neighbourhood(thisLineage),dbListTables(connections_dataBase$nearestNeighbours))){
					tmp_allNeighbours <- dbGetQuery(connections_dataBase$nearestNeighbours,paste("SELECT * FROM ",nameTable_neighbourhood(thisLineage),sep=""))$neighbours
				} else {
					# We define all the possible nearestNeighbours for this lineage since they've not been stored.
					tmp_allNeighbours <- comp_defineNeighbours(func_tmpGenotype = this_focalGenotype[1,"binaryString"], func_tmpDirection = allow_backMutations)
					# Now if the size of thisLineage is large enough we store this in the Neighbourhood shortcut database.
					if(workMatrix[thisLineage,"popSize"] >= const_hoodThresh){
						dbWriteTable(connections_dataBase$nearestNeighbours,
									nameTable_neighbourhood(thisLineage),
									data.frame("neighbours"=tmp_allNeighbours))
					}
				}
				
				# We see if this mutant's neighbourhood has been explored, if not we need to add it to the fitness landscape
				if(!as.logical(this_focalGenotype[1,"isExplored"])){
					# This uses the fitness landscape models to calculate the fitness value for all mutants in this neighbourhood.
					comp_createGenotypes(tmp_focalGenotype = this_focalGenotype[1,"binaryString"], 
										tmp_focalFitness = if(is.element(simModel,c("NK","Additive"))){ const_siteBystate_fitnessMat } else { workMatrix[thisLineage,"fitness"] }, 
										maxHamming = max_numMutations, 
										tmp_landModel = simModel, 
										tmp_relativeFitness = const_relativeFitness,
										tmpDistribution = constDist, 
										tmpParameters = const_distParameters,
										tmp_currNeighbours = tmp_allNeighbours,
										tmp_genCon = connections_dataBase$genotypeSpace,
										tmp_tableSplit = db_splitTables,
										tmp_distAsS = const_distAsS)
				}
				# Now we generate our mutants by sampling all possible neighbours with replacement
				newMutants <- table(sample(tmp_allNeighbours,
											num_birthDeaths[thisLineage,"mutants"], replace = TRUE))
				# Now for the mutants that were drawn from the neighbourhood, we go and retrieve the information from the fitness landscape database.
				# To know where to find our mutant information we find the number of mutations for the mutant, then go search those for the genotype.
				tmp_numMuts <- unlist(lapply(strsplit(names(newMutants),sepString),length))
				tmp_dbTables <- dbListTables(connections_dataBase$genotypeSpace)
				tmp_dbTables <- tmp_dbTables[unique(unlist(lapply(nameTable(unique(tmp_numMuts)),function(thisString){ which(grepl(thisString,tmp_dbTables)) })))]
				# Here we get the information for the mutant by querrying our SQL database
				# This is the actual string for the call to information
				tmp_newStrings <- gsub("[[:space:]]","",paste("\'", names(newMutants),"\'",collapse=','))
				# This is the call to get our information about mutants.  It's followed by a check that all the relevant information was found
				# If it wasn't we try to create the missing pieces.
				tmp_infoAdd <- as.matrix(dbGetQuery(connections_dataBase$genotypeSpace, paste("SELECT genotypeID,fitness,binaryString FROM ", 
																								tmp_dbTables,
																								' WHERE binaryString IN (',
																								tmp_newStrings,
																								')',collapse=" UNION ")))
				# If we didn't find some of the mutant information we try again after trying to re-fill the fitness landscape
				# This ought not to be an issue or be required.
				if(!all(is.element(names(newMutants),tmp_infoAdd[,"binaryString"]))){
					# We explore the neighbouring space of this focal lineage
					comp_createGenotypes(tmp_focalGenotype = this_focalGenotype[1,"binaryString"], 
									tmp_focalFitness = if(is.element(simModel,c("NK","Additive"))){ const_siteBystate_fitnessMat } else { workMatrix[thisLineage,"fitness"] }, 
									maxHamming = max_numMutations, 
									tmp_landModel = simModel, 
									tmp_relativeFitness = const_relativeFitness,
									tmpDistribution = constDist, 
									tmpParameters = const_distParameters,
									tmp_currNeighbours = tmp_allNeighbours,
									tmp_genCon = connections_dataBase$genotypeSpace,
									tmp_tableSplit = db_splitTables,
									tmp_distAsS = const_distAsS)
					# We go to collect the information again... it should work now.
					tmp_infoAdd <- as.matrix(dbGetQuery(connections_dataBase$genotypeSpace, paste("SELECT genotypeID,fitness,binaryString FROM ", 
																								tmp_dbTables,
																								' WHERE binaryString IN (',
																								tmp_newStrings,
																								')',collapse=" UNION ")))
				}
				for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
				# Now we should, no matter what, have information for our mutant and add it to the mutant tracking matrix
				tmp_newOrder <- unlist(lapply(names(newMutants),function(x){ which(tmp_infoAdd[,"binaryString"] == x) }))					
				newMutants <- cbind(tmp_numMuts, tmp_infoAdd[tmp_newOrder,"genotypeID"], newMutants, tmp_infoAdd[tmp_newOrder,"fitness"])
				colnames(newMutants) <- popMat_colnames
				rownames(newMutants) <- as.character(newMutants[,"genotypeID"])
				# This updates the mutant tracking matrix 
				all_newMutants <- rbind(all_newMutants, 
										cbind(newMutants,
												"progenitor"= unlist(lapply(newMutants[,"popSize"],function(thisSize){
																		paste(thisLineage,thisSize,sep=sepString)
																})) 
												) 
										)
			} # This closes out the thisLineage split for each row of num_birthDeaths which has a mutant
		} # This closes out the conditional that there are mutants to have entered this loop.
		
		# We update our tmp_stepChanges with the number of births, deaths and mutations which occured for our existing (workMatrix) lineages
		# We also report on new mutants and which lineages birthed them new mutants, this is achieved through our progenitorReport object update in the newMutants section
		# We use the progenitorReport object to easily reference which lineages 
		
		# First we report on all the existing lineages, this is the easy part. 												
		tmp_stepChanges <- data.frame(cbind(workMatrix[rownames(num_birthDeaths), 
											popMat_colnames], num_birthDeaths, 
											"progenitor"="",stringsAsFactors=FALSE)
										,stringsAsFactors=FALSE)
		# If there are no newMutants we can ignore all of this since the changes are just birth based...
		if(!is.null(all_newMutants)){
			# For reasons I can't track, if two separate progenitor lineages produce the same mutant genotypeID
			# I've found that the ID's are not always tracked properly.  This fixes that.
			all_newMutants[,"genotypeID"] <- gsub("[[:space:]]","",all_newMutants[,"genotypeID"])
			# Then we check if any of these newMutants come from different progenitors which would mean the same ID appears twice.
			tmp_repIDs <- table(all_newMutants[,"genotypeID"])
			tmp_repIDs <- names(which(tmp_repIDs > 1))
			# If there are any replicated genotypeID's then we'll want to merge those rows
			if(length(tmp_repIDs) > 0){
				tmp_removeRows <- NULL
				for(thisID in tmp_repIDs){
					tmpRows <- as.vector(which(all_newMutants[,"genotypeID"] == thisID))
					# We update the first row to be the combination of all rows
					all_newMutants[tmpRows[1],"popSize"] <- as.character(sum(as.numeric(all_newMutants[tmpRows,"popSize"])))
					all_newMutants[tmpRows[1],"progenitor"] <- paste(all_newMutants[tmpRows,"progenitor"],collapse= collapseString)
					# then flag the other rows to be removed
					tmp_removeRows <- c(tmp_removeRows, tmpRows[-1])
				}
				# Checking there are rows to remove is a sanity check but should always be true
				if(!is.null(tmp_removeRows)){
					all_newMutants <- matrix(all_newMutants[-tmp_removeRows,],ncol=ncol(all_newMutants),
												dimnames=list(rownames(all_newMutants[-tmp_removeRows]),colnames(all_newMutants)))
				}
			}
			# Now since we're about to add our all_newMutants information to the workMatrix and tmp_stepChanges data.frames, we put it in the same format
			all_newMutants <- data.frame(all_newMutants,stringsAsFactors=FALSE)
			for(thisCol in popMat_colnames){
				all_newMutants[, thisCol] <- as.numeric(all_newMutants[, thisCol])
			}													
															
			# Now we update the stepChanges information by looking at which mutants are new to the population or existing
			# For those which exist we simply update the popSize, for new ones we add rows
			tmpTypes <- is.element(all_newMutants[,"genotypeID"], tmp_stepChanges[,"genotypeID"])
			# For the TRUE responses we're updating our existing matrix so we find the similarly ordered rows
			tmpExisting <- if(any(tmpTypes)){ 
								unlist(lapply(all_newMutants[which(tmpTypes),"genotypeID"], function(x){ which(tmp_stepChanges[,"genotypeID"] == x) }))
							} else { 
								NULL 
							}
			# If anything already exists then we go through each and simply update the tmp_stepChanges row(s)
			if(!is.null(tmpExisting)){
				# Since the tmpExisting is based on the order of the which(tmpTypes), we can call the additions, or combinations, 
				# as vectors of the matrices
				
				# We update the births to be the sum of the existing births and the newly arived mutants
				tmp_stepChanges[tmpExisting,"births"] <- apply(cbind(tmp_stepChanges[tmpExisting,"births"],
																	all_newMutants[which(tmpTypes),"popSize"]),MARGIN=1,sum)
				# We update the progenitor information to be the existing information plus the newly arrived mutants.
				tmp_stepChanges[tmpExisting,"progenitor"] <- apply(cbind(tmp_stepChanges[tmpExisting,"progenitor"],
																			all_newMutants[which(tmpTypes),"progenitor"]), MARGIN = 1, function(x){ 
																				paste(x,collapse=collapseString)
																			})
				# We now remove the which(tmpTypes) and ensure it is still a data.frame of proper typed data.
				all_newMutants <- all_newMutants[-which(tmpTypes),]	
			}
			
			# For mutant genotypeID's that aren't in the population we are adding information to the population matrix.
			tmp_stepChanges <- rbind(tmp_stepChanges, 
									comp_reportPopulations(func_numMuts= all_newMutants[,"numMuts"], 
															func_genotypeID= all_newMutants[,"genotypeID"], 
															func_popSizes= all_newMutants[,"popSize"], 
															func_fitnesses= all_newMutants[,"fitness"],
															func_births= all_newMutants[,"popSize"], 
															func_deaths= rep(0,nrow(all_newMutants)), 
															func_mutants= rep(0,nrow(all_newMutants)), 
															func_progenitor= all_newMutants[,"progenitor"]))
		}
		# Now we update the popSizes in thisStep where the popSize for a lineage changes by births - deaths, mutants are already accounted for above.
		# We store our changes in a reporting matrix so that we adjust our large report_stepStates only once.
		# NOTE: We only need to update the births and deaths for those lineage which existed before mutations occured, as new mutants have popSize given by their generation
		tmp_stepChanges[rownames(num_birthDeaths),"popSize"] <- tmp_stepChanges[rownames(num_birthDeaths),"popSize"] + 
																tmp_stepChanges[rownames(num_birthDeaths),"births"] - 
																tmp_stepChanges[rownames(num_birthDeaths),"deaths"] - 
																tmp_stepChanges[rownames(num_birthDeaths),"mutants"]
		
		connections_dataBase <- sapply(names(fileName_dataBase),function(x){ 
											resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
															func_type = "connect") 
										})
		# We write out the state of the population at the end of this step.
		dbWriteTable(connections_dataBase$timeStep_States,
						name = nameTable_step(thisStep),
						value = comp_reportPopulations(func_numMuts= tmp_stepChanges[,"numMuts"], 
														func_genotypeID= tmp_stepChanges[,"genotypeID"], 
														func_popSizes= tmp_stepChanges[,"popSize"], 
														func_fitnesses= tmp_stepChanges[,"fitness"],
														func_births= tmp_stepChanges[,"births"], 
														func_deaths= tmp_stepChanges[,"deaths"], 
														func_mutants= tmp_stepChanges[,"mutants"], 
														func_progenitor= tmp_stepChanges[,"progenitor"]))
		# We close the connection.
		for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
		# We now advance our step counter to track the disturbance calls.
		counter_logSteps <- counter_logSteps + 1
	} 
	
	###################################################################################################
	########################################## END OF BODY ############################################
	###################################################################################################
	
	
	
	###################################################################################################
	################################### BEGIN OF SUMMARY TOOLS ########################################
	###################################################################################################
	
	# Track the run time value so it can be reported.
	runTime <- proc.time() - startTime
	
	# This saves the parameters, again, only difference here is that the runTime is now completed.
	if(!run_isRecycling["Parameters"] || thisRep == recycle_repStart){
		# This is the list that holds the parameter information
		runParameters <- sapply(saveParameters, function(x){ sapply(x,function(y){ eval(as.name(y)) },simplify=FALSE) },simplify="list")
		save(runParameters, runTime, connections_dataBase, fileName_dataBase,  
				file=paste(outDir, save_fileName,sep=""))
	}
	
	connections_dataBase <- sapply(names(fileName_dataBase),function(x){ 
									resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
													func_type = "connect") 
								})
	# Now for all unexplored mutational spaces around existing genotypes, we define that neighbourhood.
	tmpTables <- dbListTables(connections_dataBase$timeStep_States)
	all_lastGenotypes <- dbReadTable(connections_dataBase$timeStep_States, 
									dbListTables(connections_dataBase$timeStep_States)[which.max(unlist(lapply(strsplit(tmpTables,"_"),function(x){ as.numeric(x[2]) })))])
	# We check if there are any established lineages that need to have their local neighbourhood defined.
	establishedGenotypes <- querryEstablished(func_inMatrix= all_lastGenotypes, func_estProp = const_estProp)
	if(nrow(establishedGenotypes) > 0){
		# If a mutant is not found within the genotypeID's of the database of stepChanges, then we need to create is local nieghbourhood.
		for(thisGenotype in 1:nrow(establishedGenotypes)){
			tmp_genotypeInfo <- retrieve_binaryString(func_genotypeID = establishedGenotypes[thisGenotype,"genotypeID"],
															func_numMuts = establishedGenotypes[thisGenotype,"numMuts"])
			# We ask if this genotype has been explored
			if(!as.logical(tmp_genotypeInfo[1,"isExplored"])){
				comp_createGenotypes(tmp_focalGenotype = tmp_genotypeInfo[1,"binaryString"],
									tmp_focalFitness= if(is.element(simModel,c("NK","Additive"))){ const_siteBystate_fitnessMat } else { establishedGenotypes[thisGenotype,"fitness"] },
									maxHamming = max_numMutations, 
									tmp_landModel = simModel, 
									tmp_relativeFitness = const_relativeFitness,
									tmpDistribution = constDist, 
									tmpParameters = const_distParameters,
									tmp_genCon = connections_dataBase$genotypeSpace,
									tmp_tableSplit = db_splitTables,
									tmp_distAsS = const_distAsS)
			}
		}
	}
	for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
	
	# This is a reporting call, it will be in the stdout.
	if(toggle_forceCompletion  && thisStep != (numGenerations / size_timeStep)){ 
		print("There was a problem and the main body run did not complete, please review")
		q(save="no")
	} else if(thisStep == (numGenerations / size_timeStep)) {
		print(" Simulations completed without crashing ") 
	} else {
		print(" Simulation code reached end of loop ") 
	}
	
	######### Now the results of the run are processed. ##############
	connections_dataBase <- sapply(names(fileName_dataBase),function(x){ resetDatabase(paste(outDir, fileName_dataBase[[x]],sep=""),
																						func_type = "connect") })
	# This should be the script with the pre-processing functions.
	source(fileName_preProcessing)
	# Now we run the processing function noting that the <processedData_fileName> and <processedObjects> 
	# are actually loaded from the source file's parameters
	tmpReport <- runProcessing(func_stepsCon = connections_dataBase[["timeStep_States"]], 
								func_landscapeCon = connections_dataBase[["genotypeSpace"]],
								func_hoodCon = connections_dataBase[["nearestNeighbours"]],
								func_estProp = const_estProp,  
								func_saveFile = processedData_fileName, 
								func_processObjects = processedObjects,
								func_hoodPriority = const_hoodDepth)
	print( tmpReport )
	# Now that we've pre-processed we delete the Steps file if asked to do so, but only if we've been given a flag that the processing was completed
	if(results_removeSteps && grepl("Processing completed", tmpReport) && file.exists(processedData_fileName)){
		file.remove(paste(outDir, fileName_dataBase[["timeStep_States"]],sep=""))
	}
	for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
	# Now if we're not recycling our neighbourhood database we'll be removing it from existence
	if(!run_isRecycling["Neighbourhood"]){ file.remove(paste(outDir, fileName_dataBase[["nearestNeighbours"]],sep="")) }
	
	###################################################################################################
	#################################### END OF SUMMARY TOOLS #########################################
	###################################################################################################
	
	
	
	###################################################################################################
	################################### BEGIN OF SELFING TOOLS ########################################
	###################################################################################################
	# This is where we decide whether or not the script is replicating a new version of itself.
	# But ignore this if we're running it internally controlled by the logical toggle "externalSelfing"
	if(thisRep < maxReplicates && externalSelfing){	
		# Then we read in the selfing template script, but we check that it exsits.
		if(!file.exists(tmp_selfScript)){
			stopError("There was a problem while trying to self, as the template did not exist.  Please review")
		} else {
			# Now it depends on if we're selfing on the remote CAC server or doing a local job
			# NOTE: It's assumed that the tmp_selfScript is keyed for the proper type of job
			# This is the name for the job we want to create
			tmp_nextRep <- thisRep + 1
			tmp_jobName <- name_batchString(funcBase = save_batchString,
											func_repID = tmp_nextRep,
											func_sepString=sepString)
			# I simply update the selfing scripts
			tmpScript <- readLines(tmp_selfScript)
			# We update the currReplicate information and write out the file
			tmpScript <- updateLines(func_inLines = tmpScript,
									func_searchPattern = list("rep"=c("currReplicate=","currReplicate=([[:digit:]])+")), 
									func_values = list("rep"=paste("currReplicate=",tmp_nextRep,sep="")))
			# We write out the lines and make the script executable
			writeLines(tmpScript, con=tmp_selfScript)
			Sys.chmod(tmp_selfScript,mode="0777")
		}
		# To continue the run we also remove the stop trigger
		if (file.exists(paste(finalDir,external_stopFile,sep=""))){
			file.remove(paste(finalDir,external_stopFile,sep=""))
		}
	}
	# If we're selfing internally, but using a remote server to farm, then I'll want - every so often - to transfer completed work back "home"
	if (serverFarm && !externalSelfing) {
		# Every twenty replicates we'll make a call to transfer out files back to the original directory of the run
		if(thisRep %% 20 == 0){
			system(paste("cd ",outDir," ; cp *.RData *.sqlite *.Rout *.o ",finalDir,sep=""))
		}
	}
	
	###################################################################################################
	#################################### END OF SELFING TOOLS #########################################
	###################################################################################################
	
} # This closes the for loop of potential internal replicating as define by the vector workingReplicates 
	# and the logical toggle externalSelfing
	
####################################
######## THAT'S ALL FOLKS ##########
########       CHEERS     ##########
####################################
q(save="no")


