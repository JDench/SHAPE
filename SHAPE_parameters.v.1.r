# This is the actual file that can be directly changed by users who want to directly write new parameters for a run.
# This is the file that will contain the parameters used in common between SHAPE scripts - EXCEPTION filepathing can be independently set for the plotting
# Alternatively, and recommended, is that the SHAPE_farmerParms.v.#.r file be used to set parameters in ranges and then 
# The SHAPE_farmer.v.#.r be run to build all scripts and files for a run.
### NOTE: There are parameters that can be set in this file that cannot be set elsewhere 
###			BUT that's because it's not expect anyone should ever need to unless they're re-coding SHAPE.	
################################################################################################
####################################### PARAMETER VALUES #######################################
################################################################################################


################################################################################
############################ DATABASE AND RECORDING ############################
################################################################################
# Some of these parameters can be changed, others are required for downstream post analysis, but there is no obvious reason why a user would want/need to make changes here.
# This is the string used to separate items
sepString <- "_"
# This is a separator used to collapse our storage of progenitor and number of new mutants (births) values
collapseString <- "__:__"
# This is the maximum number of rows we want to have any one table in our database
maxRows <- 2.5e7
# This logical toggle is if we want to sub-index tables - ie limit their number of rows with maxRows
# Setting this to FALSE will mean maxRows has no real meaning, 
# NOTE:  From experience not splitting rows becomes a problem when individual tables are HUGE
#		But I've not tested issues that may arise from my UNION calls when there are too many tables....
#		Consider a genotypeDB where one numMuts has > 1e9 entries, this would mean 1000 tables when using
#		maxRows at 1e6 so set max rows quite high.  I've found that 2.5e7 is acceptably large:
# 		It's a compromise between R's "comfort" with large vectors and having only 40 tables per 1e9 genotypeIDs 
db_splitTables <- TRUE
# This is a constant naming string used for tables
string_tableNames <- "numMutations"
# This is the genotypeID initial value
tmp_newID <- 0
# this is the colnames we want placed in our population at timeStep matrices
popMat_colnames <- c("numMuts","genotypeID","popSize","fitness")
reportMat_colnames <- c(popMat_colnames,"births","deaths","mutants","progenitor")
# This is the name of temporary genotype space tables we may create or collected locally:
tmp_genotypeObject <- "tmp_genotypeSpace"

# This is the vector of constants and parameters that we want to save for future reference (see file saving in tools)
# We include some 
saveParameters <- list("Population"= c("const_focal_popValue", 
										"genomeLength", 
										"const_mutProb",
										"muts_onlyBirths",
										"numGenerations",
										"size_timeStep",
										"const_deathProb",
										"death_byDensity",
										"death_densityCorrelation",
										"death_densityCap",
										"const_ancestFitness",
										"const_estProp"),
						"Growth_Disturbance"= c("const_distType",
												"const_growthGenerations",
												"init_distPars",
												"track_distSize",
												"const_growthForm", 
												"const_growthRate",
												"scaleGrowth_byDeaths",
												"includeDrift"),
						"FitnessLandscape"= c("simModel",
												"max_numMutations",
												"const_relativeFitness",
												"allow_backMutations",
												"const_distAsS"),
						"NK_modelElements" = c("const_numInteractions",
												"const_NK_interactionMat",
												"const_siteBystate_fitnessMat"),
						"RMF_modelElements" = c("const_RMF_theta",
												"const_RMF_indWeight",
												"const_RMF_initiDistance",
												"const_RMF_globalOptima"),
						"Fixed_Landscape"=c("const_fixedFrame"),
						"DFE"= c("constDist",
								"const_distParameters"),
						"DataManagement" = c("sepString",
											"collapseString",
											"string_tableNames",
											"db_splitTables",
											"maxRows",
											"thisRep",
											"maxReplicates",
											"save_batchBase",
											"save_batchSet",
											"save_batchJob",
											"save_batchString",
											"save_batchIndex",
											"fileName_dataBase",
											"fileName_functions",
											"fileName_preProcessing"))
	
	
###################################################################################################
################################# BEGIN OF CONSTANT DEFINITIONS ###################################
###################################################################################################
# This section includes those objects/variables that the user should predefine before running the script
serverFarm <- FALSE
# This is a secondary logical toggle as to whether or not this script will delete it's Steps file after pre-processing
results_removeSteps <- TRUE
# This is a logical toggle as to whether the script replicates by external job creation or internally with for loops. 
externalSelfing <- FALSE
external_stopFile <- "someNamed.file"


####### FILE AND SYSTEM LOCATIONS #######
baseDir <- if(serverFarm){ "/global/home/yourName/"	} else if(!serverFarm) { "E:/YourPath/" }
workDir <- paste(baseDir,"A_Folder/Some_subDir/",sep="")

# This loads our functions which we're storing in a separate file.
fileName_functions <- paste(workDir,"sourceFunction_SHAPE_file.r",sep="")
if(!file.exists(fileName_functions)){
	stop(paste("Could not find the source file at location:",fileName_functions,sep=" "))
	q(save="no")
} else {
	source(file=fileName_functions)
}
# This is a string that can be used to uniquely identify the batch of jobs, something must be defined for funcBase, 
# nothing is defined for the rep as those are within a job's directory.  The jobID and setID values can be left as NULL
# and if so they will be later interpreted as such to indicated nothing was used in naming the file
#### NOTE: DO NOT INCLUDE THE sepString STRING PATTERN ANYWHERE IN THE FOLLOWING THREE OBJECTS!
save_batchBase <- "yourJob"
save_batchSet <- 1
save_batchJob <- 1
save_batchString <- name_batchString(funcBase = save_batchBase,
									func_setID = save_batchSet,
									func_jobID = save_batchJob,
									func_sepString = sepString)
									
outDir <- paste(workDir, save_batchString,"/",sep="")

# Now, if we're doing processing on a remote server with SLURM submissions, then we may have been 
# passed an outDir argument that is meant for the compute node location, in which case we'll need to save
# a new object for the final repository
finalDir <- outDir 
# This is the definition of a fileName which should offer a template from which to issue new runs
# For remote server this should be a .sh file whereas for local runs it should be the actual .r script.
tmp_selfScript <-  paste(outDir,"selfingTest_internal.r",sep="")

	

# Now inherit arguments which have been supplied at the command line, this is before the general pre-definitions thus it requires that I have hashed out their call
# in the main body of this script.
### NOTE: the command line must be called as R CMD BATCH '--args <arg1> <arg2> <arg3> ...' <filepath>.r where each each argument
# must not have space within, as they are space separated.
#### FOR THIS SCRIPT: I assume a currReplicate argument is being passed.
outsideArgs <- commandArgs(trailingOnly = TRUE)
# This checks if any arguments were passed, if not we warn that script defaults are used,
if(length(outsideArgs)==0){
    replicate(5,print("No arguments supplied, please review"))
} else {
	# If arguments are passed, then they are set as per the exact argument passed along through batch call
    for(i in 1:length(outsideArgs)){
    	# This evaluates the arguments being passed, which means as long as they were assignments the values should be written as objects.
         eval(parse(text= outsideArgs[[i]]))
    }
}

print(outDir)

# We query if the outDir has a trailing "/" or "\\" value
if(substr(outDir, nchar(outDir),nchar(outDir)) != "/"){
	outDir <- paste(outDir,"/",sep="")	
}

if(!dir.exists(outDir)){ dir.create(outDir,recursive = TRUE) }

# This is the filename for our pre-processing of results.
fileName_preProcessing <- paste(workDir,"sourceAnalysis.SHAPE_FILE.r",sep="")
# This is a logical toggle to have the simulations stop if a run does not complete
# NOTE: A run may not be expected to complete if we're simulating something such as evolutionary rescue.
toggle_forceCompletion <- FALSE

# This is a logcial to define if this run is recycled, steps likely ought not to be recycled... only the landscape
# I've not set this up to automatically recycle parameters -- EXCEPT: some NK and RMF model parameters which must be recycled to be equivalent
#### NOTE: If recycling then please view comments at the bottom as the script will resubmit it's next version and there are expectations of how 
####		the "save_batchIndex" line is to have been written.
run_isRecycling <- c("Landscape" = TRUE, "Steps" = FALSE, "Parameters"=TRUE, "Neighbourhood"=FALSE)
# We define what is the max index of our recycled calls
#### NOTE: currReplicate value should be anything but 1 when recycling is not desired.
maxReplicates <- 5
# This is a sanity check, that if we're not recycling anything then... well we don't need to replicate.
if(!any(run_isRecycling)){ assign("maxReplicates",currReplicate,pos=".GlobalEnv") }
# I introduce a recycle_repStart object so that I can track what is the startingIndex value proposed by 
# a recycled script, this allows me to have multiple runs working toward the same end batch of replicates
#### NOTE: The fitnessLandscape will still point to the 1st replicate and thus this must be initialised.....
recycle_repStart <- 1


################################################################################
########################## POPULATION AND SIMULATION ###########################
################################################################################
# We have a focal population value which is interpreted differently based on the growth form used
# if exponential this is a starting value and the value to which the population will be disturbed during perturbations
# for logistic values this is the carrying capacity, for constant this is the population size at all times.
const_focal_popValue <- 4e+06
# The size of the genotype will define the number of sites in which mutations can occur
genomeLength <- 100
# This is the probability of there being a mutation per generation of our individuals (i.e. - mu_g)
const_mutProb <- 1.11e-4
# This is a logical to define if mutations occur only in those newly born, or across the population at large
muts_onlyBirths <- FALSE
# This is a logical toggle which controls if we force rounding of values so that individuals 
# are tracked as whole, ie integer values.
track_asWhole <- FALSE
# We consider the size of the time step to review, this should be a vlue between {0,1} and represents
# the proportion of the life of an individual that passes in a step, this is because: death rate = avg(birth rate) = 1
##### NOTE: This value must follow 0 < size_timeStep <= 1
size_timeStep <- 1
# This is the number of generations that we want to simulate
numGenerations <- 10
# The death probability is a constant and gets passed as the product of our death rate (1) and our timeStep size
const_deathProb <- 0.1
# This is a logical toggle for setting if death rate should be denisty dependent
death_byDensity <- TRUE
# These next two values have no meaning if the death denisty logical toggle is FALSE....
# This is a value for affecting the strength of the the death rates' density dependence (where 1 means linear)
death_densityCorrelation <- 4
# This value determines what is the "carrying cappacity", which is basically the unit for scalling denisty
death_densityCap <- const_focal_popValue
# This is the probability of birth, it is a basal value, makes most sense to be 1, as it is passed along with fitness and size_timeStep
const_birthProb <- 1
# We can set the ancestral fitness, but if we're using a RMF model we should let the distance to optima set this value
# Thus this will actually be reset below....
# NOTE: It is most intuitive that this value be 1 since the initial genotypeFrame uses this value for fitness and not a computed one.
###		However, provided the const_relativeFitness = TRUE, this can be any value, and acts as a baseline to which the distribution's effects are added (outside RMF)
const_ancestFitness <- 1
# We define what are the states possible at each site of our genotype, this ought to be only 0,1 as for the moment I only simulate binary states.
### NOTE: The program handles this for Additive, NK and RMF models only (where it has meaning!)
const_siteStates <- c(0,1)
# This can take several values.  If it's numeric and less than or equal to 1 then is the proportion of total population a population must have reached,
# If this is a value greater than 1 then it is assumed to be the exact size of the lineage and I'll strip any decimal values SO 1.01 will == to any lineage that exists!
# It can also take the form of a type of calculation which has been coded in querryEstablished(), see the function for more, 
# if this is anything else it's assumed to be an expression that can be evaluated in the querryEstablished() function call.
# only affects reporting and tracking not the actual growth dynamics. AND I've not got much crash proofing here, so user beware with expressions.
# Implemented functions: "Desai"
const_estProp <- 1e-4
# This value is the number of individuals a lineage must have before we track its nearest neighbours in our neighbourhood reference database
# This can be any value but I've prefered to use a constant fraction divided by the mutation rate * stepSize so that the fraction controls
# what "probability" a lineage has of generating a mutant in a time step.
#### NOTE: Setting this too low will mean the reference becomes so large you may run out of disk space!  USER BEWARE!
const_hoodThresh <- ceiling(0.005 / (const_mutProb * size_timeStep))


################################################################################
############################## GROWTH FUNCTION #################################
################################################################################
# This is the suite of parameters which affect how population growth occurs, this is a combination of disturbances to size,
# and the form of growth to carrying capacity.

# This is the disturbance type to the population, I've only coded for "bottleneck" or "random",
const_distType <- "bottleneck"
# This sets the initial size of disturbance, which should be considered as the factor by which to reduce the population, to occur in the population
# as well we can include a random component, this only has meaning when the type of disturbance in "random" .. at least for now.
init_distPars <- c("factor"=100,"random"=1)
# This is initialising an object which will be updated to track each time a disturbacne occurs
track_distSize <- NULL
# This is the type of growth function that we'll employ, at present I've only coded for:
# constant (this will overide disturbance mechanics and const_growthForm), exponential (unbounded), logistic growth (bounded) and poisson
const_growthForm <- "logistic"
# We set the basal growth rate of the individuals - this represents the expected number of offspring generated by a single birth event
# this value represents (1+r) as seen in the literature of ODE's dealing with evolution or growth functions.  Another way to see this value
# is as the growth rate of a deterministic system, or the 1+growth rate of a continuous time set of equations.
const_growthRate <- 2
# If a value is set, this is the specific number of generations between disturbance events.  If left as NULL then a value
# will be calculated assuming that growth is logistic and will be based on the disturbance size and growth rate 
const_growthGenerations <- NULL
# This is a logical toggle which when true, the number of births is scaled by the number of deaths.  This is intended for functions 
# other than constant where this is assumed, since deaths are the only source of space for births.  The addition of this function is 
# mostly meant for situations where prob_d > 0 and we're trying to replciate a scenario where the amount of growth is not affected by
# deaths.  Basically implemented to help replicate analytical results when prob_b == prob_d == 1, and since my functions calculate offpsring
# as the difference between current and sequential growth predicted by the growth funcitons, this adjustment becomes necessary.
scaleGrowth_byDeaths <- TRUE
# This is a logical toggle asking if we want drift to be considered or not.  Drift will be simulated by passing the expected number of offpsring
# through a poisson distribuiton.  Otherwise these are deterministic growth functions (except for the stochastic rounding).
includeDrift <- TRUE

# This called function defines the number of steps between disturbances, it can be used in the script to allow fluctuating disturbance times
init_distSteps <- compute_distGrowth(func_distFactor = if(const_growthForm == "logistic"){init_distPars}else{0},
										func_growthType = const_growthForm, 
										func_distType = const_distType, 
										func_growthRate = const_growthRate,
										func_popSize = const_focal_popValue,
										func_focalSize = const_focal_popValue,
										func_manualGenerations = if(const_growthForm == "constant"){ 
																	numGenerations 
																} else {
																	const_growthGenerations
																},
										func_stepDivs = size_timeStep)
# If this is logistic, then we now update the amount of loss, since only now have we calculated what the factor is
### NOTE: unless loss is random, this could have been done earlier	
if(const_growthForm == "logistic"){ init_distSteps["popLost"] <- round(const_focal_popValue/init_distSteps["factor"],0) }

###################################################################################################
######################################## FITNESS LANDSCAPE ########################################
###################################################################################################

# This is the type of landscape model to be used in the simulation, we should describe of types: HoC, Additive, NK, or RMF
# See the fitnessLandscape function for more information
simModel <- "RMF"
# This is the maximum number of mutations that a single new mutant can carry, minimum of 1 for obvious reasons.
## NOTE: This scipt does not support a maxHammign distance of a single mutant being greater than 1, this is implemented for future proffing.
max_numMutations <- 1 #max(1, round(genomeLength * const_mutProb,0))
# This determine if we want to use relative fitness for our initial genotype
# NOTE: If not using relative fitness be certain that the distribution is one centered about zero
#		and such that fitness values drawn are providing values with biologically meaningful selective coefficient (s) values.
const_relativeFitness <- TRUE
# This is a logical of whether or not we allow backmutations to occur in our simulation
allow_backMutations <- TRUE
# This distribution and parameters should represent realistic biological expectations congruent with the mutation supply rate settings
# the following values work well with a skewNorm distribution: c(0.98,0.375,-8,0.01) and gives ~ 2.5% values > 1, with mean ~ 0.686, but mode ~ 0.9
# I also like the values c(0.92,0.32,-10,0.5) for skewNorm, they give a fatter beneficial tail with a good peak near 1, but more total beneficial
# Currently accepted values: Fixed,Gamma,Uniform,Normal,Chi2,beta,exp,evd,rweibull,frechet,skewNorm
constDist <- "exp"
const_distParameters <- c(100)
# This is a logical toggle the user can define so that when calculating fitness, in the fitnessLandscape function
# draws from constDist get -1 (if this is TRUE) so that the raw value is treated as relative fitness which is converter to selection coefficients.
# This is implemented to differentiate between distributions like skewNorm which can simulate the distribution of
# relative fitnesses, from use of exp() which is likely going to be used as selection coefficients.
const_distAsS <- TRUE

# This is the initial Hamming distance that our ancester should be from the theoretical global phenotype
# NOTE: This value only has meaning if the model is RMF
const_RMF_initiDistance <- 5
# From Neidhart 2014, theta ~ 0.25 was on the high order of ruggedness according to empirical fitness landscapes.
# Thus we calculate our "c" value (or independent weight) by using our defined theta, and the distribution from which we'll be drawing our random component
const_RMF_theta <- 0.35

# This is the K value for an NK landscape, it has no meaning otherwise.
const_numInteractions <- 4

# This object has no meaning unless the fitness landscape model is "Fixed", in which case this is the file name for something in the wordDir  
# for a .csv file which can be read as a data.frame containing two columns named "binaryString", "fitness".
const_fixedFrame <- NULL

###################################################################################################
################################## END OF CONSTANT DEFINITIONS ####################################
###################################################################################################