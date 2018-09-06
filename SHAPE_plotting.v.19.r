# This is a script to review a suite of SHAPE output and then compile the information stored therein

############ DEPENDENCIES ############
# This is based on the SHAPE_body.v.28.r, and sourceAnalysis.v.14.r output and associated dependencies 

############ ASSUMPTIONS #############
# This assumes that all the files you want analysed are held within the workDir or some subfolder thereof.

# Clear our workspace and any previously loaded session
rm(list=ls())

# load libraries
library(RColorBrewer)
library(RSQLite)
library(plot3D)
library(colorRamps)
library(lattice) # This is for heat map plotting.
library(abind) # This allows me to bind matrices/data.frames into arrays 
library(foreach)
# We register our parallel backend 
numCores = 4
if(.Platform$OS.type == "windows"){
	library(doSNOW)
	registerDoSNOW(makeCluster(numCores, type="SOCK"))
} else {
	library(doMC)
	registerDoMC(cores=numCores)
}

# This option will make certain that my data.frames are created as strings and not factors
options(stringsAsFactors = FALSE)

###################################################################################################
#################################### BEGIN OF PREDEFINITIONS ######################################
###################################################################################################
# Is this run locally or on a server?  This is different than the main run's serverFarm hence why it's defined here.
plot_onServer <- TRUE

# This is a control of which types of plots you'd like made:
plotTypes <- c("Muller"=FALSE,
				"walkLength"=TRUE,
				"popDemographics" = TRUE,
				"repeatability"=TRUE,
				"landscapeCompare"=TRUE,
				"entropy"=TRUE)

####### FILE AND SYSTEM LOCATIONS #######
# These are the basic and working directories in which we can find all the information of interest.
baseDir <- if(plot_onServer){ "/global/home/hpc3058/jonathan/"	} else if(!plot_onServer) { "E:/Research_Data/" }
workDir <- paste(baseDir,"SimulationRuns/multiMutants/localRun/SHAPETest/",sep="")
# This is where the output of our summary should be placed
outDir <- paste(workDir,"postAnal/",sep="")
# This creates the directories if they does not exist
if(!dir.exists(outDir)){ dir.create(outDir,recursive = TRUE) }

# We load the general functions for this script
fileName_functions <- "sourceFunctions.v.15.r"
if(file.exists(paste(workDir,fileName_functions,sep=""))){
	source(paste(workDir,fileName_functions,sep=""))
} else {
	stop("Could not find the functions source, please review.")
	Sys.sleep(10)
	q(save="no")
}

# These are file naming convention objects that, while stored in any run's parameters, may not be in the environment yet without this.
save_batchBase <- "testJobs"
save_batchSet <- TRUE
save_batchJob <- TRUE
# This is the string used to separate items
sepString <- "_"

# This is a regular expression which should define how directories will be named, and save_batchBases were built
save_jobExpression <- name_batchString(funcBase = save_batchBase,
										func_setID = if(save_batchSet){paste('[^',sepString,']+',sep="")}else{NULL},
										func_jobID = if(save_batchJob){paste('[^',sepString,']+',sep="")}else{NULL},
										func_sepString = sepString)

# These are the file suffixes, 
saveSuffixes <- c("dataBase"=".sqlite",
					"run"=".RData")
# This is the pattern of information we want to have and work with
usePattern <- if(!is.null(save_batchBase)){ save_batchBase } else {".RData"}
# These are filenames for information of collected batches of data
imageFilename <- c("saveDir"=outDir,
					"fileList"=paste("allFiles_from_", save_batchBase,".RData",sep=""),
					"parameters"=paste("jobParameters_from_", save_batchBase,".RData",sep=""),
					"walkPlot"=paste("walkPop_from_",save_batchBase,".RData",sep=""),
					"popDemographics"=paste("popDemographics_from_", save_batchBase,".RData",sep=""),
					"plotPop"=paste("plotPop_from_",save_batchBase,".RData",sep=""),
					"repeatability"=paste("repeatabilityData_from_", save_batchBase,".RData",sep=""),
					"plotRepeat"=paste("plotRepeat_", save_batchBase,".RData",sep=""),
					"plotLandscape"=paste("plotLandscape_", save_batchBase,".RData",sep=""),
					"plotEntropy"=paste("plotEntropy_", save_batchBase,".RData",sep=""))
					
processedData_filePattern <- paste("processed_runData_from_", save_batchBase,sep="")
processedObjects <- c("runDemographics","info_estLines")
# This is a string pattern which will help me define that an object is an envrionment
envString <- "envString"
# These are strings for naming our processed objects into smaller chunks
namingString <- c("popDemo"="popInfo_obj_",
					"repeatability"="repInfo_obj_")

############### Constants #################
# These are colours and plotting points characters that we will cycle through
# We will use the unqiue colours of the non-paired Qualitative sets from RColorBrewer
allCols <- unique(c(brewer.pal(9,name="Set1"),
					brewer.pal(8,name="Set2"),
					brewer.pal(9,name="Set3"),
					brewer.pal(8,name="Dark2"),
					brewer.pal(8,name="Accent")))

allPoints <- c(21:25)
# Unique point and colour combinations will be used to separate our different run parameter types as best as possible
pchCol_allCombos <- expand.grid("col"=allCols,"pch"=allPoints, stringsAsFactors = FALSE)
# This is the colour transparency to be used with plotting.
plot_colourAlpha <- 0.7
plot_lightAlpha <- 0.5
plot_jitterBuffer <- 0.25
# Where phi is the up-down angle and theta is the roation
plot3D_angles <- c("phi"=15,
					"theta"=60)

# Here are some constant string types being used by my pot analysis:
string_lineDescent <- "_->_"
string_sepLines <- "__and__"
# This is the number of decimal places you want reported values to be rounded to
const_sigFig <- 4


###################################################################################################
##################################### END OF PREDEFINITIONS #######################################
###################################################################################################


###################################################################################################
##################################### BEGIN OF FUNCTIONS ##########################################
###################################################################################################
# This allows the opacity, saturation and hue of colours to be adjusted.  The purpose is to allow
# a range of basic named colours to be dynamically manipulted into altered vectors of the colour. 
colTransparency <- function(func_cols, func_opacity=1, func_scaleSaturation = NA, func_scaleValue = NA){
	# In case the user has passed scaling values outside of the 0-1 range we simply return them to being values of 1
	funcParms = list("Saturation"= func_scaleSaturation,
					"Value"= func_scaleValue,
					"Opacity"= func_opacity)
	# We check that the opacity value is within the range of 0-1				
	for(thisParm in names(funcParms)){
		# If the user has not defined anything, then let it pass....
		if(all(!is.na(funcParms[[thisParm]]))){
			if(any(funcParms[[thisParm]] > 1)){
				funcParms[[thisParm]][which(funcParms[thisParm] > 1)] <- 1
			} else if (any(funcParms[[thisParm]] < 0)){
				funcParms[[thisParm]][which(funcParms[thisParm] < 0)] <- 0
			}
		}
	}
	
	# After extracting the RGB balues of a colour we change the values based on the opacity.
  	tmpReturn <- rgb2hsv(sapply(func_cols, col2rgb))
  	return( sapply(1:ncol(tmpReturn),function(x){
	  				hsv(h = tmpReturn[1,x], 
	  					s = if(any(is.na(funcParms[["Saturation"]]))){
	  							tmpReturn[2,x]
	  						} else {
	  							funcParms[["Saturation"]]
	  						}, 
	  					v = if(any(is.na(funcParms[["Value"]]))){
	  							tmpReturn[3,x]
	  						} else {
	  							funcParms[["Value"]]
	  						}, 
	  					alpha= funcParms[["Opacity"]]) 
	  			}) )
}

# This is a function draw a number of jitter values from a range of values
definedJitter <- function(func_returnNum, func_inRange=NULL){
	return( runif(func_returnNum,
					min=if(!is.null(func_inRange)){min(func_inRange)}else{0},
					max=if(!is.null(func_inRange)){max(func_inRange)}else{1}) )
	
}

# This is a function for standardised naming conventions
nameObject <- function(func_inString, func_inPrefix, func_splitStr = FALSE){
	if(func_splitStr){
		unlist(lapply(strsplit(func_inString, func_inPrefix),function(x){ x[length(x)] }))
	} else {
		paste(func_inPrefix, func_inString,sep="")
	}
}

# This is a function for performing ceiling to a value at a particular number of decimal places
decimalCeiling <- function(func_value, func_level=0){ 
	# Frist step is to find out how many decimal places exist within this value
	func_value <- strsplit(as.character(func_value),'.',fixed=TRUE)[[1]]
	# If there are no values after the decimal we simply return it, otherwise we find the number
	# of places after the decimal, and then round as appropriate to the func_level
	if(length(func_value) > 1){
		# We split apart the second value into individual characters and make them into elements concatenated with the integers
		func_deciList <- c(as.numeric(func_value[1]),as.numeric(strsplit(func_value[2],"")[[1]]))
		# we now increase the level value by one since our vector includes one additional element
		func_level <- func_level + 1
		# If the user is trying to ceiling at a level below the current 
		# number of decimal values contained then we perform no operation!
		if(length(func_deciList) > func_level){
			# We check if there is anything to be rounded, then do so
			if(any(func_deciList[func_level : length(func_deciList)] > 0)){
				func_deciList[func_level] <- func_deciList[func_level] + 1
				if(func_deciList[func_level] >= 10){
					func_deciList[func_level] = func_deciList[func_level] - 10
					func_deciList[func_level - 1] = func_deciList[func_level - 1] +1
				}
			}
			# we now return the value ceilinged to our level
			return( as.numeric(paste(func_deciList[1],'.',paste(func_deciList[2:func_level],collapse=""),sep="")) )
		# If we're trying to round at a number of decimal places smaller than the number of decimal values that exist
		# then we simply return the value as given
		} else {
			return( as.numeric(paste(func_value,sep='',collapse='.')) )
		}
	# If there are noa decimal values then we simply return that first portion 	
	} else {
		return( as.numeric(func_value) )	
	}
}

# This function takes a defined amount of space nad divides it into an equal number of sub-spaces based on the number of sub-groups defined
# as well as the amount of space between final groupings.  The subDivs should be a matrix with the sets of sequentially nested subdivisions 
# we want to define space for the tmpScale is the total area which can be allocated and the separator between final units of space 
# This function will return a list of lists for the ranges at each level.  This function resolves lists by sequentially calling itself through all list items
# The tmpSpace is used to define what range of values can be permitted, tmpBuffer is used to define the proportion of given space to use as trim/ or not use....
# The last argument is to querry if we want to order the values being subdivided.
define_plotSpace <- function(func_subDivs, func_tmpSpace = c(0,1),func_tmpBuffer = 0.1, func_tmpBuffer_max = 0.4, func_tmpBuffer_min=0.2, func_orderDivs = FALSE){
	# We define what are the divisions of space that we want to pass to each part, we add one so that this forces the division points to be boundaries
	# between our spaces to be defined for the current func_tmpSpace
	func_tmpDivs <- if(is.matrix(func_subDivs) || is.data.frame(func_subDivs)){
						seq(func_tmpSpace[1], func_tmpSpace[2],length.out = length(unique(func_subDivs[,1]))+1)
					} else {
						seq(func_tmpSpace[1], func_tmpSpace[2],length.out = length(unique(func_subDivs))+1)
					}
	# Now to define the actual space permissable we create space sets for all but the first and last elements, which are the boundary points.
	func_tmpDivs <- lapply(1:(length(func_tmpDivs)-1),function(func_thisSpace){ 
								# We find the space that could be occupied
								func_tmpSpace <- c(func_tmpDivs[func_thisSpace],func_tmpDivs[func_thisSpace+1])
								# So the actual space we'll return as permissible will be these boundaries, with trim of total space
								return( c(func_tmpSpace[1] + max(func_tmpBuffer_min, min(func_tmpBuffer_max,(func_tmpBuffer * (func_tmpSpace[2]-func_tmpSpace[1])))), 
											func_tmpSpace[2] - max(func_tmpBuffer_min, min(func_tmpBuffer_max,(func_tmpBuffer * (func_tmpSpace[2]-func_tmpSpace[1]))))) )
							})	
	# Using the scale of space we've been given, we'll divide it by the nunber of unique instances in the first
	# column of func_subDivs (or if it's a vector simply by that value alone)
	if(is.matrix(func_subDivs) || is.data.frame(func_subDivs)){
		# If we've been asked to order our subDivs elements, we do so
		if(func_orderDivs){ func_subDivs <- func_subDivs[order(func_subDivs[,1]),] }
		# We also define what are the rows of the func_subDivs to be passed in each sub division
		func_tmpParts <- lapply(unique(func_subDivs[,1]),function(func_thisPart){ which(func_subDivs[,1] == func_thisPart) })
		# Now the we've defined what divisions of space and parts/rows of the current data need to be subdivided we iteratively call this forward.
		func_tmpReturn <- lapply(1:length(unique(func_subDivs[,1])),function(func_thisDiv){ 
									define_plotSpace(func_subDivs = func_subDivs[func_tmpParts[[func_thisDiv]],-1], 
													func_tmpSpace  = func_tmpDivs[[func_thisDiv]],
													func_tmpBuffer = func_tmpBuffer)
								})
		return( func_tmpReturn )
	} else {
		# This mean's we're on the last element of the divisions and we can return the divisions that we've defined
		return( func_tmpDivs )
	}
}

# This calculates the entropy equation similar to Szendro 2013 - it is to calculate the degree of determinism
calcEntropy <- function(func_probVec){
	if(any(func_probVec < 0)){ stop("There was some negative numbers passed to calcEntropy...") }
	func_probVec <- func_probVec[which(func_probVec > 0)]
	func_probVec <- func_probVec
	return( -sum(func_probVec * log(func_probVec, base = exp(1))) )
}

# This quick little function is a means for me to create the strings of environments and subsequently extract information back out
# This is to be used for the names of our step report object tables
nameEnviron <- function(func_Index, funcSplit = FALSE, funcBase = envString){
	# This checks if we're asking for the name to be split or not
	if(funcSplit){
		# We return only the second piece of information given the setup of this function's naming practice.
		return( unlist(lapply(strsplit(func_Index,"_e_"),function(x){ x[2] })) )
	} else {
		return( paste(funcBase,func_Index,sep="_e_") )
	}
}

###################################################################################################
###################################### END OF FUNCTIONS ###########################################
###################################################################################################



###################################################################################################
##################################### BEGIN OF DATA LOADING #######################################
###################################################################################################


######################################################################
############### STEP 1: What files exist? ############################
######################################################################
# First step we identify the jobs that are to be worked with, this should save two objects:
# all_proccessedFiles, all_dividedFiles
if(file.exists(paste(imageFilename["saveDir"],imageFilename["fileList"],sep=""))){
	load(paste(imageFilename["saveDir"],imageFilename["fileList"],sep=""))
} else {
	all_proccessedFiles <- list.files(path=workDir,pattern= processedData_filePattern, recursive = TRUE, full.names=TRUE)
	# Now I know what the epxression is for jobs, but better yet I know how the jobs were built and so can extract all my info
	all_jobInfo <- name_batchString(funcBase = paste(save_batchBase, 
												unlist(lapply(strsplit(all_proccessedFiles, processedData_filePattern),function(thisFile){ 
																gsub(saveSuffixes["run"],"",thisFile[2])
															})),
												sep=""),
									func_setID = save_batchSet,
									func_jobID = save_batchJob,
									func_repID = TRUE,
									funcSplit = TRUE,
									func_sepString = sepString)
	# We now find unique jobs by looking at the information just after the processedData_filePattern
	# The regular expressions are so that I can have subsetted replicates of jobs yet they be considered part of a single set.
	all_uniqueJobs <- name_batchString(funcBase = save_batchBase,
										func_setID = paste('[^',sepString,']+',sep=""),
										func_jobID = unique(all_jobInfo["jobID",]),
										func_sepString = sepString)
	# We now subdivide our all_proccessedFiles into the different unique jobs
	all_dividedFiles <- sapply(all_uniqueJobs,function(x){
									return( all_proccessedFiles[which(grepl(paste(x,"(.+)",sep=sepString), all_proccessedFiles))] )
						},simplify=FALSE)
	# As a sanity check we ensure that all the all_proccessedFiles have been placed in our dividedfiles object
	if(length(all_proccessedFiles) != length(unlist(all_dividedFiles))){
		print(" There was a problem when dividing the all_proccessedFiles into divdedFiles, please review")
		Sys.sleep(10)
		q(save="no")
	}
	save(all_proccessedFiles, all_jobInfo, all_dividedFiles,file=paste(imageFilename["saveDir"],imageFilename["fileList"],sep=""))
}


######################################################################
############### STEP 2: Build some analysis specific data sets #######
######################################################################


########################## Parameters ##########################
# This is the set of information for repeatability of evolution
if(file.exists(paste(imageFilename["saveDir"],imageFilename["parameters"],sep=""))){
	load(paste(imageFilename["saveDir"],imageFilename["parameters"],sep=""))
# Then we need to find all the processedData_filePattern objects, create an environment for them,
# then load parameters and the processed information.	
} else {
	# We find all the directories which use the names of our divided files
	tmpDirs <- list.dirs(path=workDir,full.names=FALSE,recursive=FALSE)
	tmpDirs <- unique(unlist(lapply(names(all_dividedFiles), function(thisName){
						return( tmpDirs[which(grepl(thisName,tmpDirs))] )
				})))
		
	# I need to create an outer list object which stores the information for each job
	all_parmInfo <- NULL
	# Now we break the processing into chunks to be added up and collected
	tmp_chunkSize <- 100
	tmp_numChunks <- ceiling(length(tmpDirs)/tmp_chunkSize)
	for(thisChunk in 1:tmp_numChunks){
		tmpAdd <- NULL
		for(thisJob in tmpDirs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(tmpDirs))]){
			# We load the parameters for thisJob
			tmp_parmFile <- list.files(path=paste(workDir,thisJob,"/",sep=""), pattern = "Parameters")
			if(length(tmp_parmFile) == 1){
				assign("theseParms", new.env(), pos=".GlobalEnv")
				load(paste(workDir,thisJob,"/",tmp_parmFile,sep=""),envir=eval(as.name("theseParms")))
			} else {
				stop(paste("Could not find the parameters files to load for ",thisJob,sep=""))
				Sys.sleep(10)
				q(save="no")
			}
			# I now build an object which captures the runParms of interest, this should be a list object
			tmp_parmsRef <- c(theseParms$runParameters$Population,
								theseParms$runParameters$Growth_Disturbance[-c(2,3,4)],
								"distFactor"= unname(theseParms$runParameters$Growth_Disturbance$init_distPars["factor"]),
								"distSpread"= unname(theseParms$runParameters$Growth_Disturbance$init_distPars["random"]),
								"mean_realisedDilution"=mean(theseParms$runParameters$Growth_Disturbance$track_distSize[,"factor"]),
								theseParms$runParameters$FitnessLandscape,
								theseParms$runParameters$DFE,
								"db_splitTables"= theseParms$runParameters$DataManagement$db_splitTables)
			tmp_parmsRef$const_distParameters <- paste(tmp_parmsRef$const_distParameters,collapse=",")
			
			# We now assign this information into our list object, the jobs set looks for which of the names(all_dividedFiles)
			# associates with thisJob for latter grouping of data.  To assign a unique jobSet we look if there is the jobString
			# otherwise use the directory name.
			tmpReturn <- data.frame("jobName"=thisJob,
									"jobSet"= if(save_batchSet){
												# We then check which of the all_dividedFiles names could refer to this
												# We then choose the longest one that matches, as that will be the most complete
												tmpNames <- names(all_dividedFiles)[sapply(names(all_dividedFiles),grepl,"x"=thisJob)]
												tmpNames[which.max(nchar(tmpNames))]
											} else {
												thisJob
											},
									t(unlist(tmp_parmsRef)))
			# We now just reset the data types
			for(thisCol in names(tmp_parmsRef)){ mode(tmpReturn[,thisCol]) <- mode(tmp_parmsRef[[thisCol]]) } 					
			tmpAdd <- rbind(tmpAdd, tmpReturn)
		}
		all_parmInfo <- rbind(all_parmInfo, tmpAdd)
	}
	# We save the gathered information
	save(all_parmInfo,file=paste(imageFilename["saveDir"],imageFilename["parameters"],sep=""))
} # this closes out gathering the informaion for reapeatability


########################## popDemographics ##########################
# This is the set of information for repeatability of evolution
if(!file.exists(paste(imageFilename["saveDir"],imageFilename["popDemographics"],sep=""))){
	print("Starting popDemographics")
	# We create an object that tracks the names of save files, objects and jobs for each set
	all_popSets <- sapply(unique(all_parmInfo$jobSet), function(thisSet){
							tmp_thisName <- nameObject(func_inString = thisSet,
																func_inPrefix = namingString["popDemo"])
							return( list("objName"=tmp_thisName,
										"saveFile"= paste(sub(save_jobExpression,"",sub(save_batchBase,tmp_thisName,imageFilename["popDemographics"]),fixed=TRUE),sep=""),
											"setJobs" = unique(all_parmInfo$jobName[which(all_parmInfo$jobSet == thisSet)])) )
								
						}, simplify = FALSE)
	# For all_uniqueSets, we find the ones where the files don't exist so those are the only we build
	tmp_missingSets <- names(all_popSets)[unlist(lapply(all_popSets,function(thisSet){ !file.exists(paste(imageFilename["saveDir"],thisSet[["saveFile"]],sep="")) }))]
	# Now for each main job we'll gather some information and save it into it's own list object
	tmpCheck <- foreach(thisSet = tmp_missingSets,.combine="c") %dopar% {
							tmpList <- vector(mode="list",length=length(all_popSets[[thisSet]][["setJobs"]]))
							names(tmpList) <- all_popSets[[thisSet]][["setJobs"]]
							# Ok here is the level at which we'll place our foreach looping to create a list of results
							for(thisJob in all_popSets[[thisSet]][["setJobs"]]){
								# We find the job files associated to thisJob
								tmpJobs <- all_dividedFiles[[thisSet]][which(grepl(thisJob, all_dividedFiles[[thisSet]]))]
								if(length(tmpJobs) > 0){
									# I need to create an outer list object which stores the information for each job
									tmp_replicateInfo <- NULL
									# Now we break the processing into chunks to be added up and collected
									tmp_chunkSize <- 100
									tmp_numChunks <- ceiling(length(tmpJobs)/tmp_chunkSize)
									for(thisChunk in 1:tmp_numChunks){
										tmpAdd <- NULL
										# Now for each and every replicate within this job we want to grab the:
										# (1) the population wide - through time - min, mean,max fitness. - I'll return the whole demoMat
										for(thisFile in tmpJobs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(tmpJobs))]){
											# We load this file's data
											tmpString <- nameEnviron(strsplit(strsplit(thisFile,"/")[[1]][length(strsplit(thisFile,"/")[[1]])],".",fixed=TRUE)[[1]][1],
															funcBase = envString)
											assign(tmpString, new.env(), pos = ".GlobalEnv")
											load(thisFile, envir = eval(as.name(tmpString)) )
											# Ok we'll now return this information along with the demoMat
											tmpAdd <- c(tmpAdd, list("popDemographics"= eval(as.name(tmpString))$runDemographics$demoMat))
											rm(list=tmpString, pos=".GlobalEnv")
											# Clear up memory
											gc() 
										} # This closes out the for loop gathering information about the replicates for this job
										tmp_replicateInfo <- c(tmp_replicateInfo, tmpAdd)
										rm(tmpAdd)
										gc()
									}
									
									#tmp_all_popDemos <- array(sapply(tmp_replicateInfo,function(x){ x[["popDemographics"]] }),
									tmp_all_popDemos <- array(sapply(tmp_replicateInfo,function(x){ x }),
																dim=c(dim(tmp_replicateInfo[[1]]),length(tmp_replicateInfo)),
																dimnames = list(rownames(tmp_replicateInfo[[1]]),
																				colnames(tmp_replicateInfo[[1]]),
																				NULL))
									# For the population demograpihcs we collapse this into the mean values through time across all replicates
									tmp_all_popDemos <- matrix(apply(tmp_all_popDemos,MARGIN=2,function(thisCol){ 
																	# Now we return the mean and sd values for each column
																	c(apply(thisCol,MARGIN=1,mean),apply(thisCol,MARGIN=1,sd))
																			
																}),nrow=nrow(tmp_all_popDemos),ncol=ncol(tmp_all_popDemos)*2,
																dimnames=list(rownames(tmp_all_popDemos),
																				unlist(lapply(colnames(tmp_all_popDemos),function(x){ return(c(x,paste("sd",x,sep="_")))}))) ) 
									# We now assign this information into our list object
									tmpList[[thisJob]] <- tmp_all_popDemos
									rm(tmp_all_popDemos, tmp_replicateInfo)
									gc()									
								} # This closes out the logical conditional that there is at least some replicates in thisDir
							# This closes out the loop for the different jobs
							}
							# We can now save this object to a file
							assign(all_popSets[[thisSet]][["objName"]], tmpList, pos=".GlobalEnv")
							save(list = all_popSets[[thisSet]][["objName"]],file= paste(imageFilename["saveDir"],all_popSets[[thisSet]][["saveFile"]],sep=""))
							rm(list= all_popSets[[thisSet]][["objName"]])
							rm(tmpList)
							# Clear up memory
							gc()
							# We return a confirmation that the job completed 
							return(  nameObject(func_inString = all_popSets[[thisSet]][["objName"]],
												func_inPrefix = namingString["popDemo"],
												func_splitStr = TRUE)  )
						}
	# We now check if all the jobs completed
	tmp_checkCompleted <- sapply(tmp_missingSets,is.element,"set"=tmpCheck)
	if(all(tmp_checkCompleted)){
		# We now save the object which stores all the population demographics data object names and files locations
		save(all_popSets,file=paste(imageFilename["saveDir"],imageFilename["popDemographics"],sep=""))
	} else {
		stop(paste("There was a problem with the sets of: ",paste(names(tmp_checkCompleted)[which(!tmp_checkCompleted)],collapse=" ",sep=" "),sep=""))
	}
# this closes out gathering the informaion for popDemographics
} else {
	load(paste(imageFilename["saveDir"],imageFilename["popDemographics"],sep=""))
}



########################## Repeatability ##########################
# This is the set of information for repeatability of evolution
if(!file.exists(paste(imageFilename["saveDir"],imageFilename["repeatability"],sep=""))){
	print("Starting repeatability")
	# We create an object that tracks the names of save files, objects and jobs for each set
	all_repSets <- sapply(unique(all_parmInfo$jobSet), function(thisSet){
							tmp_thisName <- nameObject(func_inString = thisSet,
														func_inPrefix = namingString["repeatability"])
							return( list("objName"=tmp_thisName,
										"saveFile"= paste(sub(save_jobExpression,"",sub(save_batchBase,tmp_thisName,imageFilename["repeatability"]),fixed=TRUE),sep=""),
											"setJobs" = unique(all_parmInfo$jobName[which(all_parmInfo$jobSet == thisSet)])) )
								
						}, simplify = FALSE)
	# For all_uniqueSets, we find the ones where the files don't exist so those are the only we build
	tmp_missingSets <- names(all_repSets)[unlist(lapply(all_repSets,function(thisSet){ !file.exists(paste(imageFilename["saveDir"],thisSet[["saveFile"]],sep="")) }))]
	# Now for each main job we'll gather some information and save it into it's own list object
	tmpCheck <- foreach(thisSet = tmp_missingSets,.combine="c") %dopar% {
							library(RSQLite)
							tmpList <- vector(mode="list",length=length(all_repSets[[thisSet]][["setJobs"]]))
							names(tmpList) <- all_repSets[[thisSet]][["setJobs"]]
							# Ok here is the level at which we'll place our foreach looping to create a list of results
							for(thisJob in all_repSets[[thisSet]][["setJobs"]]){
								# We load this job's data, it'll be universal for all the sub-divided files, as per the current recycling at least...
								# It stays in the lapply call because of the file specificity.
								tmp_jobString <- nameEnviron(thisJob,
															funcBase = envString)
								assign(tmp_jobString, new.env(), pos = ".GlobalEnv")
								# We also create a parameter space for this job
								load(paste(workDir, thisJob,"/",thisJob,"_Parameters_1.RData",sep=""), envir = eval(as.name(tmp_jobString)))
								
								assign("string_tableNames",eval(as.name(tmp_jobString))$runParameters$DataManagement$string_tableNames, pos=".GlobalEnv")
											
								# This builds a connection required for fixing the absRank measurements
								connections_dataBase <- sapply(names(eval(as.name(tmp_jobString))$runParameters$DataManagement$fileName_dataBase),function(x){ 
																	resetDatabase(paste(workDir,thisJob,"/", eval(as.name(tmp_jobString))$runParameters$DataManagement$fileName_dataBase[[x]],sep=""),
																					func_type = "connect") })
								func_hoodTables <- dbListTables(connections_dataBase$nearestNeighbours)
								func_landscapeCon <- connections_dataBase$genotypeSpace
								func_landTables <- dbListTables(func_landscapeCon)
								allow_backMutations <- eval(as.name(tmp_jobString))$runParameters$FitnessLandscape$allow_backMutations
								max_numMutations <- eval(as.name(tmp_jobString))$runParameters$Population$max_numMutations
								genomeLength <- eval(as.name(tmp_jobString))$runParameters$Population$genomeLength
								db_splitTables <- eval(as.name(tmp_jobString))$runParameters$DataManagement$db_splitTables
													
								# We find the job files associated to thisJob
								tmpJobs <- all_dividedFiles[[thisSet]][which(grepl(thisJob, all_dividedFiles[[thisSet]]))]
								# I now build a list for storing the information of our replicate
								tmp_replicateInfo <- NULL
								if(length(tmpJobs) > 0){
									# Now for each and every replicate within this job we want to grab the information concerning transitions and 
									# the fitness landscape of the parent to child.  We don't define .combine to gets list return
									# I need to create an outer list object which stores the information for each job
									tmp_replicateInfo <- NULL
									# Now we break the processing into chunks to be added up and collected
									tmp_chunkSize <- 100
									tmp_numChunks <- ceiling(length(tmpJobs)/tmp_chunkSize)
									for(thisChunk in 1:tmp_numChunks){
										tmpAdd <- NULL
										# Now for each and every replicate within this job we want to grab the:
										# (1) the population wide - through time - min, mean,max fitness. - I'll return the whole demoMat
										for(thisFile in tmpJobs[(1+((thisChunk - 1)*tmp_chunkSize)):min((thisChunk*tmp_chunkSize),length(tmpJobs))]){
											tmp_fileString <- nameEnviron(strsplit(strsplit(thisFile,"/")[[1]][length(strsplit(thisFile,"/")[[1]])],".",fixed=TRUE)[[1]][1],
															funcBase = envString)
											assign(tmp_fileString, new.env(), pos = ".GlobalEnv")
											load(thisFile, envir = eval(as.name(tmp_fileString)) )
											# We extract what were the established lineages 
											tmpEstablished <- eval(as.name(tmp_fileString))$runDemographics$vec_estLineages
											# the final lineages which existed can be found by checking the last row of info_estLines[["lineDemo"]][,,as.character(tmpEstablished)]
											lineDemo_maxStep <- which.max(as.numeric(nameTable_step(rownames(eval(as.name(tmp_fileString))$info_estLines$lineDemo),funcSplit = TRUE)))
											
											# We extract what were the established lineages 
											tmp_final_estLineage <- data.frame("genotypeID"=tmpEstablished[which(eval(as.name(tmp_fileString))$info_estLines$lineDemo[lineDemo_maxStep,"isEstablished",as.character(tmpEstablished)] == 1)],
																				"binaryString"="")
											tmp_final_estLineage[,"binaryString"] <- unlist(lapply(tmp_final_estLineage[,"genotypeID"],function(x){
																										retrieve_binaryString(func_genotypeID = x,
																																func_landscapeCon = func_landscapeCon,
																																func_subNaming = db_splitTables)[,"binaryString"]
																									}))
											# The final dominant lineage is the transitioned lineage which had the largest population on the last step.  In the event of a tie
											# we look for the lineage which transitioned latest and break subsequent ties by highest fitness.  All to ensure a single return.
											tmp_transitionMat <- eval(as.name(tmp_fileString))$runDemographics$transitionMat
											tmp_final_domLineage <- cbind("popSize"= eval(as.name(tmp_fileString))$info_estLines$lineDemo[lineDemo_maxStep,"popSize",as.character(tmp_transitionMat[,"genotypeID"])],
																			tmp_transitionMat[,c("Step","fitness","genotypeID")])
											tmp_maxSized = which(tmp_final_domLineage[,"popSize"] == max(tmp_final_domLineage[,"popSize"]))
											# As a final safety, but unlikely requirement, if multiple genotypeID's transitioned at the same time, have the same popSize,
											# and have the same fitness, we simply take the first value.  This should never matter.
											tmp_final_domLineage <- if(length(tmp_maxSized) != 1){
																		tmpReturn <- which(tmp_final_domLineage[tmp_maxSized,"Step"]==max(tmp_final_domLineage[tmp_maxSized,"Step"]))
																		if(length(tmpReturn) > 1){
																			tmpReturn <- tmpReturn[which(tmp_final_domLineage[tmp_maxSized[tmpReturn],"fitness"]==max(tmp_final_domLineage[tmp_maxSized[tmpReturn],"fitness"]))][1]
																		}
																		tmp_final_domLineage[tmp_maxSized[tmpReturn],"genotypeID"]
																	} else {
																		tmp_final_domLineage[tmp_maxSized,"genotypeID"]
																	}
											# We now build this information into an object shaped as intended for use.
											tmp_final_domLineage <- c("genotypeID"=tmp_final_domLineage,
																		"binaryString"="")
											tmp_final_domLineage["binaryString"] <- tmp_final_estLineage[which(tmp_final_estLineage[,"genotypeID"] == tmp_final_domLineage["genotypeID"]),"binaryString"]
											# We now also want to track the line of descent for the dominant lineage and it's first transition
											tmp_line_of_descent <- as.character(eval(as.name(tmp_fileString))$info_estLines$end_Lines_of_Descent[[tmp_final_domLineage["genotypeID"]]])
											# We now try and extract the first transition in this line of descent
											tmp_dom_transitions <- strsplit(tmp_line_of_descent,string_lineDescent)[[1]]
											# Now if there was no transition we should be left with the WT genotype ID of "0" as the onyl element... otherwise we have a first step.
											tmp_firstTransition <- if(length(tmp_dom_transitions) == 1){
												paste(rep(tmp_dom_transitions[1],2),collapse=string_lineDescent)
											} else {
												paste(tmp_dom_transitions[1:2],collapse=string_lineDescent)
											}
											# We now add more information to the domLineage vector
											tmp_final_domLineage <- c(tmp_final_domLineage,
																	"line_ofDescent"= paste(tmp_line_of_descent,collapse=string_sepLines),
																	"firstStep"=tmp_firstTransition)
											
											# We build now the transition mat into something which shows the steps taken on the fitness landscape
											# In order to know the transitions, we look at the transition matrix and for a lineage, but starting with the first
											# that is found in the tranisition matrix which was not established on Step_0
											tmp_initialEst <- dimnames(eval(as.name(tmp_fileString))$info_estLines$lineDemo)[[3]]
											tmp_initialEst <- tmp_initialEst[which(eval(as.name(tmp_fileString))$info_estLines$lineDemo[nameTable_step(1),"isEstablished", tmp_initialEst] == 1)]
											# If there are any genotypes which were not in the initial set we go further, otherwise we're returning some null stuff
											tmp_returnMat <- NULL
											if(any(!is.element(tmp_transitionMat[,"genotypeID"],as.numeric(tmp_initialEst))) && length(tmp_initialEst) > 0){
												# Ok this means that some transition occured.  I now just want to find what transitions occured
												tmp_landscapeTopology <- eval(as.name(tmp_fileString))$info_estLines$landscapeTopology
												tmpSteps <- strsplit(rownames(tmp_landscapeTopology),string_lineDescent)
												# This is hte first row of our transitionMat at which there was a change
												tmp_minRow <- min(which(!is.element(tmp_transitionMat[,"genotypeID"],as.numeric(tmp_initialEst))))
												# Ok this is what will be assigned as the return matrix, it should inform on the
												# genotypeID_1, fitness_1, absRank, numMuts_1, transition, transitionStep
												for(thisRow in tmp_minRow:nrow(tmp_transitionMat)){
													tmpID <- tmp_transitionMat[thisRow, ]
													# We look to the landscapeTopology object to get transition information, what we are looking for 
													# is an instance where the second element (offspring) is the same as our tmpID, and the progenitor
													# is some element of the pedigree for our lineage.
													tmp_possProgenitors <- which(sapply(tmpSteps,function(x){ 
																						as.numeric(x[2]) == tmpID["genotypeID"] 
																					}))
													# If this transition is simply due to a shift in dominance, rather than a new type whic has emerged,
													# due to mutation, then there may be no progenitors.  In this case we skip it as it is not a mutational step.
													if (length(tmp_possProgenitors) == 0){ next }
													tmp_possProgenitors <- matrix(tmp_landscapeTopology[tmp_possProgenitors,],
																					nrow = length(tmp_possProgenitors),
																					dimnames=list(rownames(tmp_landscapeTopology)[tmp_possProgenitors],
																									colnames(tmp_landscapeTopology)))
													tmp_possProgenitors <- cbind(tmp_possProgenitors, t(sapply(strsplit(rownames(tmp_possProgenitors), string_lineDescent), function(x){
																												x <- as.numeric(x)
																												return( c("progenitorID"=x[1],
																														"offspringID"=x[2]) )
																										})))
													# We return now the information
													tmp_returnMat <- rbind(tmp_returnMat,
																			data.frame("progenitor_genotypeID"= tmp_possProgenitors[,"progenitorID"],
																						"offspring_genotypeID"= tmp_possProgenitors[,"offspringID"],
																						"progenitor_fitness"= tmp_possProgenitors[,"progenitor_fitness"],
																						"offspring_fitness"= tmp_possProgenitors[,"offspring_fitness"],
																						"absRank" = tmp_possProgenitors[,"absRank"],
																						"hoodSize" = tmp_possProgenitors[,"hoodSize"],
																						"hoodMin"= tmp_possProgenitors[,"hoodMin"],
																						"hoodMax"= tmp_possProgenitors[,"hoodMax"],
																						"num_altPaths"=  tmp_possProgenitors[,"num_altPaths"],
																						"relFit_altPaths" = tmp_possProgenitors[,"relFit_altPaths"],
																						"prop_maxFit" = tmp_possProgenitors[,"prop_maxFit"],
																						"progenitor_numMuts"= tmp_possProgenitors[,"progenitor_numMuts"],
																						"offspring_numMuts"= tmp_possProgenitors[,"offspring_numMuts"],
																						"transition"=rownames(tmp_possProgenitors),
																						"transitionStep"= tmpID["Step"],
																						row.names = rownames(tmp_possProgenitors),
																						stringsAsFactors = FALSE) )
												}
												# Now it's possible the same genotype transition occured multiple times which woould result in 
												# replicated rows OTHER than the transitionStep value, we search for this
												tmp_checkCols <- which(colnames(tmp_returnMat) != "transitionStep")
												tmpDuplicated <- duplicated(tmp_returnMat[,tmp_checkCols])
												if(any(tmpDuplicated)){
													# We find the paired rows for each duplicated, I setup a working reference matrix due
													# to observed issues comparing rows with my apply method.
													tmp_workMat <- matrix(gsub("[[:space:]]","",as.character(unname(unlist(tmp_returnMat[,tmp_checkCols])))),
																			ncol = ncol(tmp_returnMat[,tmp_checkCols]))
													tmpDuplicated <- lapply(which(tmpDuplicated),function(tmpRow){
																		which(apply(tmp_workMat,MARGIN=1,function(x){
																				all(x == tmp_workMat[tmpRow,tmp_checkCols])
																				}))
																	})
													# So for each unique set in the list we trim our return matrix
													tmp_removeRows <- NULL
													for(thisSet in unique(tmpDuplicated)){
														tmp_returnMat[thisSet[1],"transitionStep"] <- paste(tmp_returnMat[thisSet,"transitionStep"],collapse="_")
														tmp_removeRows <- c(tmp_removeRows,thisSet[-1])
													}
													if(!is.null(tmp_removeRows)){ tmp_returnMat <- tmp_returnMat[-tmp_removeRows,] }
												}
											} else {
												# this is the rather null case of no steps being taken
												tmp_returnMat <- data.frame("progenitor_genotypeID"= tmp_transitionMat[1,"genotypeID"],
																			"offspring_genotypeID"= NA,
																			"progenitor_fitness"=tmp_transitionMat[1,"fitness"],
																			"offspring_fitness"= NA,
																			"absRank" = NA,
																			"hoodSize" = NA,
																			"hoodMin"= NA,
																			"hoodMax"= NA,
																			"num_altPaths"=  NA,
																			"relFit_altPaths" = NA,
																			"prop_maxFit" = NA,
																			"progenitor_numMuts"= tmp_transitionMat[1,"numMuts"],
																			"offspring_numMuts"= NA,
																			"transition"=paste(tmp_transitionMat[1,"genotypeID"]),
																			"transitionStep"=0,
																			stringsAsFactors = FALSE)
												}  # This closes out building the tmp_returnMat object

											# Ok we'll now return this information
											tmpReturn <- list(list("final_estLineages"= tmp_final_estLineage,
																	"final_domLineage"= tmp_final_domLineage,
																	"transitions"= tmp_returnMat))
											names(tmpReturn)[length(tmpReturn)] <- nameEnviron(func_Index = tmp_fileString,
																								funcSplit = TRUE,
																								funcBase = envString)
											
											tmpAdd <- c(tmpAdd, tmpReturn)
											# Now we remove the environment of the file and it's data
											rm(list=c(tmp_fileString),pos=".GlobalEnv")	
										}
										tmp_replicateInfo <- c(tmp_replicateInfo, tmpAdd)
										rm(tmpAdd)
										gc()
									} # This closes out the for loop gathering information about the replicates for this job
								} # This closes out the conditional if that there is at least one file for this job
								tmpList[[thisJob]] <- tmp_replicateInfo
								# We remove the objects we created in this section
								for(thisConnection in connections_dataBase){ dbDisconnect(thisConnection) }
								# Clean up space
								rm(list=c(tmp_jobString),pos=".GlobalEnv")				
								rm(tmp_replicateInfo, connections_dataBase, func_hoodTables,func_landscapeCon,
									allow_backMutations,max_numMutations, genomeLength,db_splitTables,func_landTables)
							} # This closes out the for loop for different jobs
							
							# We can now save this object to a file
							assign(all_repSets[[thisSet]][["objName"]], tmpList, pos=".GlobalEnv")
							save(list = all_repSets[[thisSet]][["objName"]],file= paste(imageFilename["saveDir"],all_repSets[[thisSet]][["saveFile"]],sep=""))
							rm(list= all_repSets[[thisSet]][["objName"]])
							rm(tmpList)
							# Clean up space
							gc()
							# We return a confirmation that the job completed 
							return(  nameObject(func_inString = all_repSets[[thisSet]][["objName"]],
												func_inPrefix = namingString["repeatability"],
												func_splitStr = TRUE)  )
						} # This closes out the foreach loop
	
	# We now check if all the jobs completed
	tmp_checkCompleted <- sapply(tmp_missingSets,is.element,"set"=tmpCheck)
	if(all(tmp_checkCompleted)){
		# We now save the object which stores all the repeatability data object names and files locations
		save(all_repSets,file=paste(imageFilename["saveDir"],imageFilename["repeatability"],sep=""))
	} else {
		stop(paste("There was a problem with the sets of: ",paste(names(tmp_checkCompleted)[which(!tmp_checkCompleted)],collapse=" ",sep=" "),sep=""))
	}
# this closes out gathering the informaion for reapeatability
} else {
	load(paste(imageFilename["saveDir"],imageFilename["repeatability"],sep=""))
}


###################################################################################################
###################################### BEGIN OF PLOTTING ##########################################
###################################################################################################


########################################################################
############################## MULLER PLOT #############################
########################################################################
# For the Muller plots, this will assume an object I've not yet created/save
# THIS CODE IS OLD AND DESIGNED FOR A SINGLE RUN IT IS QUITE OUT OF DATE!
if(plotTypes["Muller"]){
	# First step will be to convert the object from an array with stacks into a matrix
	tmp_plotPops <- apply(info_estLines[[1]],MARGIN=3,function(thisPop){
						# From each stack we extract the popSize through time
						return( thisPop[,"popSize"] )
					})
	# We each row of this matrix will be the size of each population through time so
	# we adjust these to simply be the proportion of total population
	tmp_plotPops <- matrix(apply(tmp_plotPops,MARGIN=1,function(thisStep){ 
								return( thisStep/sum(thisStep) )
							}),byrow=TRUE,ncol=ncol(tmp_plotPops),dimnames=dimnames(tmp_plotPops))
									
	# I will now pre-process this output so that the steps are ordered top->bottom
	# This assumes that the rows are Step names.
	tmp_order <- order(as.numeric(nameTable_step(rownames(tmp_plotPops),funcSplit = TRUE)))
	tmp_plotPops <- apply(tmp_plotPops, MARGIN=2,function(thisCol){ return( thisCol[tmp_order] ) })
	
	# This object reports on the numeric step to which a row pertains
	tmp_plotSteps <- as.numeric(nameTable_step(rownames(tmp_plotPops),funcSplit = TRUE))
	
	pdf(paste(outDir,"mullerPlot.pdf",sep=""),height=4,width=6)
	# we now establish a blank plotting region into which we can plot the fate of our lineages, this uses our parameters to inform			
	plot(x=NULL,y=NULL, xlim=c(min(tmp_plotSteps),max(tmp_plotSteps)), xpd = TRUE, cex.main = 0.75, 
		ylim=c(0, 1), xlab = "Time Step", ylab = "Proportion of Popualtion Size")
	
	#### NOTE: This needs to be fixed to be robust to scenarios when the "0" population is not included......
	#### WARNING CHANGE ME WARNING CHAMGE ME 
	# In order to have th branching effect of our populations we need to define the starting point of each of our lineages
	# We do not need to define the starting population, as it will simply be the background of our plot.
	tmp_line_yStart <- sample(seq(0.05,0.95,length.out = ncol(tmp_plotPops)),replace=FALSE)
	
	# Now we go and plot each lineage separately onto the plot
	for(thisLine in 1:ncol(tmp_plotPops)){
		# We plot the lines for this population, but only need to plot the points where the lineage exists
		# I subtract the values by 1 since we have Step_0 as part of the steps included.
		thesePoints <- which(tmp_plotPops[,thisLine] > 0)
		# I now find the number of steps between each of these points
		tmp_diffPoints <- diff(thesePoints)
		# We also compute the breakPoints as being those differences > 1
		tmp_breakPoints <- which(tmp_diffPoints > 1)
		# I can create now a vector to hold the sets I'm about to create
		tmp_pointSets <- vector("mode"="list", length=length(tmp_breakPoints))
		# If there are no break points we simply use all of thesePoints
		if(length(tmp_breakPoints) == 0){
			tmp_pointSets[[1]] <- thesePoints
		# Otherwise, starting with the first of thesePoints, and between breaks, we create chunks between the diffPoint indexes which are not break points
		} else {
			tmp_pointSets[[1]] <- thesePoints[1]:thesePoints[tmp_breakPoints[1]]
			# Now if there is more than 1 break point we perform this subsetting
			if(length(tmp_breakPoints) > 1){
				for(thisSet in 2: length(tmp_breakPoints)){
					# This works by the logic that when we have a break, we want the start the next indexing at the thesePoints after the break: hence the +1
					# And carry it through to the next break point.  In this way if break points are consecutive we get that the previous index + 1 == next index
					# and so a singleton step gets properly considered.
					tmp_pointSets[[thisSet]] <- thesePoints[tmp_breakPoints[thisSet-1]+1]:thesePoints[tmp_breakPoints[thisSet]]
				}
			}
			# And lastly we check if the last break was the last point, if not we add it as a final set
			if(max(tmp_breakPoints) != length(thesePoints)){
				tmp_pointSets[[length(tmp_pointSets)+1]] <- thesePoints[(max(tmp_breakPoints)+1):length(thesePoints)]
			}
		}
		
		# Ok now that we have sets of continuous step points where the lineage existed, we plot these polygons
		for(thisSet in 1:length(tmp_pointSets)){
			# Now we will be plotting a polygon which uses the proportion of the population, centered about the start, and scaled within the region of 0,1
			# We first define what our values should be, setting this in a matrix. This then allows us to scale values
			tmp_theseProp <- matrix(c(tmp_line_yStart[thisLine] + (tmp_plotPops[tmp_pointSets[[thisSet]],thisLine]/2),
										tmp_line_yStart[thisLine] - (tmp_plotPops[tmp_pointSets[[thisSet]],thisLine]/2)),ncol=2)
			# We'll need to adjust time steps which have values outside our range							
			for(thisStep in union(which(tmp_theseProp[,1] > 1),which(tmp_theseProp[,2] < 0))){
				# We find which value is the offender and adjust the other
				if(tmp_theseProp[thisStep,1] > 1){
					tmp_theseProp[thisStep,] <- c(1,tmp_theseProp[thisStep,2] - (tmp_theseProp[thisStep,1] - 1))
				} else {
					tmp_theseProp[thisStep,] <- c(tmp_theseProp[thisStep,1] - tmp_theseProp[thisStep,2],0)	
				}
			}
			
			# Now we have a polygon that we can draw, I subtract the value 1 because Step_0 is included!
			polygon(x = c(tmp_plotSteps[tmp_pointSets[[thisSet]]], tmp_plotSteps[rev(tmp_pointSets[[thisSet]])]),
					y = c(tmp_theseProp[,1],rev(tmp_theseProp[,2])),
					col = colTransparency(allCols[thisLine], plot_colourAlpha), 
					border= allCols[thisLine], lwd = 0.1)
		}
		
	}
	dev.off()
}



########################################################################
############################ WalkLength PLOT ###########################
########################################################################

if(plotTypes["walkLength"]){
	# I now build a summary of the number of mutations in the final established lines, but also time and size of the first step.
	if(!file.exists(paste(imageFilename["saveDir"],imageFilename["walkPlot"],sep=""))){
		walk_plotInfo <- NULL
		# This uses the repSets information as that is what holds our final_estLineages compiled objects
		# The names of sets for popDemographics and repeatability should be the same so this ought to work.
		for(thisSet in names(all_repSets)){
			# We load the object which contains the information for thisSet, popDemopgraphics and repeatability.
			tmpLoad_repObj <- c("name"=all_repSets[[thisSet]][["objName"]],
							"file"=paste(imageFilename["saveDir"],all_repSets[[thisSet]][["saveFile"]],sep=""))
			load(tmpLoad_repObj["file"])
			tmp_repData <- eval(as.name(tmpLoad_repObj["name"]))
				
			# As a sanity check we ensure that the loaded object has all the same jobs we expect
			if(!all(sapply(names(tmp_repData),is.element,"set"=all_repSets[[thisSet]][["setJobs"]]))){
				stop(paste("There was a problem finding all the jobs associated with the set: ",thisSet,sep=""))
			}
			
			# I grab some parameter information for thisJob
			tmp_jobParms <- all_parmInfo[which(all_parmInfo$jobSet==thisSet)[1],c("const_mutProb","simModel","constDist","const_distParameters","distFactor")]
			# Now I grab the number of mutations that exist within the final established lineages, but also the timing and fitness increase of the first step
			for(thisJob in all_repSets[[thisSet]][["setJobs"]]){
				#Sanity check
				if(nrow(tmp_jobParms) != 1){
					stop(paste("There was a problem in trying to find parameters for this job: ",thisJob," please review",sep=""))
				}
				tmp_changeInfo <- NULL
				tmp_mutCounts <- NULL
				for(thisRun in tmp_repData[[thisJob]]){
					# This is how many mutations exist in the final dominant lineage
					tmp_domStrings <- strsplit(as.character(thisRun$final_domLineage["binaryString"]),sepString)
					tmp_mutCounts <- c(tmp_mutCounts,unname(mean(sapply(tmp_domStrings,length))))
					# We also find what is the first transition that occured in this community
					if(!any(is.na(thisRun$transitions[1,]))){
						tmp_transitionInfo <- thisRun$transitions[1,c("progenitor_fitness","offspring_fitness","transitionStep")]
							tmp_changeInfo <- cbind(tmp_changeInfo,
													c("stepSize"=mean(tmp_transitionInfo[,"offspring_fitness"]- tmp_transitionInfo[,"progenitor_fitness"]),
													"transitionTime"=mean(tmp_transitionInfo[,"transitionStep"])))
					}
				}
				# In the event there were no transitions...
				if(!is.matrix(tmp_changeInfo)){ tmp_changeInfo <- matrix(c(0,NA),nrow=2,dimnames=list(c("stepSize","transitionTime"),NULL)) }
				walk_plotInfo <- rbind(walk_plotInfo,
									cbind(tmp_jobParms,
											data.frame("combDFE"=paste(tmp_jobParms[,c("constDist","const_distParameters")],collapse="__"),
														"meanWalk"=mean(tmp_mutCounts,na.rm=TRUE),
														"sdWalk"= sd(tmp_mutCounts,na.rm=TRUE),
														"mean_stepSize"=mean(tmp_changeInfo["stepSize",],na.rm=TRUE),
														"sd_stepSize"=sd(tmp_changeInfo["stepSize",],na.rm=TRUE),
														"mean_transitionTime"=mean(tmp_changeInfo["transitionTime",],na.rm=TRUE),
														"sd_transitionTime"=sd(tmp_changeInfo["transitionTime",],na.rm=TRUE),
														"set"=thisSet,
														"job"=thisJob,
														stringsAsFactors = FALSE)))
			}
			rm(list=c(tmpLoad_repObj["name"],"tmp_repData"))
		}
		rownames(walk_plotInfo) <- NULL
		# I update my dilution to be rounded
		walk_plotInfo$distFactor <- factor(round(walk_plotInfo$distFactor,0),
										levels=unique(round(walk_plotInfo$distFactor,0)[order(round(walk_plotInfo$distFactor,0),decreasing=FALSE)]))
		# I set the mutation rate as factor
		walk_plotInfo$const_mutProb <- factor(walk_plotInfo$const_mutProb, levels= unique(walk_plotInfo$const_mutProb)[order(unique(walk_plotInfo$const_mutProb),decreasing=FALSE)])
		# We can now save this list of matrices
		save(walk_plotInfo,
			file=paste(imageFilename["saveDir"],imageFilename["walkPlot"],sep=""))
	} else {
		load(paste(imageFilename["saveDir"],imageFilename["walkPlot"],sep=""))
	}	
	
	# Ok, now for each unique, and variable, parameter type, we define a plotting character and colour
	# The start and end arguments of rainbow are set to 0 for red and 4/6 for blue - to define the colour range.
	walk_plotCols <- colTransparency(func_cols = rainbow(length(levels(walk_plotInfo$distFactor)),start=0,end=4/6),
									func_opacity = 1,
									func_scaleSaturation = 0.8,
									func_scaleValue = 0.8) 
	names(walk_plotCols) <- levels(walk_plotInfo$distFactor)
	
	# I'll now setup my pdf file for output and then build a blank plot before adding the lines
	# Each column will be for a landscape model, each row will be a combDFE
	tmp_plotDims <- c("cols"= length(unique(walk_plotInfo$simModel)),
						"rows"= length(unique(walk_plotInfo$combDFE)))
	# Here is the range of x I want plotted, it's based on the levels of distFactor
	walk_plotExpand = 3
	walk_xLim <- seq(1,length(unique(walk_plotInfo$const_mutProb)) * walk_plotExpand,length.out = length(unique(walk_plotInfo$const_mutProb)))
	# I now define the plotSpace for my x-Axis
	walk_plotSpots <- define_plotSpace(func_subDivs = walk_plotInfo[,c("const_mutProb","distFactor")], 
										func_tmpSpace = range(walk_xLim),
										func_tmpBuffer = 0.05, 
										func_tmpBuffer_max = 0.2, 
										func_tmpBuffer_min=0.01, 
										func_orderDivs = TRUE)
	# We make three plots, one for each of the walk length, 1st step size, and transition time
	walk_plotColumns <- matrix(c("meanWalk","sdWalk","mean_stepSize","sd_stepSize","mean_transitionTime","sd_transitionTime"),ncol=3,
							dimnames=list(c("mean","sd"),c("Walk Length","First Step Size","First Transition Time")))
	for(thisPlot in colnames(walk_plotColumns)){
		pdf(paste(outDir,"plot_",gsub("[[:space:]]","_",thisPlot),".pdf",sep=""),height=4 * tmp_plotDims["rows"] + 1,width=6* tmp_plotDims["cols"])
		# The number of panels to be plotted is based on the col and row product + 1 for our legend
		tmp_numPanels <- prod(tmp_plotDims)
		# We have out matrix of all parameters of interest but so that we can plot by DFE and simModel we parse in this way at first
		layout(matrix(c(seq(1,tmp_numPanels),
							rep(tmp_numPanels+1,tmp_plotDims["cols"])), 
						nrow=tmp_plotDims["rows"] + 1, byrow=TRUE), 		
				heights=c(rep(4,tmp_plotDims["rows"]),1))
		par(mai=rep(0.75,4))
		# The yLim is based on the max and min values 
		walk_yLim <- walk_plotInfo[,walk_plotColumns[c("mean","sd"),thisPlot]]
		walk_yLim <- walk_yLim[which(apply(walk_yLim,MARGIN=1,function(x){ !any(is.na(x)) })),]
		walk_yLim <- c(max(0,min(walk_yLim[,1] - 1.96 * walk_yLim[,2])),
						decimalCeiling(max(walk_yLim[,1] + 1.96 * walk_yLim[,2]),2))

		for(thisDFE in unique(walk_plotInfo$combDFE)){
			for(thisModel in unique(walk_plotInfo$simModel)){
				plot(x = NULL, y = NULL, xlim = c(min(walk_xLim),max(walk_xLim)), ylim= walk_yLim,
					xlab="Mutation Rate",ylab= paste("Mean ",thisPlot,sep=""),
					axes = FALSE, cex.lab = 1.5, cex.axis = 1.5,
					main = paste(paste("Landscape Model: ",thisModel,sep=""),
								paste("DFE: ",thisDFE,sep=""),
								sep="\n"))
				axis(side=1,at = sapply(walk_plotSpots,function(x){ mean(unlist(x)) }), labels = unique(walk_plotInfo$const_mutProb)[order(unique(walk_plotInfo$const_mutProb),decreasing = FALSE)],
					cex.axis = 1.5)
				axis(side=2,at = sapply(seq(min(walk_yLim),max(walk_yLim),length.out=7),decimalCeiling,"func_level"=2),
					cex.axis = 1.5)
				# Now we'll be plotting separate information for each set of jobs
				for(thisSet in unique(walk_plotInfo[intersect(which(walk_plotInfo$combDFE == thisDFE),which(walk_plotInfo$simModel == thisModel)),"set"])){
					thisCombo <- which(walk_plotInfo[,"set"] == thisSet)
					thisDil <- unique(walk_plotInfo[thisCombo,"distFactor"])
					thisRate <- unique(walk_plotInfo[thisCombo,"const_mutProb"])
					# The X position is given by the divided space object
					tmpX <- walk_plotSpots[[which(levels(walk_plotInfo$const_mutProb) == thisRate)]][[which(levels(walk_plotInfo$distFactor) == thisDil)]]

					polygon(x = c(tmpX, rev(tmpX)),
							y = c(rep(mean(walk_plotInfo[thisCombo,walk_plotColumns["mean",thisPlot]],na.rm=TRUE) + 1.96 * sd(walk_plotInfo[thisCombo,walk_plotColumns["sd",thisPlot]],na.rm=TRUE),2),
									rev(rep(mean(walk_plotInfo[thisCombo,walk_plotColumns["mean",thisPlot]],na.rm=TRUE) - 1.96 * sd(walk_plotInfo[thisCombo,walk_plotColumns["sd",thisPlot]],na.rm=TRUE),2))),
							col = colTransparency(walk_plotCols[as.character(thisDil)],plot_colourAlpha), border = walk_plotCols[as.character(thisDil)])
					lines(x=tmpX,
							y = rep(mean(walk_plotInfo[thisCombo,walk_plotColumns["mean",thisPlot]],na.rm=TRUE),2),
							lwd = 2, col = walk_plotCols[as.character(thisDil)],
							cex = 3)
				}
			}
		}
		# This will now be plotted in the extra space of the 
		par(mai=c(0,0,0,0))
		plot.new()
		legend("center",legend = names(walk_plotCols),
				fill = walk_plotCols,
				title = "Dilution Factor", bty = "n", horiz = TRUE,
				cex = 2)
		dev.off()
	}
} # This closes the plotting
		
	
	
########################################################################
############################# popDemo PLOT #############################
########################################################################

# This will be used to track aspects of our populations through time, these will include:
# (i) the number of competing lineages (a metric of clonal interference)
# (ii) the mean fitness
if(plotTypes["popDemographics"]){
	# I now build a summary of the population infromation extracted for each of the sets
	if(!file.exists(paste(imageFilename["saveDir"],imageFilename["plotPop"],sep=""))){
		working_parmInfo <- all_parmInfo
		working_parmInfo$combDFE <- apply(working_parmInfo[,c("constDist","const_distParameters")],MARGIN=1,paste,"collapse"="__")
		working_parmInfo$distFactor <- round(working_parmInfo$distFactor,0)
		# This is not done as an array specifically so that jobs with different numbers of generations can be compared with this script.
		extract_popCols <- c("meanFit","sd_meanFit","numLines","sd_numLines","numEstablished","sd_numEstablished","numMutants")
		# This will define how many segments of time a run will be broken into, this number should be large enough that a run is considered as many pieces 
		consider_popPeriods <- 100
		# As a safety we change this value to be the minimum of what is input OR the number of generations in a job
		consider_popPeriods <- min(consider_popPeriods,min(working_parmInfo$numGenerations))
		# This value (+1) will be multiplied by the window size, to give some overlap, so it can be any value > 0, but ideally < 1
		consider_windowOverlap <- 0.25
		# I now need to find out how many unique combinations there are for 
		plot_parmInterest <- c("numGenerations","const_mutProb","simModel","combDFE","distFactor")
		plot_jobCombos <- unique(working_parmInfo[,plot_parmInterest])
		dimnames(plot_jobCombos) <- list(paste("parmCombo",1:nrow(plot_jobCombos),sep="_"),plot_parmInterest)
		
		pop_plotInfo <- array(NA,dim=c(consider_popPeriods, length(extract_popCols), nrow(plot_jobCombos)),
								dimnames=list(NULL, extract_popCols, rownames(plot_jobCombos)))
		
		# These plots should be the first page of the output
		for(thisCombo in 1:nrow(plot_jobCombos)){
			# These will be the names of jobs which are for these parameters
			tmpJobs <- intersect(intersect(intersect(intersect(which(working_parmInfo$combDFE == plot_jobCombos[thisCombo,"combDFE"]),
																	which(working_parmInfo$numGenerations == plot_jobCombos[thisCombo,"numGenerations"])),
														which(working_parmInfo$const_mutProb == plot_jobCombos[thisCombo,"const_mutProb"])),
											which(working_parmInfo$distFactor == plot_jobCombos[thisCombo,"distFactor"])),
								which(working_parmInfo$simModel == plot_jobCombos[thisCombo,"simModel"]))
			# As long as there is one group, then we can work!
			if(length(tmpJobs) > 0){
				tmpPlot <- NULL
				for(thisSet in unique(working_parmInfo$jobSet[tmpJobs])){		
					# We load the object which contains the information for thisSet
					tmp_loadObj <- c("name"=all_popSets[[thisSet]][["objName"]],
									"file"=paste(imageFilename["saveDir"],all_popSets[[thisSet]][["saveFile"]],sep=""))
					load(tmp_loadObj["file"])
					tmpData <- eval(as.name(tmp_loadObj["name"]))
					# As a sanity check we ensure that the loaded object has all the same jobs we expect
					if(!all(sapply(names(tmpData),is.element,"set"=all_popSets[[thisSet]][["setJobs"]]))){
						stop(paste("There was a problem finding all the jobs associated with the set: ",thisSet,sep=""))
					}
					for(thisJob in intersect(working_parmInfo$jobName[tmpJobs],names(tmpData))){
						# We add the popDemographics information to our current tmpPlot information
						tmpPlot <- abind(tmpPlot,tmpData[[thisJob]][,extract_popCols],along = 3)
					}
					rm(list=c(tmp_loadObj["name"],"tmpData"))
				}
				# Ok, now we go and create a sliding window that extracts the mean, through our sets, at
				# a number of time points throughout our information
				tmp_stepSize <- ceiling(dim(tmpPlot)[1]/consider_popPeriods * (1+consider_windowOverlap))
				tmpMids <- seq(1,dim(tmpPlot)[1],length.out=consider_popPeriods+1)
				tmpMids <- round((tmpMids[-1] + tmpMids[-length(tmpMids)])/2,0)
				# Ok, so now we will step through our data, getting the mean value, and placing that within our stored output
				tmpAdd <- matrix(NA,ncol=dim(pop_plotInfo)[2],nrow=dim(pop_plotInfo)[1],dimnames=list(NULL,dimnames(pop_plotInfo)[[2]]))
				for(thisRow in 1:length(tmpMids)){
					tmpAdd[thisRow,] <- apply(tmpPlot[max(1,tmpMids[thisRow] - tmp_stepSize):min(nrow(tmpPlot),tmpMids[thisRow] + tmp_stepSize),,],MARGIN=2,mean)
				}
				pop_plotInfo[,,rownames(plot_jobCombos)[thisCombo]] <- tmpAdd
			}
		}
		# We can now save this list of matrices
		save(pop_plotInfo,working_parmInfo,consider_popPeriods,consider_windowOverlap,plot_jobCombos,
				file=paste(imageFilename["saveDir"],imageFilename["plotPop"],sep=""))
	} else {
		load(paste(imageFilename["saveDir"],imageFilename["plotPop"],sep=""))
	}	
	# Ok, now for each unique, and variable, parameter type, we define a plotting character and colour
	# I'll divide mutation rate by colour, Dilution Factor by hue, DFE + parms by plotting character, and number of generation by character size
	# I have no thoughts about how to divide this more... sorry?!?  More specific scripting required for that.
	# The start and end arguments of rainbow are set to 0 for red and 4/6 for blue - to define the colour range.
	pop_plotCols <- rainbow(length(unique(working_parmInfo$distFactor)),start=0,end=4/6)
	names(pop_plotCols) <- as.character(unique(round(working_parmInfo$distFactor,0))[order(unique(round(working_parmInfo$distFactor,0)),decreasing=FALSE)])
	pop_plotCols <- matrix(sapply(pop_plotCols,function(thisCol){ 
								return( colTransparency(func_cols = thisCol,
														func_opacity = 1,
														func_scaleSaturation = rep(0.8,length(unique(working_parmInfo$const_mutProb))),
														func_scaleValue = rep(0.8,length(unique(working_parmInfo$const_mutProb)))) ) 
							}),ncol = length(pop_plotCols),
							dimnames=list(as.character(unique(working_parmInfo$const_mutProb)),names(pop_plotCols)))
	
	pop_cexPch = seq(1,min(c(length(unique(working_parmInfo$numGenerations)),3)),length.out=length(unique(working_parmInfo$numGenerations)))
	names(pop_cexPch) <- as.character(unique(working_parmInfo$numGenerations))
	
		
	# I'll now setup my pdf file for output and then build a blank plot before adding the lines
	# Each column will be for a landscape model, each row will be a combDFE, each page will be 
	# a different plot, being either for the mean fitness, or the number of lines and/or estLines
	tmp_plotDims <- c("cols"= length(unique(working_parmInfo$combDFE)),
						"rows"= length(unique(working_parmInfo$const_mutProb)))
	# I want to create a version with and without the SD polygons
	for(this_plotType in c("noSD","wSD")){
		# If we're not plotting the SD, then we don't need it considered in our y_lim plotting values hense the alternative evalutations
		# Now as a measure of consistency I'll go and get the min/maximum values for my meanFit + 1.96 * sd, this is done by the first loop
		# the second loops take the min and max values across the vectors returned
		pop_yLim <- list("fitness"=apply(array(pop_plotInfo[,c("meanFit","sd_meanFit"),],
												dim=c(dim(pop_plotInfo)[1],2,dim(pop_plotInfo)[3])),
											MARGIN=3,function(y){
												if(any(y == -Inf)){
														tmpReplace <- which(y == -Inf,arr.ind=TRUE)
														for(tmpRow in 1:nrow(tmpReplace)){
															y[tmpReplace[tmpRow,]] <- -1
														}
													}
												tmpReturn <- apply(y,MARGIN=1,function(x){
																c(x[1] - if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0},
																x[1] + if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0})
															})
												return( c("min"=max(0,min(tmpReturn[1,])),
															"max"=max(tmpReturn[2,])) )
											}),
						"numLineages"=apply(array(log10(pop_plotInfo[,c("numLines","sd_numLines"),]),
												dim=c(dim(pop_plotInfo)[1],2,dim(pop_plotInfo)[3])),
											MARGIN=3,function(y){
													if(any(y == -Inf)){
														tmpReplace <- which(y == -Inf,arr.ind=TRUE)
														for(tmpRow in 1:nrow(tmpReplace)){
															y[tmpReplace[tmpRow,]] <- -1
														}
													}
													tmpReturn <- apply(y,MARGIN=1,function(x){
																	c(x[1] - if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0},
																	x[1] + if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0})
																})
													return( c("min"=max(0,min(tmpReturn[1,])),
																"max"=max(tmpReturn[2,])) ) 
												}),
						"numEstablished"=apply(array(log10(pop_plotInfo[,c("numEstablished","sd_numEstablished"),]),
												dim=c(dim(pop_plotInfo)[1],2,dim(pop_plotInfo)[3])),
											MARGIN=3,function(y){
													if(any(y == -Inf)){
														tmpReplace <- which(y == -Inf,arr.ind=TRUE)
														for(tmpRow in 1:nrow(tmpReplace)){
															y[tmpReplace[tmpRow,]] <- -1
														}
													}
													tmpReturn <- apply(y,MARGIN=1,function(x){
																	c(x[1] - if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0},
																	x[1] + if(this_plotType == "wSD"){ ifelse(is.na(x[2]),0,1.96*x[2]) }else{0})
																})
													return( c("min"=max(0,min(tmpReturn[1,])),
																"max"=max(tmpReturn[2,])) )
												}),
						"numMutants"=apply(array(log10(pop_plotInfo[,c("numMutants"),]),
												dim=c(dim(pop_plotInfo)[1],1,dim(pop_plotInfo)[3])),
											MARGIN=3,function(y){
													if(any(y == -Inf)){
														y[which(y == -Inf,arr.ind=TRUE)] <- -1
													}
													return( c("min"=max(0,min(y)),
																"max"=max(y)) )
												}))
		pop_yLim <- sapply(pop_yLim,function(thisSet){ c(max(0,min(thisSet[1,])),max(thisSet[2,])) },simplify=FALSE)
		# Now I round and I'll use log10 values for the number of lineages
		pop_yLim[["fitness"]] <- c(round(pop_yLim[["fitness"]][1]*0.95,2),
									decimalCeiling(pop_yLim[["fitness"]][2],2))
		pop_yLim[["numLineages"]] <- c(round(pop_yLim[["numLineages"]][1]*0.95,2),
										decimalCeiling(pop_yLim[["numLineages"]][2],2))
		pop_yLim[["numEstablished"]] <- c(round(pop_yLim[["numEstablished"]][1]*0.95,2),
										decimalCeiling(pop_yLim[["numEstablished"]][2],2))
		pop_yLim[["numMutants"]] <- c(round(pop_yLim[["numMutants"]][1]*0.95,2),
										decimalCeiling(pop_yLim[["numMutants"]][2],2))								
		# I also build a series of y-labels
		pop_yLabs <- c("fitness"="Community Mean Relative Fitness",
						"numLineages"="Num. of Genotypes (Log[10])",
						"numEstablished"="Num. of Est. Genotypes (Log[10])",
						"numMutants"="Num. of Mutants per Generation")
		
		# Here is the range of x I want plotted, it's based on the length of the simulation and so is normalised for all runs.
		pop_xLim <- c(0,1)
		
		for(this_plotResponse in names(pop_yLim)){
			pdf(paste(outDir,"plot_popDemo_",this_plotResponse,"_",this_plotType,".pdf",sep=""),height=4 * tmp_plotDims["rows"] + 1,width=6* tmp_plotDims["cols"])
			# The layout is setup such that each page will have a number of rows equal to the the number of mutation rates,
			# and a number of columns equal to the comboDFE values.  Each page will be a separate simModel
			# NOTE: This form means plots will be output by ROW so Column Type should be the inner most loop prior to plot call
			tmp_numPanels <- prod(tmp_plotDims)
			# We have out matrix of all parameters of interest but so that we can plot by DFE and simModel we parse in this way at first
			for(thisModel in unique(plot_jobCombos$simModel)){
				layout(matrix(c(seq(1,tmp_numPanels),
									rep(tmp_numPanels+1,tmp_plotDims["cols"])), 
								nrow=tmp_plotDims["rows"] + 1, byrow=TRUE), 		
						heights=c(rep(4,tmp_plotDims["rows"]),1))
				par(mai=rep(0.75,4))
				for(thisRate in unique(as.character(plot_jobCombos$const_mutProb))){
					for(thisDFE in unique(plot_jobCombos$combDFE)){
						plot(x = NULL, y = NULL, xlim = pop_xLim, ylim= pop_yLim[[this_plotResponse]],
							xlab="Simulation Time",ylab= pop_yLabs[this_plotResponse],
							cex.lab = 1.5, cex.axis = 1.5,
							main = paste(paste("Landscape Model: ",thisModel,sep=""),
										paste("DFE: ",thisDFE,sep=""),
										paste("Mutation Rate: ",thisRate,sep=""),
										sep="\n"))
						# Now I go and find all the jobs that share a particular set of parameters, and I'll be returning their mean value
						for(thisCombo in rownames(plot_jobCombos)[intersect(intersect(which(plot_jobCombos$combDFE == thisDFE),
																			which(plot_jobCombos$simModel == thisModel)),
																			which(as.character(plot_jobCombos$const_mutProb) == thisRate))]){
							thisDil <- as.character(plot_jobCombos[thisCombo,"distFactor"])
							thisGens <- as.character(plot_jobCombos[thisCombo,"numGenerations"])
							# The Y values used will depend on the type of response to be plotted
							tmpPlot <- pop_plotInfo[,if(this_plotResponse == "fitness"){
															c("meanFit","sd_meanFit")
														} else if(this_plotResponse == "numLineages"){
															c("numLines","sd_numLines" )
														} else if(this_plotResponse == "numEstablished"){
															c("numEstablished","sd_numEstablished")
														} else if(this_plotResponse == "numMutants"){
															c("numMutants","numMutants")
														},thisCombo]
							if(this_plotResponse == "numMutants"){ tmpPlot[,2] == 0 }
							# If not dealing with the fitness values I'll be plotting in log10, so I make the change and update and Inf to be 0
							if(this_plotResponse != "fitness"){
								tmpPlot <- log10(tmpPlot)
								if(any(grepl("-Inf",as.character(tmpPlot)))){
									tmpPlot[which(grepl("-Inf",as.character(tmpPlot)))] <- -1
								}
							}
							tmpX <- seq(min(pop_xLim),max(pop_xLim),length.out=nrow(tmpPlot))

							if(this_plotType == "wSD"){
								polygon(x = c(tmpX, rev(tmpX)),
										y = c(tmpPlot[,1]+ 1.96 * tmpPlot[,2],
												rev(tmpPlot[,1] - 1.96 * tmpPlot[,2])),
										col = colTransparency(pop_plotCols[thisRate,thisDil],plot_lightAlpha), border = pop_plotCols[thisRate,thisDil])
							}
							lines(x=tmpX,
									y = tmpPlot[,1],
									lwd = pop_cexPch[thisGens] * ifelse(this_plotType == "wSD",1,2.5), col = pop_plotCols[thisRate,thisDil],
									cex = 3)
						}
					}
				}
				# This will now be plotted in the extra space of the 
				par(mai=c(0,0,0,0))
				plot.new()
				legend("center",legend = colnames(pop_plotCols),
						fill = pop_plotCols[1,],
						title = "Dilution Factor", bty = "n", horiz = TRUE,
						cex = 3)
			}
			dev.off()
		} # this closes the loop for the type of plotting (wSD / noSD)
	} # This closes the loop for the type of response plotted	
}






########################################################################
########################### Repeatability PLOT #########################
########################################################################

# For the repeatability make certain you have the all_repSets object in your space
# This is going to generate two plots, the first is a 3D plot where we can visualise how dilution
# factor may be affecting repeatability within the rankFitness - relFitness space of genotypes
# The rankFitness and relFitness space is a great comparative index for genotypes since it integrates
# the impact of fitness landscape and distribution of fitness effect and makes them comparable.
# The second plot will simply look at the top (x) most repeatable steps of a landscape and see how that
# pathway is affected by dilution factor
if(plotTypes["repeatability"]){
	# This builds a processed version of the repeatability data for the convenience of plotting and stats (which is part of another script...)
	if(!file.exists(paste(imageFilename["saveDir"],imageFilename["plotRepeat"],sep=""))){
		# If there is a need to build the object then first requirement is to load all the repeatability objects
		# The information of object name and location are stored in the "all_repSets" object which should be in memory
		
		# This is our final data object for the repeatability of all transitions, by rank
		rank_plotInfo <- NULL
		domLineage_plotInfo <- NULL
		# These are the names of columns which whould not be set to numeric
		colnames_nonNumeric <- c("transition","set","job","rep", "mutsChange", "landscapeModel","DFE","domMutations","domMut_repeatability")
		
		for(thisSet in names(all_repSets)){
			# We load the object which contains the information for thisSet
			tmp_loadObj <- c("name"=all_repSets[[thisSet]][["objName"]],
							"file"=paste(imageFilename["saveDir"],all_repSets[[thisSet]][["saveFile"]],sep=""))
			load(tmp_loadObj["file"])
			tmpData <- eval(as.name(tmp_loadObj["name"]))
			# As a sanity check we ensure that the loaded object has all the same jobs we expect
			if(!all(sapply(names(tmpData),is.element,"set"=all_repSets[[thisSet]][["setJobs"]]))){
				stop(paste("There was a problem finding all the jobs associated with the set: ",thisSet,sep=""))
			}
			for(thisJob in all_repSets[[thisSet]][["setJobs"]]){
				# I grab some parameter information for thisJob
				tmp_jobParms <- all_parmInfo[which(all_parmInfo$jobName==thisJob),c("const_mutProb","simModel","constDist","const_distParameters")]
				#Sanity check
				if(nrow(tmp_jobParms) != 1){
					stop(paste("There was a problem in trying to find parameters for this job: ",thisJob," please review",sep=""))
				}
				tmpTransitions <- NULL
				tmp_domBinary <- NULL
				for(thisInfo in names(tmpData[[thisJob]])){
					tmpTransitions <- rbind(tmpTransitions,
											data.frame(cbind(tmpData[[thisJob]][[thisInfo]][["transitions"]],
															"set"=thisSet,
															"job"= thisJob,
															"rep"=thisInfo)) ) 
					# Now we deal with the dominant lineage mutation information, by using the binary string for that lineage
					# We'll also return the absRank and fitness information
					tmp_trans_refRow <- which(tmpTransitions[,"offspring_genotypeID"] == as.numeric(tmpData[[thisJob]][[thisInfo]][["final_domLineage"]]["genotypeID"]))[1] 
					# If the dominant lineage is WT, then there won't be a tmp_trans_refRow value (it will be NA), 
					# if neither of these are the case then we return what information exists, else there have been no mutations
					if(!is.na(tmp_trans_refRow) && !all(unname(tmpData[[thisJob]][[thisInfo]][["final_domLineage"]]) == c("0",""))){
						tmp_domBinary <- rbind(tmp_domBinary,
												data.frame(matrix(tmpData[[thisJob]][[thisInfo]][["final_domLineage"]], nrow = 1,
																dimnames=list(tmpData[[thisJob]][[thisInfo]][["final_domLineage"]]["genotypeID"],
																				names(tmpData[[thisJob]][[thisInfo]][["final_domLineage"]]))),
															"fitness"=tmpTransitions[tmp_trans_refRow,"offspring_fitness"],
															"absRank"=tmpTransitions[tmp_trans_refRow,"absRank"],
															"set"=thisSet,
															"job"= thisJob,
															"rep"=thisInfo) )
					}
				}
				tmpTransitions <- data.frame(apply(tmpTransitions,MARGIN=2,unlist))
				# I go and set all columns except for transition and set,job,rep back to numeric
				for(thisCol in setdiff(colnames(tmpTransitions), colnames_nonNumeric)){ tmpTransitions[,thisCol] <- as.numeric(tmpTransitions[,thisCol]) }
				# We remove any rows with NA values and check there is anything left
				tmpTransitions <- tmpTransitions[which(apply(tmpTransitions,MARGIN=1,function(x){ !any(is.na(x)) })),]
				if(nrow(tmpTransitions) != 0){
					# I now build a data.frame of all the individual transitions
					rank_plotInfo <- rbind(rank_plotInfo,
											 t(sapply(unique(tmpTransitions[,"transition"]),function(thisTransition){ 
													# We find the row(s) which relate to this transition
													tmpRows <- which(tmpTransitions[,"transition"] == thisTransition)
													# We want to return the information concerning the parent and child fitness, numMuts, rank and dilution factor,
													# repeatability calculates the number of unique reps in tmprows vs. the number of replicates for thisJob
													return( data.frame("transition"= thisTransition,
																		"parentFitness"= tmpTransitions[tmpRows[1],"progenitor_fitness"],
																		"childFitness"= tmpTransitions[tmpRows[1],"offspring_fitness"],
																		"parent_numMuts"= tmpTransitions[tmpRows[1],"progenitor_numMuts"],
																		"child_numMuts"= tmpTransitions[tmpRows[1],"offspring_numMuts"],
																		"stepType"= tmpTransitions[tmpRows[1],"offspring_numMuts"] - tmpTransitions[tmpRows[1],"progenitor_numMuts"],
																		"mutsChange"= paste(tmpTransitions[tmpRows[1],c("progenitor_numMuts","offspring_numMuts")],
																							collapse= string_lineDescent),
																		"absRank"= tmpTransitions[tmpRows[1],"absRank"],
																		"Dilution" = all_parmInfo[which(all_parmInfo$jobName == thisJob),"distFactor"],
																		"repeatability"= length(unique(tmpTransitions[tmpRows,"rep"]))/length(tmpData[[thisJob]]),
																		"set"= thisSet,
																		"job"=thisJob,
																		"mutationRate"=unname(tmp_jobParms["const_mutProb"]),
																		"landscapeModel"=unname(tmp_jobParms["simModel"]),
																		"DFE"=paste(tmp_jobParms[1,c("constDist","const_distParameters")],collapse="__"),
																		stringsAsFactors = FALSE) )
								
												})) )
				}
				# We now also return repeatability information for the domLineage_plotInfo
				if(!is.null(tmp_domBinary)){
					if(nrow(tmp_domBinary) != 0){
						# I now build a data.frame of the information for the number of different end dominant genotypes,
						# the number of unique mutations found across all these genotypes, and the repeatability of each mutation.
						tmp_domMuts <- strsplit(tmp_domBinary[,"binaryString"],sepString)
						tmp_unique_domMuts <- unique(unlist(tmp_domMuts))
						tmp_unique_domMuts <- tmp_unique_domMuts[order(as.numeric(tmp_unique_domMuts))]
						# Now I find the repeatability of each of these mutations by looking at how many end genotypes posessed it
						tmp_unique_domRepeatability <- sapply(tmp_unique_domMuts,function(thisMut){
															length(which(sapply(tmp_domMuts,function(thisDom){ is.element(thisMut,thisDom) })))/length(tmpData[[thisJob]])
														})
						# I build something similar but for genotypeID
						tmp_unique_domID <- table(tmp_domBinary[,"genotypeID"])/length(tmpData[[thisJob]])
						tmp_unique_domID <- tmp_unique_domID[order(tmp_unique_domID,decreasing=TRUE)]
						
						tmp_jobFrame <- data.frame("numGenotypes"=length(unique(tmp_domBinary[,"genotypeID"])),
													"numMutations"=length(tmp_unique_domMuts),
													"domMutations"=paste(names(tmp_unique_domRepeatability),sep="",collapse=sepString),
													"domMut_repeatability"=paste(tmp_unique_domRepeatability,sep="",collapse=sepString),
													"domIDs"=paste(names(tmp_unique_domID),sep="",collapse=sepString),
													"domID_repeatability"=paste(tmp_unique_domID,sep="",collapse=sepString),
													"meanFitness"=mean(tmp_domBinary[,"fitness"], na.rm=TRUE),
													"sdFitness"=sd(tmp_domBinary[,"fitness"], na.rm=TRUE),
													"mean_absRank"=mean(tmp_domBinary[,"absRank"], na.rm=TRUE),
													"sd_absRank"=sd(tmp_domBinary[,"absRank"], na.rm=TRUE),
													"Dilution" = all_parmInfo[which(all_parmInfo$jobName == thisJob),"distFactor"],
													"mutationRate"=unname(tmp_jobParms["const_mutProb"]),
													"landscapeModel"=unname(tmp_jobParms["simModel"]),
													"DFE"=paste(tmp_jobParms[1,c("constDist","const_distParameters")],collapse="__"),
													"set"= thisSet,
													"job"=thisJob,						
													stringsAsFactors = FALSE)
						# And we past it onto the existing information from other runs
						domLineage_plotInfo <- rbind(domLineage_plotInfo,tmp_jobFrame)
					}
				}
			}
			# We can now remove the loaded object from memory
			rm(list=tmp_loadObj["name"])
		}
		# We now fix up the output... again.
		rank_plotInfo <- suppressWarnings(data.frame(apply(rank_plotInfo,MARGIN=2,unlist),fix.empty.names=FALSE,check.names=FALSE))
		rownames(rank_plotInfo) <- rownames(rank_plotInfo) <- NULL
		for(thisCol in setdiff(colnames(rank_plotInfo), colnames_nonNumeric)){ rank_plotInfo[,thisCol] <- as.numeric(rank_plotInfo[,thisCol]) }
		
		domLineage_plotInfo <- suppressWarnings(data.frame(apply(domLineage_plotInfo,MARGIN=2,unlist),fix.empty.names=FALSE,check.names=FALSE))
		rownames(domLineage_plotInfo) <- rownames(domLineage_plotInfo) <- NULL
		for(thisCol in setdiff(colnames(domLineage_plotInfo), colnames_nonNumeric)){ domLineage_plotInfo[,thisCol] <- as.numeric(domLineage_plotInfo[,thisCol]) }
		
		# I'll now also build some summary information concerning the values between runs.  
		# We divide by set and abs rank, calculating mean and SD of repeatability
		rank_summaryStats <- sapply(unique(rank_plotInfo[,"set"]),function(thisSet){
										# We define the rows of this set 
										tmpRows <- which(rank_plotInfo[,"set"] == thisSet)
										# I'll subdivide by the type of change in mutations
										sapply(unique(rank_plotInfo[tmpRows,"mutsChange"]),function(thisChange){
											# These are the rows we're working with
											tmp_subRows <- tmpRows[which(rank_plotInfo[tmpRows,"mutsChange"] == thisChange)]
											# This is for convenience
											tmpRanks <- unique(rank_plotInfo[tmp_subRows,"absRank"])[order(unique(rank_plotInfo[tmp_subRows,"absRank"]))]
											# Now I build a matrix that looks at all the ordered rank values and gets the 
											# mean and sd of repeatability measured
											tmpReturn <- matrix(,nrow=length(tmpRanks),ncol=3,dimnames=list(as.character(tmpRanks),c("numJobs","Mean","SD")))
											for(thisRank in tmpRanks){
												tmp_useRows <- tmp_subRows[which(rank_plotInfo[tmp_subRows,"absRank"]==thisRank)]
												tmpReturn[as.character(thisRank),] <-  c("numJobs"=length(unique(rank_plotInfo[tmp_useRows,"job"])),
																							"Mean"= mean(rank_plotInfo[tmp_useRows,"repeatability"]),
																							"SD"= sd(rank_plotInfo[tmp_useRows,"repeatability"]))
											}
											# If there was only a single instance then SD will be NA, we switch these to zero values
											if(any(is.na(tmpReturn[,"SD"]))){ tmpReturn[which(is.na(tmpReturn[,"SD"])),"SD"] <- 0 }
											return( tmpReturn )
										},simplify="list")
									},simplify="list")
		
		save(rank_plotInfo, domLineage_plotInfo, rank_summaryStats,file=paste(imageFilename["saveDir"],imageFilename["plotRepeat"],sep=""))
	} else {
		# This should hold an object called rank_plotInfo
		load(paste(imageFilename["saveDir"],imageFilename["plotRepeat"],sep=""))	
	}
	
	##########################################################################################
	################################## Heat Map Plots ########################################
	##########################################################################################
	
	########################### By mutational step #########################
	# We use the formula method for plotting but we'll transform out data a bit for clarity due to the 
	# dispersion of data within the considered ranges
	tmpData <- rank_plotInfo
	tmp_rankLimit <- 50
	tmp_rankDivide <- 5
	tmp_maxRepeatability <- 1.0
	tmp_repeatabilityGroups <- 51
	tmpData <- tmpData[which(tmpData[,"absRank"] <= tmp_rankLimit),]
	tmpData <- tmpData[which(tmpData[,"repeatability"] <= tmp_maxRepeatability),]
	tmpData <- tmpData[order(tmpData[,"Dilution"]),]
	tmpData[,"Dilution"] <- factor(round(tmpData[,"Dilution"],0),levels=unique(round(tmpData[,"Dilution"],0))[order(unique(round(tmpData[,"Dilution"],0),decreasing=FALSE))])
	tmpData[,"mutationRate"] <- factor(tmpData[,"mutationRate"],levels=unique(tmpData[,"mutationRate"])[order(unique(tmpData[,"mutationRate"]),decreasing=FALSE)])
	for(thisCol in c("DFE","landscapeModel","mutsChange")){ tmpData[,thisCol] <- factor(tmpData[,thisCol]) }
	
	# I will sub-divide these plots on a single page by the mutation probability, and DFE, but then also
	# by different pages given the change type.  These next three values should be the columns of rank_plotInfo used to separate
	plot_pageTypes <- "mutsChange"
	plot_colTypes <- "mutationRate"
	plot_rowTypes <- "DFE"
	
	plot_numCols <- length(unique(tmpData[,plot_colTypes]))
	plot_numRows <- length(unique(tmpData[,plot_rowTypes]))
	
	# Now a heatmap plot so that we can better explore the interaction between our factors
	pdf(paste(outDir,"plot_repeatability_heatmap_byChange.pdf",sep=""),width=6*plot_numCols,height=4*plot_numRows)
	
	for(thisPage in unique(tmpData[,plot_pageTypes])){
		tmp_pageRows <- which(tmpData[,plot_pageTypes] == thisPage)
		if(length(tmp_pageRows) > 0){
		for(thisModel in unique(tmpData[tmp_pageRows,"landscapeModel"])){
			listPlots <- list()
			tmp_modelRows <- tmp_pageRows[which(tmpData$landscapeModel[tmp_pageRows] == thisModel)]
			if(length(tmp_modelRows) > 0){
			for(thisRow in levels(tmpData[,plot_rowTypes])){
				for(thisCol in levels(tmpData[,plot_colTypes])){
					tmp_colRows <- tmp_modelRows[intersect(which(as.character(tmpData[tmp_modelRows,plot_rowTypes]) == thisRow),
															which(as.character(tmpData[tmp_modelRows,plot_colTypes]) == thisCol))]
					# Make certain there are at least some rows to work with....
					if(length(tmp_colRows) > 0){
						listPlots <- c(listPlots,
										list(levelplot(repeatability ~ absRank + Dilution, data = tmpData[tmp_colRows,], 
													col.regions=matlab.like(tmp_repeatabilityGroups), xlab = "Rank Fitness", ylab = "Dilution Factor",
													xlim=c(0,tmp_rankLimit), cut = tmp_repeatabilityGroups, at = seq(0,tmp_maxRepeatability,length.out=tmp_repeatabilityGroups),
													scales = list("x"=list(at = c(1,seq(tmp_rankDivide,tmp_rankLimit,by=tmp_rankDivide))),
																	"y"=list(at = as.numeric(unique(tmpData[,"Dilution"])), labels= levels(tmpData[,"Dilution"]))),
													main = paste(paste("Change: ",thisPage,sep=""),
																paste("Model: ",thisModel,sep=""),
																paste("DFE: ",thisRow,sep=""),
																paste("mutRate: ",thisCol,sep=""),
																sep="\n"),
													cex.main = 0.4)))
					}
				}
			}}
			# We now plot for thisPage
			for(thisPlot in 1:length(listPlots)){
				# The split works on a vector of (colPos,rowPos,numCols,numRows)
				print(listPlots[[thisPlot]],
						split=c(thisPlot - floor((thisPlot-1)/plot_numCols)*plot_numCols,ceiling(thisPlot/plot_numCols),plot_numCols, plot_numRows),
						more= thisPlot < length(listPlots))
			}
		}}			
	}
	dev.off()
	
	
	########################### By mutations in dominant genotype #########################
	# We use the formula method for plotting but we'll transform out data a bit for clarity due to the 
	# dispersion of data within the considered ranges
	
	# I now go and parse out the repeatability information for these dominant types
	tmpData <- NULL
	domPlot_keepCols <- c("numGenotypes","numMutations","Dilution","mutationRate","landscapeModel","DFE","set","job")
	for(thisRow in 1:nrow(domLineage_plotInfo)){
		tmpReps <- as.numeric(unlist(strsplit(domLineage_plotInfo[thisRow,"domMut_repeatability"],sepString)))
		tmpAdd <- suppressWarnings(cbind(data.frame(domLineage_plotInfo[thisRow,domPlot_keepCols]),
													"repeatability"=tmpReps))
		tmpAdd <- cbind(tmpAdd[order(tmpAdd[,"repeatability"],decreasing=TRUE),],
						"absRank" = 1:length(tmpReps))
		tmpData <- rbind(tmpData, tmpAdd)
	}
	
	tmp_rankLimit <- 30
	tmp_rankDivide <- 5
	tmp_maxRepeatability <- 1.0
	tmp_repeatabilityGroups <- 51
	tmpData <- tmpData[which(tmpData[,"absRank"] <= tmp_rankLimit),]
	tmpData <- tmpData[which(tmpData[,"repeatability"] <= tmp_maxRepeatability),]
	tmpData <- tmpData[order(tmpData[,"Dilution"]),]
	tmpData[,"Dilution"] <- factor(round(tmpData[,"Dilution"],0),levels=unique(round(tmpData[,"Dilution"],0))[order(unique(round(tmpData[,"Dilution"],0),decreasing=FALSE))])
	tmpData[,"mutationRate"] <- factor(tmpData[,"mutationRate"],levels=unique(tmpData[,"mutationRate"])[order(unique(tmpData[,"mutationRate"]),decreasing=FALSE)])
	for(thisCol in c("DFE","landscapeModel")){ tmpData[,thisCol] <- factor(tmpData[,thisCol]) }
	
	# I will sub-divide these plots on a single page by the mutation probability, and DFE, but then also
	# by different pages given the landscapeModel which ought to affect final genotype.  These next three values should be the columns of rank_plotInfo used to separate
	plot_pageTypes <- "landscapeModel"
	plot_colTypes <- "mutationRate"
	plot_rowTypes <- "DFE"
	
	plot_numCols <- length(unique(tmpData[,plot_colTypes]))
	plot_numRows <- length(unique(tmpData[,plot_rowTypes]))
	
	# Now a heatmap plot so that we can better explore the interaction between our factors
	pdf(paste(outDir,"plot_repeatability_heatmap_domMuts.pdf",sep=""),width=6*plot_numCols,height=4*plot_numRows)
	
	for(thisPage in unique(tmpData[,plot_pageTypes])){
		tmp_pageRows <- which(tmpData[,plot_pageTypes] == thisPage)
		if(length(tmp_pageRows) > 0){
			listPlots <- list()
			for(thisRow in levels(tmpData[,plot_rowTypes])){
				for(thisCol in levels(tmpData[,plot_colTypes])){
					tmp_colRows <- tmp_pageRows[intersect(which(as.character(tmpData[tmp_pageRows,plot_rowTypes]) == thisRow),
															which(as.character(tmpData[tmp_pageRows,plot_colTypes]) == thisCol))]
					# Make certain there are at least some rows to work with....
					if(length(tmp_colRows) > 0){
						listPlots <- c(listPlots,
										list(levelplot(repeatability ~ absRank + Dilution, data = tmpData[tmp_colRows,], 
													col.regions=matlab.like(tmp_repeatabilityGroups), xlab = "Rank Order Repeatability", ylab = "Dilution Factor",
													xlim=c(0,tmp_rankLimit), cut = tmp_repeatabilityGroups, at = seq(0,tmp_maxRepeatability,length.out=tmp_repeatabilityGroups),
													scales = list("x"=list(at = c(1,seq(tmp_rankDivide,tmp_rankLimit,by=tmp_rankDivide))),
																	"y"=list(at = as.numeric(unique(tmpData[,"Dilution"])), labels= levels(tmpData[,"Dilution"]))),
													main = paste("Landscape: ",thisPage,"\n mutRate: ",thisCol,"\n DFE: ",thisRow,sep=""),
													cex.main = 0.5)))
					}
				}
			}
			# We now plot for thisPage
			for(thisPlot in 1:length(listPlots)){
				# The split works on a vector of (colPos,rowPos,numCols,numRows)
				print(listPlots[[thisPlot]],
						split=c(thisPlot - floor((thisPlot-1)/plot_numCols)*plot_numCols,ceiling(thisPlot/plot_numCols),plot_numCols, plot_numRows),
						more= thisPlot < length(listPlots))
			}
		}			
	}
	dev.off()
	
}



########################################################################
######################## landscapeCompare PLOT #########################
########################################################################

# This plot is to compare the distribution of fitness effects depending on the simModel and DFE (function + parameters)
if(plotTypes["landscapeCompare"]){
	if(!file.exists(paste(imageFilename["saveDir"],imageFilename["plotLandscape"],sep=""))){
		# We find all the directories which use the names of our divided files
		tmpDirs <- list.dirs(path=workDir,full.names=FALSE,recursive=FALSE)
		tmpDirs <- unique(unlist(lapply(names(all_dividedFiles), function(thisName){
							return( tmpDirs[which(grepl(thisName,tmpDirs))] )
					})))
			
		# I need to create an outer list object which stores the information related to each number of mutations
		# of interest.  In this case I'll only concern myself with say the first four levels, but could change this....
		### NOTE: For well explored landscapes, higher number of mutations may result in vectors too long for R.
		all_mutLevels <- seq(1,5,by=1)
		# I also define the max and min expected fitness values, and then the number of breaks
		# We also define what percent of the highest valued fitnesses we want used for the highFit sets
		landscape_histParms <- list("normFit"=c("min"=0,"max"=2,"breaks"=100,"perc"=1),
									"highFit"=c("min"=0.9,"max"=2,"breaks"=100,"perc"=0.1),
									"veryHFit"=c("min"=1,"max"=2,"breaks"=100,"perc"=0.01))
								
		
		# I now define the number of unique simModel and DFE combinations that exist
		landscapeTypes <- expand.grid("simModel"=unique(all_parmInfo$simModel),
								"combDFE"=unique(apply(all_parmInfo[,c("constDist","const_distParameters")],MARGIN=1,paste,"collapse"="__")),
								stringsAsFactors=FALSE)
		attr(landscapeTypes,"out.attrs") <- NULL
		all_landScapeInfo <- sapply(paste("numMuts",all_mutLevels,sep="_"),function(x){
										sapply(unique(landscapeTypes$simModel),function(y){
											sapply(unique(landscapeTypes$combDFE),function(z){
												NA
											},simplify=FALSE)
										},simplify=FALSE)
									},simplify=FALSE)
		# I want my output grouped by number of mutations so we use this as the outer organising element
		for(this_mutLevel in all_mutLevels){
			# Because the directory names reflect the jobName, I can use the all_parmInfo$jobName info to group my jobs
			for(thisLandscape in 1:nrow(landscapeTypes)){
				tmp_landDirs <- tmpDirs[which(is.element(tmpDirs,all_parmInfo$jobName[intersect(intersect(which(all_parmInfo$constDist == strsplit(landscapeTypes[thisLandscape,"combDFE"],"__")[[1]][1]),
																											which(all_parmInfo$const_distParameters == strsplit(landscapeTypes[thisLandscape,"combDFE"],"__")[[1]][2])),
																								which(all_parmInfo$simModel == landscapeTypes[thisLandscape,"simModel"]))] ))]
				# We then create a list of all the fitness effects, as measured from the landscape file, in each directory
				# For memory reasons we save the information as tables which are used to manually make histogram plots
				tmpFitnesses <- NULL
				# I use a counter so that every ten sets that have been combined we re-group
				tmp_groupCounter <- 0
				for(thisDir in tmp_landDirs){
					# For this directory we find the landscape file
					tmpFile <- list.files(path=paste(workDir,thisDir,sep=""),pattern="Landscape_1.sqlite")
					# If there is not exactly one file then we've had a problem...
					if(length(tmpFile) != 1){
						stop(paste("Problem with respect to finding the file for ",thisDir,sep=""))
					}
					# Otherwise we return the fitness value for the appropriate number of mutations found in this file
					tmpCon <- resetDatabase(func_conName = paste(workDir,thisDir,"/",tmpFile,sep=""))
					tmpTable <- dbListTables(tmpCon)
					tmp_splitTables <- all_parmInfo[which(all_parmInfo[,"jobName"] == thisDir),"db_splitTables"]
					# This search differs based on whether or not we were splitting tables (object db_splitTables)
					tmpTable <- tmpTable[which(grepl(nameTable(this_mutLevel,
																func_baseString = "numMutations",
																func_sepString = "_",
																func_subNaming= tmp_splitTables),tmpTable))]
					tmpReturn <- if(length(tmpTable) > 0){
									unname(unlist(dbGetQuery(tmpCon,paste(paste("SELECT fitness FROM ",
																			tmpTable,
																			collapse=" UNION ",
																			sep=""),
																		" ORDER BY fitness ",
																		sep=""))))
								} else {
									c(0,0)
								}
					dbDisconnect(tmpCon)
					tmpFitnesses <- c(tmpFitnesses, table(tmpReturn) )
					tmp_groupCounter <- tmp_groupCounter + 1
					if(tmp_groupCounter >= 25){
						tmpFitnesses <- sapply(unique(names(tmpFitnesses)),function(thisName){ 
											return( sum(tmpFitnesses[which(names(tmpFitnesses) == thisName)]) )
										})
						tmp_groupCounter <- 0
					}
				}
				# We now trim our vector to those values within our range
				tmpAdd <- sapply(landscape_histParms,function(theseParms){
									# This subsets our vector into elements of interest and uses cuts to extract those values within
									# the percent highest values we care for.
									theseFit <- tmpFitnesses[intersect(which(as.numeric(names(tmpFitnesses)) >= theseParms["min"]),
																		which(as.numeric(names(tmpFitnesses)) <= theseParms["max"]))]
									theseFit <- theseFit[ceiling(length(theseFit) * (1-theseParms["perc"])):length(theseFit)]
									# Now I need to cut the named values into "breaks" number of ranges
									tmpSeq <- seq(theseParms["min"],theseParms["max"],length.out = theseParms["breaks"]+1)
									tmpRanges <- cut(as.numeric(names(theseFit)),breaks = tmpSeq)
									# I now build a matrix which will represent the lower and upper bound for my rectangle
									# as well as the count.  This is based on values collected.
									theseFit <- sapply(levels(tmpRanges),function(thisLevel){ sum(theseFit[which(tmpRanges == thisLevel)]) })
									return( list("hist"= theseFit,
												"splits"=tmpSeq) )
								},simplify=FALSE)
				all_landScapeInfo[[paste("numMuts",this_mutLevel,sep="_")]][[landscapeTypes[thisLandscape,"simModel"]]][[landscapeTypes[thisLandscape,"combDFE"]]] <- tmpAdd
			}
		}
		save(all_landScapeInfo,landscape_histParms,landscapeTypes,
			file=paste(imageFilename["saveDir"],imageFilename["plotLandscape"],sep=""))
	} else {
		load(paste(imageFilename["saveDir"],imageFilename["plotLandscape"],sep=""))
	}
	
	# Now all we're looking to do is plot the different distributions, I'll do the density for clarity
	# The colouring I'll do will be based on the combDFE and simModel information
	land_plotCols <- colTransparency(rainbow(length(all_landScapeInfo),start=0,end=4/6),plot_lightAlpha)
	names(land_plotCols) <- names(all_landScapeInfo)

	# I'll be doing a new plot for each of the levels of mutations scored
	plotDims <- c(length(unique(landscapeTypes$simModel)),length(unique(landscapeTypes$combDFE)))
	######## HISTORGRAM PLOTS #######
	for(thisSet in names(landscape_histParms)){
		pdf(paste(outDir,"plot_landscape_",thisSet,"_hist.pdf",sep=""),width=6*plotDims[2],height=4*plotDims[1]+1)
		layout(matrix(c(seq(1, prod(plotDims)),
							rep(prod(plotDims)+1,plotDims[2])), 
					nrow=plotDims[1] + 1, byrow=TRUE),
				heights=c(rep(4,plotDims[1]),1))
		par(mai=rep(0.75,4))
		
		for(thisModel in names(all_landScapeInfo[[1]])){
			for(thisDFE in names(all_landScapeInfo[[1]][[thisModel]])){
				# I build back together the ploting information
				tmp_plotData <- array(,dim=c(landscape_histParms[[thisSet]]["breaks"],3,length(all_landScapeInfo)),
										dimnames=list(NULL,c("x0","x1","y"),names(all_landScapeInfo)))
				tmpCounts <- vector(mode="numeric",length=length(all_landScapeInfo))
				names(tmpCounts) <- names(all_landScapeInfo)
				# We recollect the plot data which includes the range of values and counts of grouped values	
				for(this_numMuts in names(all_landScapeInfo)){
					tmpSplits <- all_landScapeInfo[[this_numMuts]][[thisModel]][[thisDFE]][[thisSet]]$splits
					tmpCounts[this_numMuts] <- sum(all_landScapeInfo[[this_numMuts]][[thisModel]][[thisDFE]][[thisSet]]$hist)
					tmp_plotData[,,this_numMuts] <- matrix(c(tmpSplits[-length(tmpSplits)],
															tmpSplits[-1],
															if(tmpCounts[this_numMuts] != 0){
																all_landScapeInfo[[this_numMuts]][[thisModel]][[thisDFE]][[thisSet]]$hist/tmpCounts[this_numMuts]
															} else {
																all_landScapeInfo[[this_numMuts]][[thisModel]][[thisDFE]][[thisSet]]$hist
															}),
															ncol=3)
				}
				plot(x= NULL, y=NULL, xlim = c(min(tmp_plotData[,"x0",]),max(tmp_plotData[,"x1",])),
						ylim = c(0,max(tmp_plotData[,"y",])), xlab = "Fitness", ylab = "Proportion",
						main = paste(paste("Model: ",thisModel,sep=""),
										paste("DFE: ",thisDFE,sep=""),
										paste("Prop Fit: ",landscape_histParms[[thisSet]]["perc"],sep=""),
										sep="\n"))
				# We recollect the plot data which includes the range of values and counts of grouped values	
				for(this_numMuts in names(all_landScapeInfo)){
					for(thisBar in 1:nrow(tmp_plotData)){
						polygon(x= c(rep(tmp_plotData[thisBar,"x0",this_numMuts],2),rep(tmp_plotData[thisBar,"x1",this_numMuts],2)),
								y=c(0,rep(tmp_plotData[thisBar,"y",this_numMuts],2),0),
								col=land_plotCols[this_numMuts], border = "black")
					}
				}
				legend("topright",legend=tmpCounts,fill=land_plotCols,title="Counts",bty="n")
			}
		}
		par(mai=c(0,0,0,0))
			plot.new()
			legend("center",legend = sapply(strsplit(names(all_landScapeInfo),"_"),function(x){ x[2] }),
					col = land_plotCols,
					title = "Number of Mutations", bty = "n", horiz = TRUE,
					cex = 3, lwd = 4)
		dev.off()
	}
	
}



########################################################################
############################ Entropy PLOT ##############################
########################################################################


# This will use Szendro 2013's calculation of entropy to measure the determinism of evolution, divided by each parameter combinations
if(plotTypes["entropy"]){
	# This builds a processed version of the entropy data for the convenience of plotting and stats (which is part of another script...)
	if(!file.exists(paste(imageFilename["saveDir"],imageFilename["plotEntropy"],sep=""))){
		# This is our final data object for the repeatability of all transitions, by rank
		entropy_plotInfo <- NULL
		for(thisSet in names(all_repSets)){
			# We load the object which contains the information for thisSet
			tmp_loadObj <- c("name"=all_repSets[[thisSet]][["objName"]],
							"file"=paste(imageFilename["saveDir"],all_repSets[[thisSet]][["saveFile"]],sep=""))
			load(tmp_loadObj["file"])
			tmpData <- eval(as.name(tmp_loadObj["name"]))
			# As a sanity check we ensure that the loaded object has all the same jobs we expect
			if(!all(sapply(names(tmpData),is.element,"set"=all_repSets[[thisSet]][["setJobs"]]))){
				stop(paste("There was a problem finding all the jobs associated with the set: ",thisSet,sep=""))
			}
			for(thisJob in all_repSets[[thisSet]][["setJobs"]]){
				# I grab some parameter information for thisJob
				tmp_jobParms <- all_parmInfo[which(all_parmInfo$jobName==thisJob),c("const_mutProb","simModel","constDist","const_distParameters")]
				#Sanity check
				if(nrow(tmp_jobParms) != 1){
					stop(paste("There was a problem in trying to find parameters for this job: ",thisJob," please review",sep=""))
				}
				# now for each replicate associated to thisJob, we need to go get the line of descent for the final_domLineage
				tmpAdd <- foreach(thisFile = names(tmpData[[thisJob]]),.combine="rbind")%dopar%{
					# We identify what was the final established lineage for this replicate run
					tmp_domLineage <- tmpData[[thisJob]][[thisFile]]$final_domLineage[c("genotypeID","binaryString")]
					
					tmpFile <- list.files(path=paste(workDir,thisJob,"/",sep=""),
											pattern=paste(thisFile,saveSuffixes["run"],sep=""),
											full.names=TRUE)
					# We make certain that we've identified a single file, and if so we load it to grab the line of descent information
					if(length(tmpFile) > 1){
						stop(paste("Problem finding a unique file to load data from for: ",thisFile,sep=""))
					#} else if (length(tmpFile) == 0){
					#	tmpAdd <- tmpAdd[-which(rownames(tmpAdd) == thisFile),]
					} else if (length(tmpFile) == 1){
						tmp_fileString <- nameEnviron(tmpFile,funcBase = envString)
						assign(tmp_fileString, new.env(), pos = ".GlobalEnv")
						# This should load objects called info_estLines and runDemographics
						load(tmpFile, envir = eval(as.name(tmp_fileString)) )
						# What we need is the info_estLines$end_Lines_of_Descent for the tmpDom lineage, checking that the information is consistent...
						if(!is.element(tmp_domLineage["genotypeID"],names(eval(as.name(tmp_fileString))$info_estLines$end_Lines_of_Descent))){
							stop(paste("Problem with consistency of the line of descent and listed dominant lineage for: ",thisFile,sep=""))
						}
						tmp_line_of_descent <- as.character(eval(as.name(tmp_fileString))$info_estLines$end_Lines_of_Descent[[tmp_domLineage["genotypeID"]]])
						# We now try and extract the first transition in this line of descent
						tmp_dom_transitions <- strsplit(tmp_line_of_descent,string_lineDescent)
						# Now if there was no transition we should be left with the WT genotype ID of "0" as the only element... otherwise we have a first step.
						tmp_firstTransition <- unlist(lapply(tmp_dom_transitions,function(x){ 
													if(length(x) == 1){
														paste(rep(x,2),collapse=string_lineDescent)
													} else {
														paste(x[1:2],collapse=string_lineDescent) 
													} }))
						#tmpAdd[thisFile,] <- unname(c(tmp_domLineage,tmp_line_of_descent,tmp_firstTransition))
						rm(list=c(tmp_fileString),pos=".GlobalEnv")
						gc()
						tmp_mergeInfo <- t(sapply(1:length(tmp_line_of_descent),function(x) { 
												unname(c(tmp_domLineage,tmp_line_of_descent[x],tmp_firstTransition[x], thisFile))
											},USE.NAMES=FALSE))
						#tmpAdd <- rbind(tmpAdd,tmp_mergeInfo)
						return( tmp_mergeInfo )
					} else {
						return( NULL )
					}
				}
				# We build up the information concerning the entropy values, this first requires that we work our probabilities
				tmp_probInfo <- lapply(1:(ncol(tmpAdd)-1),function(thisCol){
									tmpReturn <- table(tmpAdd[,thisCol])
									return( unname(tmpReturn/length(unique(tmpAdd[,ncol(tmpAdd)]))) )
								})
				tmp_entropyInfo <- sapply(tmp_probInfo,calcEntropy,simplify=FALSE)
				names(tmp_entropyInfo) <- paste("entropy",c("genotypeID","binaryString","line_ofDescent","firstStep"),sep="_")
				
				entropy_plotInfo <- rbind(entropy_plotInfo,
											data.frame(tmp_entropyInfo,
														"Dilution" = all_parmInfo[which(all_parmInfo$jobName == thisJob),"distFactor"],
														"set"= as.numeric(strsplit(thisJob,"[[:punct:]]")[[1]][2]),
														"job"=as.numeric(strsplit(thisJob,"[[:punct:]]")[[1]][3]),
														"mutationRate"=unname(tmp_jobParms["const_mutProb"]),
														"landscapeModel"=unname(tmp_jobParms["simModel"]),
														"DFE"=paste(tmp_jobParms[1,c("constDist","const_distParameters")],collapse="__"),
														"jobName"=thisJob,
														stringsAsFactors = FALSE),
											make.row.names=FALSE)
			}
			# We can now remove the loaded object from memory
			rm(list=tmp_loadObj["name"])
			gc()
		}
		# I round the Dilution information for the plotInfo
		entropy_plotInfo[,"Dilution"] <- round(entropy_plotInfo[,"Dilution"],0)
		
		# Now we save the information for later use.
		save(entropy_plotInfo,file=paste(imageFilename["saveDir"],imageFilename["plotEntropy"],sep=""))
		
	} else {
		# This should hold an object called rank_plotInfo
		load(paste(imageFilename["saveDir"],imageFilename["plotEntropy"],sep=""))	
	}
	
	########### Barplots ###########
	# Now we can setup the plot which will simply be to report on the entropy meaures, as divided across parameter combinations
	tmp_plotData <- entropy_plotInfo
	tmp_plotData$Dilution <- factor(tmp_plotData$Dilution,levels=unique(tmp_plotData$Dilution[order(tmp_plotData$Dilution,decreasing = FALSE)]))
	tmp_plotData$mutationRate <- factor(tmp_plotData$mutationRate,levels=unique(tmp_plotData$mutationRate[order(tmp_plotData$mutationRate,decreasing = FALSE)]))
	
	entropy_lightAlpha <- 1
	entropy_plotCols <- colTransparency(rainbow(length(unique(tmp_plotData$Dilution)),start=0,end=4/6),entropy_lightAlpha)
	names(entropy_plotCols) <- as.character(levels(tmp_plotData$Dilution))
	
	# I'll be doing a new plot for each of the levels of mutations scored
	tmp_plotDims <- c("cols"= length(unique(tmp_plotData$landscapeModel)),
						"rows"= length(unique(tmp_plotData$DFE)))
	# Here is the range of x I want plotted, it's based on the levels of distFactor
	entropy_plotExpand = length(unique(tmp_plotData$mutationRate))
	entropy_xLim <- seq(1,length(unique(tmp_plotData$mutationRate)) * entropy_plotExpand,length.out = length(unique(tmp_plotData$mutationRate)))
	# I now define the plotSpace for my x-Axis
	entropy_plotSpots <- define_plotSpace(func_subDivs = tmp_plotData[,c("mutationRate","Dilution")], 
										func_tmpSpace = range(entropy_xLim),
										func_tmpBuffer = 0.05, 
										func_tmpBuffer_max = 0.2, 
										func_tmpBuffer_min=0.01, 
										func_orderDivs = TRUE)
	
	pdf(paste(outDir,"plot_entropy_hist.pdf",sep=""),width=8*tmp_plotDims[2],height=4*tmp_plotDims[1]+1)
	layout(matrix(c(seq(1, prod(tmp_plotDims)),
						rep(prod(tmp_plotDims)+1,tmp_plotDims[2])), 
				nrow=tmp_plotDims[1] + 1, byrow=TRUE),
			heights=c(rep(4,tmp_plotDims[1]),1))
	# Now this plot will be separated into different pages for the different entropy values I'm going to plot.
	plot_entropyNames <- c("entropy_genotypeID","entropy_line_ofDescent","entropy_firstStep")
	for(thisEntropy in plot_entropyNames){
		par(mai=c(rep(0.75,3),0.15))
		# We define the y limits given the maximum entropy value + a bit
		entropy_yLim <- c(0,decimalCeiling(max(tmp_plotData[,thisEntropy])*1.05,2))
	
		
		for(thisModel in unique(tmp_plotData$landscapeModel)){
			for(thisDFE in unique(tmp_plotData$DFE)){
				plot(x = NULL, y = NULL, xlim = c(min(entropy_xLim),max(entropy_xLim)), ylim= entropy_yLim,
					xlab="",ylab= "",
					axes = FALSE, cex.lab = 1.75, cex.axis = 1.5,
					main = paste(paste("Landscape Model: ",thisModel,sep=""),
								paste("DFE: ",thisDFE,sep=""),
								paste("Entropy Type: ",sub("entropy_","",thisEntropy),sep=""),
								sep="\n"))
				axis(side=1,at = sapply(entropy_plotSpots,function(x){ mean(unlist(x)) }), labels = levels(tmp_plotData$mutationRate),
					cex.axis = 1.75)
				mtext("Mutation Rate",side=1,line = 4, cex = 1.5)
				axis(side=2,at = sapply(seq(min(entropy_yLim),max(entropy_yLim),length.out=7),decimalCeiling,"func_level"=2),
					cex.axis = 1.75, las = 2)
				mtext("Entropy",side=2,line = 5, cex = 1.5)
				# Now we'll be plotting separate information for each mutation rate and Dilution factor associated with this set
				for(thisRate in levels(tmp_plotData$mutationRate)){
					for(thisDilution in levels(tmp_plotData$Dilution)){
						# We look for rows which associate to all the subsetting criteria and check that there is at least one such row
						tmpRows <- intersect(intersect(intersect(which(as.character(tmp_plotData$mutationRate) == thisRate),
																which(as.character(tmp_plotData$Dilution) == thisDilution)),
														which(tmp_plotData$landscapeModel == thisModel)),
											which(tmp_plotData$DFE == thisDFE))
						if(length(tmpRows) > 0){
							# We find the x positions based on the defined plot space call
							tmpX <- entropy_plotSpots[[which(levels(tmp_plotData$mutationRate) == thisRate)]][[which(levels(tmp_plotData$Dilution) == thisDilution)]]

							tmpY <- tmp_plotData[tmpRows,thisEntropy]
							polygon(x = c(tmpX, rev(tmpX)),
									y = c(rep(mean(tmpY) + 1.96 * sd(tmpY),2),
											rev(rep(mean(tmpY) - 1.96 * sd(tmpY),2))),
									col = colTransparency(entropy_plotCols[as.character(thisDilution)],plot_colourAlpha), 
									border = entropy_plotCols[as.character(thisDilution)])
							lines(x=tmpX,
									y = rep(mean(tmpY),2),
									lwd = 2, col = entropy_plotCols[as.character(thisDilution)],
									cex = 3)
						}
					}
				}
			}
		}
		# This will now be plotted in the extra space of the 
		par(mai=c(0,0,0,0))
		plot.new()
		legend("center",legend = names(entropy_plotCols),
				fill = entropy_plotCols,
				title = "Dilution Factor", bty = "n", horiz = TRUE,
				cex = 2)
	}
	dev.off()
	
	########### Lineplot ###########
	# Now we can setup the plot which will simply be to report on the entropy meaures, as divided across parameter combinations
	tmp_plotData <- entropy_plotInfo
	
	entropy_lightAlpha <- 0.4
	entropy_plotCols <- c("red3","cornflowerblue","purple")#colTransparency(rainbow(length(unique(tmp_plotData$mutationRate)),start=0,end=4/6))
	names(entropy_plotCols) <- as.character(unique(tmp_plotData$mutationRate))
	
	# I'll be doing a new plot for each of the levels of mutations scored
	tmp_plotDims <- c("cols"= length(unique(tmp_plotData$landscapeModel)),
						"rows"= length(unique(tmp_plotData$DFE)))
	
	pdf(paste(outDir,"plot_entropy_lines.pdf",sep=""),width=8*tmp_plotDims[2],height=4*tmp_plotDims[1]+1)
	layout(matrix(c(seq(1, prod(tmp_plotDims)),
						rep(prod(tmp_plotDims)+1,tmp_plotDims[2])), 
				nrow=tmp_plotDims[1] + 1, byrow=TRUE),
			heights=c(rep(4,tmp_plotDims[1]),1))
	# Now this plot will be separated into different pages for the different entropy values I'm going to plot.
	plot_entropyNames <- c("entropy_genotypeID","entropy_line_ofDescent","entropy_firstStep")
	entropy_xLim <- unique(tmp_plotData$Dilution)[order(unique(tmp_plotData$Dilution))]
	for(thisEntropy in plot_entropyNames){
		par(mai=c(rep(0.75,3),0.15))
		# We define the y limits given the maximum entropy value + a bit
		entropy_yLim <- c(0,decimalCeiling(max(tmp_plotData[,thisEntropy])*1.05,2))
	
		for(thisModel in unique(tmp_plotData$landscapeModel)){
			for(thisDFE in unique(tmp_plotData$DFE)){
				plot(x = NULL, y = NULL, xlim = c(min(log2(entropy_xLim)),max(log2(entropy_xLim))), ylim= entropy_yLim,
					xlab="",ylab= "",
					axes = FALSE, cex.lab = 1.75, cex.axis = 1.5,
					main = paste(paste("Landscape Model: ",thisModel,sep=""),
								paste("DFE: ",thisDFE,sep=""),
								paste("Entropy Type: ",sub("entropy_","",thisEntropy),sep=""),
								sep="\n"))
				axis(side=1,at = log2(entropy_xLim), labels = entropy_xLim,
					cex.axis = 1.75)
				mtext("Dilution Factor",side=1,line = 4, cex = 1.5)
				axis(side=2,at = sapply(seq(min(entropy_yLim),max(entropy_yLim),length.out=7),decimalCeiling,"func_level"=2),
					cex.axis = 1.75, las = 2)
				mtext("Entropy",side=2,line = 5, cex = 1.5)
				# Now we'll be plotting separate information for each mutation rate
				for(thisRate in unique(tmp_plotData$mutationRate)){
					# We look for rows which associate to all the subsetting criteria and check that there is at least one such row
					tmpRows <- intersect(intersect(which(as.character(tmp_plotData$mutationRate) == thisRate),
													which(tmp_plotData$landscapeModel == thisModel)),
										which(tmp_plotData$DFE == thisDFE))
					if(length(tmpRows) > 0){
						# We find the x positions based on the defined plot space call
						tmpX <- unique(tmp_plotData[tmpRows,"Dilution"])[order(unique(tmp_plotData[tmpRows,"Dilution"]))]
						# The y values will be the mean and SD of the rows, by dilution
						tmpY <- sapply(tmpX,function(x){ 
										tmp_subRows <- tmpRows[which(tmp_plotData[tmpRows,"Dilution"] == x)]
										return( c("mean"=mean(tmp_plotData[tmp_subRows,thisEntropy]),
													"sd"=sd(tmp_plotData[tmp_subRows,thisEntropy])) )
									})
						polygon(x = c(log2(tmpX), rev(log2(tmpX))),
								y = c(tmpY["mean",] + 1.96 * tmpY["sd",],
										rev(tmpY["mean",] - 1.96 * tmpY["sd",])),
								col = colTransparency(entropy_plotCols[as.character(thisRate)],entropy_lightAlpha), 
								border = entropy_plotCols[as.character(thisRate)])
						lines(x=log2(tmpX),
								y = tmpY["mean",],
								lwd = 2, col = entropy_plotCols[as.character(thisRate)],
								cex = 3)
					}
				}
			}
		}
		# This will now be plotted in the extra space of the 
		par(mai=c(0,0,0,0))
		plot.new()
		legend("center",legend = names(entropy_plotCols),
				fill = entropy_plotCols,
				title = "Mutation Rate", bty = "n", horiz = TRUE,
				cex = 2)
	}
	dev.off()
}



###################################################################################################
######################################## END OF PLOTTING ##########################################
###################################################################################################

# Quit without saving our session.
q(save="no")
