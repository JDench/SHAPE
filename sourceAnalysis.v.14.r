# This is a script which acts as the source information for functions to be called for processing bdMut_ series information
# that is being produced by another script.  Thus this contains nothing but information concerning the functions and is dependent
# upon another script for inheriting information such as file-pathing etc...

############ DEPENDENCIES ############
# This uses the output object created by SHAPE_body.v.28.r or higher and is called from within that run itself.
# This also depends on certain functions that are included in sourceFunctions.v.15.r that should be within the global environment.
# This also assumes the following objects exist within the workspace: allow_backMutations, max_numMutations, const_focal_popValue, collapseString, sepString, size_timeStep

# load libraries
library(abind) # This allows us to bind arrays together
# These would have been compiler, and RSQLite if not called within the context of a job run....

# This option will make certain that my data.frames are created as strings and not factors
options(stringsAsFactors = FALSE)
###################################################################################################
#################################### BEGIN OF PREDEFINITIONS ######################################
###################################################################################################
# This is a logical toggle as to whether we want this script to place its output in a newly named directory
# or to simply use the outDir object.
use_newDir <- FALSE
####### FILE AND SYSTEM LOCATIONS #######
postDir <- if(use_newDir){
				paste(workDir, "postAnal","/",sep="")
			} else {
				outDir
			}
if(!dir.exists(postDir)){ system(paste("mkdir -p ", postDir,sep="")) }

# These are filenames for information of collected batches of data
processedData_fileName <- paste(postDir,"processed_runData_from_", save_batchString,"_", thisRep,".RData",sep="")
processedObjects <- c("runDemographics","info_estLines")
# This is an object to control which strains we get deep neighbourhood information for
# It should be one of "none","limited","priority","full"
# setting this higher will cost more and more in post analysis runtime.
const_hoodDepth <- "limited"

############## CONSTANTS ################
# This is a string used to compress strings of information such as the evolutionary trajectory
string_lineDescent <- "_->_"
# This is the number of decimal places you want reported values to be rounded to
const_sigFig <- 4

###################################################################################################
##################################### END OF PREDEFINITIONS #######################################
###################################################################################################



###################################################################################################
##################################### BEGIN OF FUNCTIONS ##########################################
###################################################################################################
# this is a convenience script for building dataframe lines for our pedigree, it return a data.frame 
# which can be used to create or expand a pedigree data.frame.
buildPedigree <- function(func_focalID){
	# The pedigree frame is simply a list of lists
	return( sapply(as.character(func_focalID),function(thisID){  
					list()
				},simplify=FALSE) )
}

# This function looks through the func_pedigreeFrame (passed in full as func_pedigreeAll), and func_lineageDemographics (passed in full as func_demoArray) 
# which should include all parents that have ever been a parent to a lineage which needs parenthood, tracked with values >= 0 for as long as it's important 
# for focalID lineage definition!
findParent <- function(func_focalGenotype, func_startStep, func_stepMatrix, func_progenitorList, 
						func_demoArray, func_pedigreeAll){
	# We look for the first step at which this lineage did not exist, but earlier than the starting step 
	func_missingSteps <- rownames(func_stepMatrix)[which(func_stepMatrix[,"popSize"] == 0)]
	func_missingSteps <- as.numeric(nameTable_step(func_missingSteps,funcSplit = TRUE))
	# We keep only the steps which are after our start step
	func_missingSteps <- func_missingSteps[which(func_missingSteps < func_startStep)]
	# If there are no missing steps then this genotype has always existed, so we return it as the parent
	if(length(func_missingSteps) == 0){
		return( func_focalGenotype )
	# this means there is something to be done
	} else {
		# We look for the parent(s) of our genotype on the max step
		func_fociParents <- unlist(lapply(names(func_progenitorList),function(thisProgenitor){ 
								# If this porgenitor gave a mutant at this step we return it, else NULL
								if(is.element(max(func_missingSteps+1), func_progenitorList[[thisProgenitor]])){
									return( thisProgenitor )
								} else {
									return( NULL )
								}
							}))
		# Then we pass the max step which is missing the foci along with the parent(s) down the recursion
		# This assumes our func_stepMatrix and func_progenitorList will always have the information for all
		# parental types
		return( paste(unlist(lapply(func_fociParents,function(thisProgenitor){
									findParent(func_focalGenotype = as.numeric(thisProgenitor),
													func_startStep = max(func_missingSteps),
													func_stepMatrix = func_demoArray[,, thisProgenitor],
													func_progenitorList = func_pedigreeAll[[thisProgenitor]],
													func_demoArray = func_demoArray,
													func_pedigreeAll = func_pedigreeAll)
						})), 
						func_focalGenotype,
						sep= string_lineDescent) )
	}

}

# This is a function that steps forward through time and extracts demongraphic information at the level of the population
# It will extract information such as the Fitness, Number of Lineages, and Transitions between dominant genotypes.
# Most important it will also return the information related to which lineages will eventually establish in the population
# which as a piece of information will be critical for lineage specific information extraction.
extract_popDemographics <- function(func_stepsCon, func_estValue, func_landscapeCon, func_hoodCon){
	# We find the ordered series of step tables that we'll be calling in this extaction process
	func_allTables <- unlist(dbListTables(func_stepsCon))
	# We now order that set to an ascending series
	func_allTables <- func_allTables[order(as.numeric(nameTable_step(func_allTables,funcSplit = TRUE)), decreasing = FALSE)]
	# Now we establish some reporting lists the first is the general fitness and the next the number of lineages information
	# the next being for the information regarding transition timing, the last for tracking any lineage that is every established
	tmpReturn <- list("demoMat"= matrix(-1,nrow = length(func_allTables), ncol = 8, 
										dimnames = list(func_allTables,c("minFit","meanFit","maxFit","sdFit","numLines","numEstablished","numMutants","popSize"))),
						"transitionMat" = matrix(-1,ncol=5,dimnames=list(NULL,c("Step","genotypeID","numMuts","fitness","transitionGens"))),
						"vec_estLineages" = 0,
						"vec_final_estLineages" = 0,
						"Hists" = vector(mode="list",length=length(func_allTables)))
	# We set the names of the histograms
	names(tmpReturn[["Hists"]]) <- func_allTables
	
	# We now go forward through each of the time steps and we extract information
	for(thisStep in func_allTables){
		# We go through each of the steps, extracting the pertinent information
		tmpData <- dbGetQuery(func_stepsCon, paste("SELECT numMuts, genotypeID, popSize, fitness, mutants FROM ",thisStep,sep=""))
		# Any genotype with fewer than 0 individuals will be recorded as having effectively zero
		tmpData[which(tmpData[,"popSize"] <= 0),"popSize"] <- 1e-1
		# We extract all established lineages 
		tmp_estLineages <- querryEstablished(func_inMatrix = tmpData, func_estProp = func_estValue)
		# If there are any estblished line(s) we will be updating this information
		if(nrow(tmp_estLineages) > 0){
			# We now keep only those which are no already in the vec_estLineages																					
			tmpReturn[["vec_estLineages"]] <- unique(c(tmpReturn[["vec_estLineages"]], tmp_estLineages$genotypeID))
			# If we're on the final step, we'll record which lineages were established at this point.
			if(thisStep == func_allTables[length(func_allTables)]){
				tmpReturn$vec_final_estLineages <- tmp_estLineages$genotypeID
			}
		}
		
		# Now we check if there has been a transition between dominant genotypes, this is recorded by looking
		# if there are any dominant lineages, and thereafter by looking at the last row
		tmp_maxLine <- tmpData[which(tmpData$popSize == max(tmpData$popSize)),]
		# We find if there are any established lineage(s) not previously recorded at the same step, meaning we need to find which
		# rows of our transition matrix relate to the last time step(s)
		tmp_transRow <- which(tmpReturn[["transitionMat"]][,"Step"] == max(tmpReturn[["transitionMat"]][,"Step"]))
		# If any of the max sized lineage(s) in this step is not the same as the current max then we update our information about the max lineages in thisStep
		if(any(!is.element(tmp_maxLine$genotypeID,tmpReturn[["transitionMat"]][tmp_transRow,"genotypeID"]))){
			tmpReturn[["transitionMat"]] <- rbind(tmpReturn[["transitionMat"]],
													matrix(c(rep(as.numeric(nameTable_step(thisStep,funcSplit = TRUE)),nrow(tmp_maxLine)),
																tmp_maxLine$genotypeID,
																tmp_maxLine$numMuts,
																tmp_maxLine$fitness,
																rep(as.numeric(nameTable_step(thisStep,funcSplit = TRUE)) * size_timeStep - 
																	ifelse(tmpReturn[["transitionMat"]][nrow(tmpReturn[["transitionMat"]]),"Step"] != -1,tmpReturn[["transitionMat"]][nrow(tmpReturn[["transitionMat"]]),"Step"],0)* size_timeStep,
																	nrow(tmp_maxLine)) ),
															nrow=nrow(tmp_maxLine)))
		}	
		
		# We insert the information into our demoFrame
		tmpReturn[["demoMat"]][thisStep,] <- c(min(tmpData$fitness),
												weighted.mean(tmpData$fitness,tmpData$popSize),
												max(tmpData$fitness),
												if(nrow(tmpData) == 1){ 0 } else { sd(tmpData[,"fitness"]* tmpData[,"popSize"]/sum(tmpData[,"popSize"])) },
												nrow(tmpData),
												nrow(tmp_estLineages),
												sum(tmpData$mutants),
												sum(tmpData$popSize))
												
		# We now create a histogram of the population fitness, with counts considering popSize
		# I actually do the fitness histogram in a "cheaty" fashion 
		tmp_fitHist <- hist(tmpData$fitness,plot=FALSE)
		# Now using the counts and breaks that exist (where breaks is always +1 longer...) we sum the popSize
		# Noting that the lower bound of the first instance of breaks is included but is exclusive for the others.
		tmp_fitHist$counts <- log10(if(length(tmp_fitHist$counts) > 1){
									sapply(1:length(tmp_fitHist$counts),function(x){ 
										if(x == 1){
											sum(tmpData$popSize[intersect(which(tmpData$fitness >= tmp_fitHist$breaks[x]),
																			which(tmpData$fitness <= tmp_fitHist$breaks[x+1]))])
										} else {
											sum(tmpData$popSize[intersect(which(tmpData$fitness > tmp_fitHist$breaks[x]),
																			which(tmpData$fitness <= tmp_fitHist$breaks[x+1]))])
										}
									})
								# This means there is one fitness bin so all popSize are the count value...	
								} else {
									sum(tmpData$popSize)
								})
		# Now we build the histogram of population sizes
		tmp_sizeHist <- hist(log10(tmpData$popSize),plot = FALSE)
		# then update the counts to be a matter of log2
		tmp_sizeHist$counts <- log2(tmp_sizeHist$counts)
		# We remove any infinite values from the counts of either histogram
		tmp_fitHist$counts[which(tmp_sizeHist$counts == "-Inf")] <- -1
		tmp_sizeHist$counts[which(tmp_sizeHist$counts == "-Inf")] <- -1
			
		tmpReturn[["Hists"]][[thisStep]] <- list("fitness"= tmp_fitHist,
												"lines"= tmp_sizeHist)
		
	}
	# We just drop that initialising row for transitionMat
	tmpReturn[["transitionMat"]] <- matrix(tmpReturn[["transitionMat"]][-1,],ncol=ncol(tmpReturn[["transitionMat"]]), dimnames= dimnames(tmpReturn[["transitionMat"]]))
	
	# We now return this information
	return( tmpReturn )
}

# This is a function to extract the lineage specific information.  This info will be mostly through time style of information
# but will also include information about it's line of descent, growth pressures pre-establishment, and popSize.  The refMatrix 
# will be querried for information regarding the numMuts and fitness values of genotypes, but is not required but is also required 
# to know what are transitions which occure, hence why it too must be supplied.  It should be the runDemographics[["transitionMat"]] object.
# We include the argument of hoodExplore with possible values of "none", "limited","priority", and "full".
# NOTE: That use of limited requires that you pass a func_refMatrix of expected shape (has a "genotypeID" column)!
extractInfo_focalID <- function(func_focalID, func_estValue, func_stepsCon, func_landscapeCon,func_hoodCon, func_refMatrix,
								func_descentSep = string_lineDescent, func_hoodExplore = "priority"){
	# We find the ordered series of step tables that we'll be calling in this extaction process
	func_allTables <- unlist(dbListTables(func_stepsCon))
	# We now order that set to a decreasing series
	func_allTables <- func_allTables[order(as.numeric(nameTable_step(func_allTables,funcSplit = TRUE)), decreasing = TRUE)]
	# I also order the focalID's
	func_focalID <- func_focalID[order(func_focalID)]
	
	# We establish the reporting information objects 
	func_lineageDemographics <- array(-1,dim=c(length(func_allTables),5,length(func_focalID)),
										dimnames=list(func_allTables,c("Step","popSize","isEstablished","births","mutsIn"),as.character(func_focalID)))
	func_pedigreeFrame <- buildPedigree(func_focalID[order(func_focalID)])
	
	# Now we simply step back through each of the time steps and grab information pertinent to our focal lineages
	for(thisTable in func_allTables){
		# I start by getting the table index, this will be used when we need to pass information forward or backward
		tmpIndex <- which(func_allTables == thisTable)
		# We load the step information table
		tmpData <- dbGetQuery(func_stepsCon, paste("SELECT * FROM ",thisTable," WHERE genotypeID IN (",
														paste(unique(names(func_pedigreeFrame)),collapse=","),
														")",sep=""))
		# I also querry which line(s) are established
		tmp_estLines <- querryEstablished(func_inMatrix = tmpData, func_estProp = func_estValue)
		
		# My next step should be to record the popSize of that lineage in this step, is it established, any births,
		# but also, and this is the trickier part, the information concerning number of mutants which comes from progenitor info
		for(func_thisMatrix in unique(names(func_pedigreeFrame))){ 
			func_dataRow <- which(tmpData[,"genotypeID"] == func_thisMatrix)
			# If there is some information but this genotype is not yet in our matrix then we need update our array
			if(!is.element(func_thisMatrix, dimnames(func_lineageDemographics)[[3]])){
				func_lineageDemographics <- abind(func_lineageDemographics, 
													array(-1,dim=c(length(func_allTables),5,length(func_thisMatrix)),
													dimnames=list(func_allTables,c("Step","popSize","isEstablished","births","mutsIn"), func_thisMatrix)) )
	
			}
			# Now we actually update our demographics but if the lineage does not exist in this step we return zero values.
			if(length(func_dataRow) > 0){
				# I pre-sum the number of mutants which were fed into this lineage as I'll be subtracting those from births so these values are separated
				func_tmpSum <- sum(unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], collapseString)[[1]],function(x){
																				as.numeric(strsplit(x,sepString)[[1]][2])
																		})),na.rm=TRUE)
				func_lineageDemographics[thisTable,,as.character(func_thisMatrix)] <- c(as.numeric(nameTable_step(thisTable,funcSplit = TRUE)),
																						tmpData[func_dataRow,"popSize"],
																						is.element(func_thisMatrix, tmp_estLines[,"genotypeID"]),
																						tmpData[func_dataRow,"births"] - func_tmpSum,
																						func_tmpSum)
			} else {
				func_lineageDemographics[thisTable,,as.character(func_thisMatrix)] <- c(as.numeric(nameTable_step(thisTable,funcSplit = TRUE)),
																							rep(0,4))
			}
			# Now we can work on the pedigree for this lineage, we start by checking if there are any parents for the main focal lineage
			# But only if there are rows to work with
			if(length(func_dataRow) > 0){
				# We have parents if the progenitor column does not contain a blank or WT value
				if(!is.element(tmpData[func_dataRow,"progenitor"],c("","WT"))){
					# These are the genotypeIDs which contributed to our focalID
					func_tmpParents <- unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], collapseString)[[1]],function(x){
																				tmpReturn <- as.numeric(strsplit(x,sepString)[[1]][1])
																				# Now because a lineage may exist and receive new mutants this will create
																				# situations where we get NA values returned so we return null for those parents
																				if(is.na(tmpReturn) || is.null(tmpReturn)){
																					return( NULL )
																				} else {
																					return( tmpReturn )
																				}
																			}))
					# Lastly we see if we need to update the pedigree data.frame for this focalID by querrying which of it's parents are not yet represented
					func_missingParents <- func_tmpParents[which(!is.element(func_tmpParents, names(func_pedigreeFrame)))]
					if(length(func_missingParents) > 0){
						func_pedigreeFrame <- c(func_pedigreeFrame,
												buildPedigree(func_focalID = func_missingParents)) 
					}
				}
			}
		}
		# Now we update our tmpData records 
		tmpData <- dbGetQuery(func_stepsCon, paste("SELECT * FROM ",thisTable," WHERE genotypeID IN (",
														paste(unique(names(func_pedigreeFrame)),collapse=","),
														")",sep=""))
		# Now for each of the focal lineages' pedigreeID's we ask ...
		for(func_thisID in unique(names(func_pedigreeFrame))){
			func_dataRow <- which(tmpData[,"genotypeID"] == as.numeric(func_thisID))
			# If we can find the row for thisID, then we want to update it provided we find some "parents"
			if (length(func_dataRow) > 0 && !is.element(tmpData[func_dataRow,"progenitor"],c("","WT"))){
				# Then we want to update the lineage of this genotype, this requires us to first define parents
				func_tmpParents <- unlist(lapply(strsplit(tmpData[func_dataRow,"progenitor"], collapseString)[[1]],function(x){
																			tmpReturn <- as.numeric(strsplit(x,sepString)[[1]][1])
																			# Now because a lineage may exist and receive new mutants this will create
																			# situations where we get NA values returned so we return null for those parents
																			if(is.na(tmpReturn) || is.null(tmpReturn)){
																				return( NULL )
																			} else {
																				return( tmpReturn )
																			}
																		}))													
				# Right, then we're going to store parents is all the steps at which they provide mutant(s), 
				# this is a convenience to reduce the length of lists, but makes later extractions perhaps more involved...
				# For each of the parents, we check if they're in the the list of thisID, and if not add them
				for(func_thisParent in as.character(func_tmpParents)){
					func_recordStep <- as.numeric(nameTable_step(thisTable,funcSplit = TRUE))
					if(is.element(func_thisParent,names(func_pedigreeFrame[[func_thisID]]))){
						func_pedigreeFrame[[func_thisID]][[func_thisParent]] <- c(func_pedigreeFrame[[func_thisID]][[func_thisParent]], func_recordStep)
					} else {
						func_pedigreeFrame[[func_thisID]][[func_thisParent]] <- func_recordStep
					}	
				}
			} # This is the conditional that there is something to update and that a progenitor is one of those things...
		} # This closes the loop of looking for focalID genotypes
	} # This closes out the for loop for the time step tables										
	
	# I define which genotypes existed on the first step, this is so that they can be ignored when they appear as root genotypes. 
	func_startGenotypes <- unname(unlist(dbGetQuery(func_stepsCon, paste("SELECT genotypeID FROM ",nameTable_step(0),sep=""))))
	# Right, now what we want is the lineage for all focalID genotypes which exist on the last step AND 
	# any lineages which were focal but not a first genotype.  These last would have been 
	# focalID's and so must be in the pedigree information already.
	func_endFoci <- intersect(func_focalID, unname(unlist(dbGetQuery(func_stepsCon, paste("SELECT genotypeID FROM ", func_allTables[1],sep="")))))
	
	func_nonEndFoci <- setdiff(func_focalID, c(func_startGenotypes, func_endFoci))
	# This calls our function of finding a parent, is calls itself recursively to work through the matrix.
	func_endLineages <- list()
	for(thisFoci in func_endFoci){
		func_endLineages[[as.character(thisFoci)]] <- findParent(func_focalGenotype = thisFoci,
																func_startStep = as.numeric(nameTable_step(rownames(func_lineageDemographics)[which(func_lineageDemographics[,"popSize",as.character(thisFoci)] > 0)[1]],funcSplit = TRUE)),   
																func_stepMatrix = func_lineageDemographics[,,as.character(thisFoci)],
																func_progenitorList = func_pedigreeFrame[[as.character(thisFoci)]],
																func_demoArray = func_lineageDemographics,
																func_pedigreeAll = func_pedigreeFrame)
	}
	func_nonendLineages <- list()
	if(length(func_nonEndFoci) > 0){
		for(thisFoci in func_nonEndFoci){
			func_nonendLineages[[as.character(thisFoci)]] <- findParent(func_focalGenotype = thisFoci,
																	func_startStep = as.numeric(nameTable_step(rownames(func_lineageDemographics)[which(func_lineageDemographics[,"popSize",as.character(thisFoci)] > 0)[1]],funcSplit = TRUE)),   
																	func_stepMatrix = func_lineageDemographics[,,as.character(thisFoci)],
																	func_progenitorList = func_pedigreeFrame[[as.character(thisFoci)]],
																	func_demoArray = func_lineageDemographics,
																	func_pedigreeAll = func_pedigreeFrame)
		}
	}
	
	# We need to know which fitness landscape tables exist 
	func_landTables <- dbListTables(func_landscapeCon)
	# Great, now we start building a matrix which displays neighbourhood information for all unique transitions.
	# We find a unique transition by looking at the lineages at all intermediate positions between two values.
	func_uniqueTransitions <- unique(unlist(lapply(strsplit(as.character(unique(unlist(c(func_endLineages, func_nonendLineages)))), string_lineDescent),function(thisLineage){
									# If there is only a single element, there has been no transition...
									if(length(thisLineage) > 1){
										return( paste(thisLineage[-length(thisLineage)], thisLineage[-1],sep= string_lineDescent) )
									} else {
										NULL
									}
								})))
	# Before going to get he landscape topology information I'd like to extract the numMuts and fitness value(s) for unique
	# genotypes that are involved in transitions. 
	func_uniquetransitionIDs <- unique(unlist(strsplit(as.character(unique(unlist(c(func_endLineages, func_nonendLineages)))), string_lineDescent)))
	func_refInfo <- matrix(c(as.numeric(func_uniquetransitionIDs),rep(-1,length(func_uniquetransitionIDs) *2)),
							ncol=3,dimnames=list(func_uniquetransitionIDs,c("genotypeID","numMuts","fitness")))
	# Ok, now if we've got a func_refMatrix object lets populate what we can into our func_refInfo
	if(!is.null(func_refMatrix)){
		# we define the unique genotypeID in the ref object
		func_tmpIDs <- intersect(unique(func_refMatrix[,"genotypeID"]), as.numeric(func_uniquetransitionIDs))
		func_refInfo[as.character(func_tmpIDs),c("numMuts","fitness")] <- func_refMatrix[unlist(lapply(func_tmpIDs,function(x){ 
																										which(func_refMatrix[,"genotypeID"] == x)[1] 
																									})),
																							c("numMuts","fitness")]
	}
	# Next if we're missing the fitness OR numMuts for any of these genotypeIDs we'll go and grab that from the fitness landscape
	func_missIDs <- unname(unlist(apply(func_refInfo,MARGIN=1,function(thisLine){
							if(any(thisLine[c("genotypeID","numMuts","fitness")] == -1)){
								return( thisLine["genotypeID"] )
							} else {
								return( NULL )
							}		
						})))
	# If we're not missing the information for any IDs then we can skip the following set of calls
	if(length(func_missIDs) > 0){
		# So for each of these ID's we go and query for the fitness and or numMuts information as required
		# NOTE: numMuts will be informed by binaryString length after strsplit
		# I make a single database query which may mean I'm replicating information but such is life...
		func_updateInformation <- dbGetQuery(func_landscapeCon, paste("SELECT genotypeID,binaryString,fitness FROM ", 
													func_landTables,
													' WHERE genotypeID IN(', 
													paste(func_missIDs,collapse=","),
													')',
												collapse=" UNION "))
		# Now we update the binaryString information to be genome length
		func_updateInformation[,"binaryString"] <- unlist(lapply(strsplit(func_updateInformation[,"binaryString"],sepString),length))
		colnames(func_updateInformation)[which(colnames(func_updateInformation) == "binaryString")] <- "numMuts"
		# We now update the information for any genotypes that were missing information
		for(thisCol in c("numMuts","fitness")){
			func_refInfo[as.character(func_missIDs[order(func_missIDs)]),thisCol] <- func_updateInformation[order(func_updateInformation[,"genotypeID"]), thisCol]
	
		}
	}

	# So we now look at the neighbourhood of the parental types and find where the fitness of the offspring fits within this
	# That will tell us what was the rank of the transition.  I won't gather information about the extent of neighbourhood 
	# exploration, I have no distinct use for this information and its calculation is costly with respect to time.
	func_hoodTables <- dbListTables(func_hoodCon)
	# If there have been no transitions we return a NULL value, else the data.frame expected
	func_rankMat <- if(is.null(func_uniqueTransitions)){
						NULL
					} else {
						t(sapply(func_uniqueTransitions,function(thisTransition){
							# The first element is the progenitor while the second element is the offspring.
							# We grab the numMuts and fitness information for each from our func_refInfo
							tmp_transitionInfo = func_refInfo[strsplit(thisTransition, string_lineDescent)[[1]],]
							# We gather the mutatation range of the parental type
							tmp_mutsRange <-  max(0,(tmp_transitionInfo[1,"numMuts"] - if(allow_backMutations){ max_numMutations }else{ 0 })): 
												min(genomeLength ,(tmp_transitionInfo[1,"numMuts"] + (max_numMutations)))
							# We record the binary stirng information of the parental genotype so that it is included within the search
							tmp_parent_binaryString <- retrieve_binaryString(func_genotypeID = tmp_transitionInfo[1,"genotypeID"], 
																			func_numMuts = tmp_transitionInfo[1,"numMuts"],
																			func_landscapeCon = func_landscapeCon)[1,"binaryString"]
							# We now find the neighbours of the parental type
							tmp_parentalHood <- unique(c(tmp_parent_binaryString,
														if(is.element(nameTable_neighbourhood(tmp_transitionInfo[1,"genotypeID"]), func_hoodTables)){
															dbGetQuery(func_hoodCon,paste("SELECT * FROM ",nameTable_neighbourhood(tmp_transitionInfo[1,"genotypeID"]),sep=""))$neighbours
														} else {
															# We define all the possible nearestNeighbours for this lineage since they've not been stored.
															# This includes calling for the binary string of this genotype
															comp_defineNeighbours(func_tmpGenotype = tmp_parent_binaryString, 
																					func_tmpDirection = allow_backMutations)
														}))
							# Knowing the neighbours we can gather the fitness of each, this is done as we have their binary strings
							tmp_tmpStrings <- gsub("[[:space:]]","",paste("\'", tmp_parentalHood,"\'",collapse=','))
							tmp_querryTables <- func_landTables[unique(unlist(lapply(nameTable(tmp_mutsRange),function(tmpTable){ 
																							which(grepl(tmpTable, func_landTables)) 
																						})))]
							tmp_hoodFit <- unlist(lapply(tmp_querryTables,function(thisTable){
																unlist(dbGetQuery(func_landscapeCon, paste('SELECT fitness FROM ',
																thisTable, 
																' WHERE binaryString IN (', 
																tmp_tmpStrings,
																')',
																sep=""))) 
														}), use.names=FALSE)
							tmp_hoodFit <- tmp_hoodFit[order(tmp_hoodFit, decreasing = TRUE)]
							tmp_childRank <- which(tmp_hoodFit == tmp_transitionInfo[2,"fitness"])[1]
							tmp_altPaths <- which(tmp_hoodFit == tmp_transitionInfo[1,"fitness"])[1] - 1
										### REMOVED: That way if there are problems they ought to be highlighted....
										# This is an unnecessary sanity check in-case somehow the offspring was not found in the parental hood
										# Which ought not to be possible but it doesn't affect the interpretation of functional transitions.
										#tmp_hoodFit <- c(tmp_transitionInfo[2,"fitness"], tmp_hoodFit)
										
							return(	c("absRank"= tmp_childRank,
										"hoodSize"=length(tmp_parentalHood),
										"hoodMin"=min(tmp_hoodFit),
										"hoodMax"=max(tmp_hoodFit),
										"progenitor_numMuts"= tmp_transitionInfo[1,"numMuts"],
										"progenitor_fitness"= tmp_transitionInfo[1,"fitness"],
										"offspring_numMuts"= tmp_transitionInfo[2,"numMuts"],
										"offspring_fitness"= tmp_transitionInfo[2,"fitness"],
										"progenitorID"=tmp_transitionInfo[1,"genotypeID"],
										"offspringID"=tmp_transitionInfo[2,"genotypeID"],
										"num_altPaths"=  tmp_altPaths,
										"relFit_altPaths" = if(tmp_altPaths > 0 && tmp_altPaths >= tmp_childRank){
																calc_relativeFitness(tmp_hoodFit[1:tmp_altPaths])[tmp_childRank] 
															} else {
																0
															},
										"prop_maxFit" = if(tmp_altPaths > 0 && tmp_altPaths >= tmp_childRank){
															(tmp_altPaths - (tmp_childRank -1)) / tmp_altPaths 
														} else {
															0
														}) )	
										
										
						}))
					}
	# We can now return our information
	return( list("lineDemo"=func_lineageDemographics,
				"linePedigree"=func_pedigreeFrame,
				"landscapeTopology"=func_rankMat,
				"end_Lines_of_Descent"=func_endLineages,
				"transition_Lines_of_Descent"=func_nonendLineages) )
}

# This is a function to actually perform the processing, it expects to be informed of the connection to a Steps file
# the proportion of total population that constitutes establishment (ie that we care to track it)
# NOTE: The default values placed here are based upon expected values from my script which generates the data, these can change as required.
runProcessing <- function(func_stepsCon = connections_dataBase[["timeStep_States"]], func_landscapeCon = connections_dataBase[["genotypeSpace"]],
							func_hoodCon = connections_dataBase[["nearestNeighbours"]],
							func_estProp = const_estProp,  func_saveFile = processedData_fileName, func_processObjects = processedObjects,
							func_hoodPriority = "priority"){
	# We check to see if there already exists a file which would have the same name as the output from this file, if not then we process.
	if(file.exists(func_saveFile)){
		return( paste("Found a pre-processed file of same expected name: ", func_saveFile," did not perform processing of run",sep="") )
	} else {
		# For each of these sapply calls, we'll nest a sapply call if a particular run was replicated
		
		for(thisObject in processedObjects){ assign(thisObject,NULL,pos=".GlobalEnv") }
		# Next is to extract demographic data from each of our sets
		runDemographics <- extract_popDemographics(func_stepsCon = func_stepsCon, func_estValue = func_estProp,
													func_landscapeCon = func_landscapeCon, func_hoodCon = func_hoodCon)
		# Now we're extracting information regarding lines which have established in our population
		info_estLines <- extractInfo_focalID(func_focalID = runDemographics[["vec_estLineages"]],
												func_estValue = func_estProp,
												func_stepsCon = func_stepsCon,
												func_landscapeCon = func_landscapeCon,
												func_hoodCon = func_hoodCon,
												func_refMatrix = runDemographics[["transitionMat"]],
												func_descentSep = string_lineDescent,
												func_hoodExplore = func_hoodPriority)
		# We now save our objects
		save(list= func_processObjects, file= func_saveFile)
		# We report having completed
		return( paste("Processing completed for :", func_saveFile,sep="") )
	}
}

###################################################################################################
###################################### END OF FUNCTIONS ###########################################
###################################################################################################

