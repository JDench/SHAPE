# This script is used to build a complete experiment with SHAPE.  
# It takes the parameters defined in SHAPE_farmerParms.v.#.ranges, then builds all directories and necessary files as outlined.
# Running this script 

rm(list=ls())

# This should be the explicit filepath for the directory in which all the template files of SHAPE exist.
### NOTE: The last character of this must be a forward slash  "/" 
templateDir <- "C:/Users/Jonathan/Documents/Programming/MyScripts/SHAPE/" #"/global/home/hpc3058/jonathan/SHAPE/"

# This should be an integer value which defines the maximum number of independent jobs you want to have bundled in a single submission script
# It applies to both local and remote server applications.  Setting this to 1 means you'll have to indpendently submit each job by running the .sh
# while setting this too high may risk overwhelming your computational resources.  
maxJobs <- 3

# This is an optinal string which can be used to uniquely define your parameter input file name, leave as NULL if you've not renamed it
input_fileName <- NULL

# This is the explicit filepath for the R installation, it can be for the remote server or your local installation.
appLocation <- "/global/home/hpc3058/bin/R-3.4.2/bin/R"
# these are common arguments passed to the R call for executing a batch run (unless R changes this can remain as is).
commonArgs <- "CMD BATCH"


# This is information concerning a remote server's submission system.  I'm basing this off of the SLURM system of the Center for Advanced Computing
# Queen's University computing platform.  If your system is different you may need to tweak this.  Sorry?

# This should be a vector of arguments passed for job submissions on a remote server 
# The example here would call 1 core with 8 Gb RAM and a wall time of 14 days and an outFile be named
# You can add more arguments if your server requires this, they'll get used.
### BUT where the job's name MUST be identified as ---  fakeJob   ----  and the output log as  --- fakeOut ---, you can change the argument queues
submitArgs <- c("number_ofCores"='-c 1',
				"memory"='--mem=8192',
				"jobName"='-J fakeJob',
				"wallTime"='-t 14-00:00:00',
				"fileOut"='-o fakeOut')
# I also assume your remote server will create a local directory for your job once submitted, and that there will be some environmental
# parameter you can call to ID that file path, place it here
remoteLocation <- "$TMPDISK"



###################################################################################################
########################## YOU LIKELY DON'T NEED TO MAKE CHANGES BELOW THIS LINE ##################
###################################################################################################

# These are arguments that SHAPE requires to be passed,  see the SHAPE_body file for details if you care.
passedArgs <- '"--args currReplicate=1 outDir=\'fake_serverPath/fakeDir/\'"'
# This is a string used as the name for a trigger queue to stop external selfing
external_stopFile <- "stop_trigger.file"


###################################################################
############################### BODY ##############################
###################################################################
# This is a reference object that links that farmerParms input names to the object names used in SHAPE
inputReference <- list("serverFarm"=							"serverFarm <-",
						"results_removeSteps"=					"results_removeSteps <-",
						"externalSelfing"=						"externalSelfing <-",
						"toggle_forceCompletion"=				"toggle_forceCompletion <-",
						"baseDir_local"=						c("baseDir <-",'else if(!serverFarm)','\"[^\"]*\"'),
						"baseDir_server"=						c("baseDir <-",'if(serverFarm)','\"[^\"]*\"'),
						"workDir"=								c("workDir <-",'paste(baseDir','\"[^\"]*\"'),
						"save_batchBase"=						"save_batchBase <-",
						"fileName_functions"=					c("fileName_functions <-",'\"[^\"]*\"'),
						"fileName_preProcessing"=				c("fileName_preProcessing <-",'\"[^\"]*\"'),
						"maxReplicates"=						"maxReplicates <-",
						"uniqueReplicates"=						"save_batchSet <-",
						"disturbanceType"=						"const_distType <-",
						"disturbanceSize"=						c("init_distPars <-",'\"factor\"=','([[:digit:]])+'),
						"disturbanceSpread"=					c("init_distPars <-",'\"random\"=','([[:digit:]])+'),
						"generations_betweenDisturbances"=		"const_growthGenerations <-",
						"numGenerations"=						"numGenerations <-",
						"genomeLength"=							"genomeLength <-",
						"targetNumber"=							"const_focal_popValue <-",
						"probabilityBirth"=						"const_birthProb <-",
						"growthRate"=							"const_growthRate <-",
						"probabilityDeath"=						"const_deathProb <-",
						"death_byDensity"=						"death_byDensity <-",
						"death_densityFactor"=					"death_densityCorrelation <-",
						"deathRate_sizeCap"=					"death_densityCap <-",
						"scale_births_by_deaths"=				"scaleGrowth_byDeaths <-",
						"mutationRate"=							"const_mutProb <-",
						"mutations_only_at_birth"=				"muts_onlyBirths <-",
						"allow_backMutations"=					"allow_backMutations <-",
						"growthModel"=							"const_growthForm <-",
						"simModel"=								"simModel <-",
						"const_ancestFitness"=					"const_ancestFitness <-",
						"constDist"=							"constDist <-",
						"const_distParameters"=					"const_distParameters <-",
						"const_distAsS"=						"const_distAsS <-",
						"const_RMF_initiDistance"=				"const_RMF_initiDistance <-",
						"const_RMF_theta"=						"const_RMF_theta <-",
						"const_numInteractions"=				"const_numInteractions <-",
						"const_fixedFrame"=						"const_fixedFrame <-",
						"stochasticBirths"=						"includeDrift <-",
						"track_only_integer_individuals"=		"track_asWhole <-",
						"size_ofEstablishment"=					"const_estProp <-",
						"trackingThreshhold" =					"const_hoodThresh <-")				
					
# This is a reference of the input parameters that are not used in making factorial combinations for an experiment
notExpanded_reference <- c("serverFarm",
							"results_removeSteps",
							"externalSelfing",
							"toggle_forceCompletion",
							"baseDir_local",
							"baseDir_server",
							"workDir",
							"save_batchBase",
							"fileName_functions",
							"fileName_preProcessing",
							"maxReplicates",
							"uniqueReplicates")
# This is a reference of those run parameters that used 1:1 to make a parameter set.
comboReference <- c("simModel",
					"const_ancestFitness",
					"constDist",
					"const_distParameters",
					"const_distAsS")
# This is a reference of parameters that only need to be made factorial given particular conditions
conditionReference <- list("RMF"=c("const_RMF_theta",
									"const_RMF_initiDistance"),
							"NK"=c("const_numInteractions"),
							"Fixed"=c("const_fixedFrame"))

# This is a standard critical error response
stopError <- function(func_message){
		print(func_message)
		Sys.sleep(20)
		q(save="no")
}
							
# These are the string literal expressions that can help this script find the appropriate templates in the templateDir folder
fileName_templates <- c("functions"="sourceFunctions.v.",
						"body"="SHAPE_body.v.",
						"parameter_input"=if(is.null(input_fileName)){"SHAPE_farmerParms"}else{input_fileName},
						"parameter_output"="SHAPE_parameters.v.",
						"plotting"="SHAPE_plotting.v.",
						"run_reset"="SHAPE_farmingReset.v.",
						"serverSubmit"="SHAPE_serverTemplate.v.",
						"localSubmit"="SHAPE_localTemplate.v.")
# Step 1: find the template files and report any that are missing
allFiles <- list.files(path = templateDir, recursive = FALSE)
for(thisFile in names(fileName_templates)){
	tmpFile <- allFiles[which(grepl(fileName_templates[thisFile], allFiles,fixed=TRUE))]
	# We check if there is a unique file 
	if(length(tmpFile) == 1){
		fileName_templates[thisFile] <- paste(templateDir,tmpFile,sep="")
		if(thisFile == "functions"){ source(fileName_templates[thisFile]) }
	} else {
		stopError(paste("Could not find a unique file like ",fileName_templates[thisFile]," in ",templateDir," please review and ensure the template exists",sep=""))
	}
}

# Step 2: read the input parameters and begin by building the experiment's folders
inputParms <- readLines(fileName_templates["parameter_input"])
inputParms <- inputParms[intersect(intersect(which(!grepl("#",inputParms,fixed=TRUE)),
												which(!grepl("[[:space:]]",substr(inputParms,
																					sapply(nchar(inputParms),function(x){ min(1,x) }),
																					sapply(nchar(inputParms),function(x){ min(1,x) })) ))),
												which(nchar(inputParms) > 0))]
# We now parse the input into values, this involves supressing all leading and trailing spaces.
inputParms <- lapply(strsplit(inputParms," <-"),function(x){
						# There is one condition where we split by c(), the rest by commas
						return( c(x[1],
								if(x[1] != "const_distParameters"){
									tmpReturn <- unlist(strsplit(x[2],',',fixed=TRUE))
									# We now go and try to evaluate 
									tmpReplace <- suppressWarnings(as.numeric(tmpReturn))
									tmpReturn[which(!is.na(tmpReplace))] <- tmpReplace[which(!is.na(tmpReplace))]
									tmpSpaces <- gregexpr("[[:space:]]",tmpReturn)
									for(thisItem in 1:length(tmpSpaces)){
										if(tmpSpaces[[thisItem]][1] != -1){
											# This removes leading spaces
											tmpStart <- if(tmpSpaces[[thisItem]][1] == 1){
															tmpMax <- which(diff(tmpSpaces[[thisItem]]) != 1)
															if(length(tmpMax) == 0){ 
																tmpMax <- 0 
															} else {
																tmpMax <- tmpSpaces[[thisItem]][max(tmpMax)] - 1
															}
															2+tmpMax
														} else {
															1
														}
											# This removes trailing spaces
											tmpEnd <- if(tmpSpaces[[thisItem]][length(tmpSpaces[[thisItem]])] == nchar(tmpReturn[thisItem])){
															tmpMax <- which(diff(tmpSpaces[[thisItem]]) != 1)
															if(length(tmpMax) == 0){ 
																tmpMax <- nchar(tmpReturn[thisItem]) 
															} else {
																tmpMax <- tmpSpaces[[thisItem]][max(tmpMax)+1]
															}
															tmpMax -1
														} else {
															nchar(tmpReturn[thisItem])
														}
											tmpReturn[thisItem] <- substr(tmpReturn[thisItem],tmpStart,tmpEnd)
										}
									}
									tmpReturn
								} else {
									tmpReturn <- unlist(strsplit(gsub("[[:space:]]","",x[2]),'c(',fixed=TRUE))
									tmpReturn <- tmpReturn[which(nchar(tmpReturn) > 0)]
									unlist(lapply(tmpReturn,function(x){ 
											paste('c(',
													if(substr(x,nchar(x),nchar(x)) == ','){
														substr(x,1,nchar(x)-1)
													} else {
														x
													},
													sep="") 
										}))
								} )) 
					})
# I now name the inputParameters given the first item in their list
names(inputParms) <- unlist(lapply(inputParms,function(x){x[1]}))
for(thisList in 1:length(inputParms)){ inputParms[[thisList]] <- inputParms[[thisList]][-1] }
# Now the inputParms of I've allowed the user to input "targetNumber" if they want the "const_focal_popValue"
if(inputParms$deathRate_sizeCap == "targetNumber"){
	inputParms$deathRate_sizeCap <- gsub("[[:space:]]","",gsub("<-","",inputReference[["targetNumber"]]))
}

# This ensures that file pathing objects have terminal slashes
for(thisObject in c("baseDir_local","baseDir_server","workDir")){
	tmpString <- inputParms[[thisObject]]
	if(substr(tmpString,nchar(tmpString)-1,nchar(tmpString)-1) != "/"){
		inputParms[[thisObject]] <- paste(substr(tmpString,1,nchar(tmpString)-1),
										"/",
										substr(tmpString,nchar(tmpString),nchar(tmpString)),
										sep="")
	}
}

# I also grab the parameters output template file so I can get the object sepString
tmpLines <- readLines(fileName_templates["parameter_output"],warn=FALSE)
tmpValue <- tmpLines[which(grepl("sepString <-",tmpLines))[1]]
assign("sepString",trimQuotes(gsub("[[:space:]]","",sub("sepString <-","",tmpValue))),pos=".GlobalEnv")
rm(tmpLines,tmpValue)

# Now we assess if this is a server run and build the appropriate directories, getting the true path requires 
# removing the leading and trailing escape character and quotation marks
tmpDir <- c(inputParms[[if(grepl("TRUE",inputParms[["serverFarm"]])){
								"baseDir_server"
							} else {
								"baseDir_local"
							}]],
				inputParms[["workDir"]])
tmpDir <- paste(trimQuotes(tmpDir),collapse="")
if(!dir.exists(tmpDir)){
	dir.create(tmpDir, recursive = TRUE)
}   
# I also create a template directory 
if(!dir.exists(paste(tmpDir,"templates/",sep=""))){
	dir.create(paste(tmpDir,"templates/",sep=""), recursive = TRUE)
}
# We also copy over the source functions and preProcessing files
file.copy(from = paste(templateDir,allFiles[which(is.element(allFiles,trimQuotes(unlist(inputParms[c("fileName_functions","fileName_preProcessing")]))))],sep=""), 
			to = tmpDir, 
			overwrite = TRUE)

# This is the point where I build a matrix of the parameter combinations 
parameterCombos <- buildCombos(inputParms)
# We write out the parameter combos
write.csv(parameterCombos, file = paste(tmpDir,trimQuotes(inputParms$save_batchBase),"_parameterCombos.csv",sep=""), 
			append = (inputParms$starting_parameterCombination != 1), quote = FALSE, row.names = FALSE)
# Now for each possible parameter combination, and unique statistical replicate, we build directories
tmp_newDirs <- paste(tmpDir,
				name_batchString(funcBase = trimQuotes(inputParms$save_batchBase),
									func_setID = unlist(lapply(1:as.numeric(inputParms$uniqueReplicates),rep,"times"=nrow(parameterCombos))),
									func_jobID = parameterCombos$save_batchJob,
									func_sepString = sepString),
				sep="")
for(thisDir in tmp_newDirs){ 
	if(!dir.exists(thisDir)){
		dir.create(thisDir, recursive = TRUE) 
	}
}
  
# I now copy in templates for all template files, updating the plotting and parameter output as required
for(thisTemplate in names(fileName_templates)){
	if(thisTemplate == "plotting"){
		writePlotting(func_infile = fileName_templates[thisTemplate],func_inParms = inputParms, 
						func_outFile = sub(templateDir,tmpDir,fileName_templates[thisTemplate]))
	} else if (thisTemplate == "parameter_output") {
		writeParameters(func_infile = fileName_templates[thisTemplate], func_inParms = inputParms, 
						func_inCombos = parameterCombos, func_outDir = tmpDir, 
						func_bodyScript = fileName_templates["body"])
	} else if(as.logical(inputParms$serverFarm) && thisTemplate == "serverSubmit" ||
				!as.logical(inputParms$serverFarm) && thisTemplate == "localSubmit") {
		write_subScript(func_subScipt = fileName_templates[thisTemplate], 
						func_outDir = tmpDir, func_inCombos = parameterCombos, 
						func_inParms = inputParms, func_maxJobs= maxJobs)
	}
	# Regardless we'll be writting out the basic template 
	# We copy the file and overwrite any that are in that directory and share a name
	file.copy(fileName_templates[thisTemplate],paste(tmpDir,"templates/",sep=""),overwrite = TRUE)
} 
# I make all .sh files executable
Sys.chmod(list.files(path=tmpDir,pattern='.sh',recursive=TRUE,full.names=TRUE),
		mode="0777")

# I now save an image of this workSpace so that future scripts can load into what has been built here
fileName_farmingImage <- paste(trimQuotes(inputParms$save_batchBase),"_farmerImage.RData",sep="")
save.image(paste(tmpDir,fileName_farmingImage,sep=""))

# I now update the run reset file
tmpLines <- readLines(fileName_templates["run_reset"])
tmpLines <- updateLines(func_inLines = tmpLines,
						func_searchPattern = list("expDir"="experimentDir <-",
													"save_batchString"="saved_batchString <-",
													"farmerImage"="fileName_farmerImage <-",
													"uniqueSets"="startedSets <-",
													"uniqueJobs"="startedJobs <-"),
						func_values = list("expDir"=addQuotes(tmpDir),
											"save_batchString"=inputParms$save_batchBase,
											"farmerImage"=addQuotes(fileName_farmingImage),
											"uniqueSets"=paste("1:",inputParms$uniqueReplicates,sep=""),
											"uniqueJobs"=paste("1:",nrow(parameterCombos),sep="")))
writeLines(tmpLines,con = sub(templateDir,tmpDir,fileName_templates["run_reset"]))


# End
q(save="no")
