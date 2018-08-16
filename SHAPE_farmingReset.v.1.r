# This is for use with remote server's that are based on the SLURM submissin system, such as on the CAC Queens University servers.
# It uses the farmerImage file created in the experiment's working directory.
#### NOTE: This script has no current application on a system which does not have a SLURM job system.
rm(list=ls())
##############################################################################################################
####################################### CHANGE ELEMENTS WITHIN HERE ##########################################
##############################################################################################################
# Enter your unique ID for the server, such that it can be querried for your jobs.... (as string)
userID <- "hpc3058"
# From the initiation of the experiment, with the farmer, these are the number of ALL
# unique statistical replicate sets, and parameter combinations
### CAUTION: These will initially be set to assume all jobs have been initiated, only YOU
###                   can meaningfully reset these to reflect what has actually been started.
startedSets <- 1:3
startedJobs <- 1:32

##############################################################################################################
##############################################################################################################
##############################################################################################################

saved_batchString <- "basicSHAPE"

# This is your filepath location for finding the experiment and assocaited files.
experimentDir <- "/global/home/hpc3058/jonathan/SimulationRuns/multiMutants/basicSHAPE/"
fileName_farmerImage <- "basicSHAPE_farmerImage.RData"
if(!file.exists(paste(experimentDir,fileName_farmerImage,sep=""))){
	stop("Could not find the reference information required for performing the farming reset.  Pelase review.")
	Sys.sleep(5)
	q(sve="no")
} else {
	load(paste(experimentDir,fileName_farmerImage,sep=""))
}
# this indicates where we're storing information with respect to the jobs that have been previously IDed as complete
completedFile <- paste(experimentDir,"farmedOut_completedJobs.RData",sep="")
completedJobs <- NULL
if(file.exists(completedFile)){
	load(completedFile)
}

tmp_expectedJobs <- name_batchString(funcBase = saved_batchString, 
									func_setID = unlist(lapply(startedSets,rep,"times"=length(startedJobs))),
									func_jobID = rep(startedJobs,length(startedSets)))
# This is the regular expression for our batch string
regexpr_batchString <- name_batchString(funcBase = saved_batchString, 
									func_setID = '([[:digit:]])+',
									func_jobID = '([[:digit:]])+')
									
# We querry the server to know which jobs are running
tmp_fileName <- paste(experimentDir,"tmp_querryInfo_",userID,".txt",sep="")
system(paste("squeue -u ",userID," > ",tmp_fileName,sep=""))
queeryTable <- readLines(tmp_fileName)
if(file.exists(tmp_fileName)){ file.remove(tmp_fileName) }
# These are the jobs that should exist.
existJobs <- unlist(lapply(queeryTable,function(x){ 
					tmpReturn <- regexpr(regexpr_batchString,x)
					# We querry for our batch string expression was found, meaning a non -1 answer
					if(tmpReturn == -1){
						return( NULL )
					} else {
						# If something exists we substr out our information
						return( substr(x,tmpReturn, tmpReturn+attr(tmpReturn,"match.length")-1) )
					} }))
# We know which jobs are missing based on those that we expect vs those that are running and complete.
missJobs <- setdiff(tmp_expectedJobs,existJobs)
missJobs <- setdiff(missJobs, completedJobs)

# This was the maximum number of reps requested
maxReps <- as.numeric(inputParms$maxReplicates)

# Now that we know where the jobs are run, then we go an list all the directories that exist
allDirs <- list.dirs(path=experimentDir, recursive=FALSE, full.names=FALSE)
# This is the string which defines processed information of a run
completionString <- "processed_runData_from_"
removeSuffix <- '.RData'
# This will just help me track the types of restarts I've had to do.
restartConditions <- list("Continue"=list("Name"="",
										"Number"=0),
							"Completed"=list("Name"="",
										"Number"=0))
new_completedJobs <- NULL

# We run through the missing jobs
for(thisJob in missJobs){
	# I separate the set and job information
	jobInfo <- name_batchString(funcBase = thisJob, 
									func_setID = TRUE,
									func_jobID = TRUE,
									funcSplit = TRUE)
	# We go an look for a directory that conforms to exactly this set and job information
	jobDir <- paste(allDirs[which(allDirs == thisJob)],"/",sep="")
	# Make certain we've found exactly one job directory
	if(length(jobDir) != 1){
		stop(paste("There was a problem defining a unique job direcotry for ",thisJob,sep=""))
		next
	}
	
	# Look for output data signifying a job is completed.
	repOuts <- list.files(path = paste(experimentDir,jobDir,"/",sep=""), pattern= completionString)

	# Also try and find the job submission script
	tmpScript <- list.files(path = paste(experimentDir,jobDir,sep=""), pattern = name_subScript(thisJob))
	# Make certain we've found exactly one script
	if(length(tmpScript) != 1){
		stop(paste("There was a problem defining a unique submission script for ",thisJob,sep=""))
		break
	}
	
	if(length(repOuts) == 0){
		# This means the job never completed one run.  We delete any files that may exist for Steps
		removeFiles <- 	c(list.files(path=paste(experimentDir,jobDir,"/",sep=""), pattern = "_Neighbourhood_"),
							list.files(path=paste(experimentDir,jobDir,"/",sep=""), pattern = "_Steps_"))
		if(length(removeFiles) > 0){
			file.remove(paste(experimentDir,jobDir,"/", removeFiles,sep=""))
		}
		# Reset the job.
		system(paste("sbatch ", experimentDir, jobDir, tmpScript,sep=""))
	} else {
		# Extract the replicate values of jobs
		tmpOuts <- as.numeric(name_batchString(funcBase = gsub(removeSuffix,"",gsub(completionString,"",repOuts),fixed=TRUE), 
												func_setID = TRUE,
												func_jobID = TRUE,
												func_repID = TRUE,
												funcSplit = TRUE)["repID",])
		
		# If the maximum output exists the job is complete, otherwise we update the job submission script and resubmit
		# after removing any Neighbourhood and Steps files with stems greater than the maxOut
		if(max(tmpOuts) >= maxReps){
			new_completedJobs <- c(new_completedJobs, thisJob)
			restartConditions[["Completed"]][["Number"]] <- restartConditions[["Completed"]][["Number"]] + 1
			restartConditions[["Completed"]][["Name"]] <- c(restartConditions[["Completed"]][["Name"]],thisJob)
		} else {
			# Here is where we query and Step and Neighbourhood files and remove them, but either all or just the final
			# depending on the inputParameters about intermediate files
			removeFiles <- 	c(list.files(path=paste(experimentDir,jobDir,"/",sep=""), pattern = "_Neighbourhood_"),
								list.files(path=paste(experimentDir,jobDir,"/",sep=""), pattern = "_Steps_"))
			if(!as.logical(inputParms$results_removeSteps)){
				removeFiles <- removeFiles[which(grepl(paste(max(tmpOuts)+1,'.sqlite',sep=""),removeFiles))]
			}
			if(length(removeFiles) > 0){
				file.remove(paste(experimentDir,jobDir,"/", removeFiles,sep=""))
			}
			# We load the script, update to the next replicate, add a line to copy the Landscape and parameters files, then submit.
			tmpScript_fileName <- paste(experimentDir,jobDir,tmpScript,sep="")
			# This would be to advanced the currReplicate, but since the current one failed, we can just re-submit...
			tmpLines <- readLines(tmpScript_fileName)
			# This is the line which includes the call argument
			tmp_callLine <- which(grepl("currReplicate=",tmpLines))
			tmpLines[tmp_callLine] <- sub('currReplicate=[[:digit:]]+',paste('currReplicate=',max(tmpOuts)+1,sep=""),tmpLines[tmp_callLine])
			
			# We now write out this submission script into the experimentDir
			writeLines(tmpLines, con= tmpScript_fileName)
			Sys.chmod(tmpScript_fileName, mode="0777")
			system(paste("sbatch ", tmpScript_fileName,sep="")) 
				
			restartConditions[["Continue"]][["Number"]] <- restartConditions[["Continue"]][["Number"]] + 1
			restartConditions[["Continue"]][["Name"]] <- c(restartConditions[["Continue"]][["Name"]],thisJob)
		}
	}
}
# We report on the conditions of restarts and then exit.
print(restartConditions)

completedJobs <- unique(c(completedJobs, new_completedJobs))
save(completedJobs,file= completedFile)


if(all(is.element(tmp_expectedJobs,completedJobs))){
	print(paste("Completed all jobs in ",saved_batchString,sep=""))
}


q(save="no")

