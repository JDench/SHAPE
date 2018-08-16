# SHAPE
Simulated Haploid Asexual Population Evolution

Thank you for choosing to use my tool SHAPE (Simulated Haploid Asexual Population Evolution).

This set of scripts was designed to help study evolutionary dynamics at the level of both 
population demographics and genetics.  The scripts allow users to define a range of parameters 
which are then built into a fully factorial experiment that can be run.

You don't need to be able to code to work with SHAPE, you need only know how to open an R script 
and edit text held within.  Just know that wherever you've copied the set of scripts for SHAPE is 
considered the template directory and you should not add new files into that folder.  When
you run SHAPE it will build a new directory and copy all the necessary files for you to run
your experiment.

Two important definitions when reading through SHAPE's guidelines and comments are the terms: 
filepath - which means the full named path on your system that leads to a file or directory. 
string literal - which means a set of words, without special characters, held between double quote 
		e.g.  "someWords"
	NOTE: I recommend you avoid using underscores (ie: "_") in anything except your filepaths

To define parameters simply open the "SHAPE_farmerParms.v.#.r" file and change parameters.  Read
the comments above each parameter assignment to better understand what and how it might be used.
Once you've defined your experiment's parameters open the "SHAPE_farmer.v.#.r" and change at least
the lines with:		templateDir <-  
		and			maxJobs <-
		and			appLocation <-
which are near the top.  The first line is the filepath to the directory where you've copied
the template scripts for SHAPE.  The second controls how many concurrent jobs you accept to have 
running.  Note that if you're running this on your local computer you will see best performance by 
setting this number not higher than the number of true (not hyperthreaded) processors for your 
system minus 1, and likely subtracting that number by one wouldn't hurt either.  The third is the 
full filepath and executable filename for your installed version of R.

The R software environment is required to run the tool and once you've updated the two "farmer" 
script files you can simply run "SHAPE_farmer.v.#.r" which should build your experiment in the 
directory that you defined in the parameters file.  Go to that directory and you can run the: 
"submit_<your-defined-experiment-named>_#.sh" 
files to get your experiment going unless you're not using a bash shell.  All the automated 
farming of experimental control assumes a bash shell so without this you have to manually
go into each subfolder within the root experiment directory and run the scripts named:
"SHAPE_body_<your-defined-experiment-named>_.r"
The farmer parameters allow users to define if a job is to be run on a local machine or a 
remote server with a SLURM job system, such as that of the Center for Advanced Computing at 
Queens University.

SHAPE is not going to finish quickly (ie: within a day) unless you've designed a very small experiment.  
I cannot give a good sense of the time that will be required except to note that the length of your genome 
and mutation supply rate are both going to be the biggest factors that slow runs.  As an example I've 
run an experiment with SHAPE that had 162 parameter combinations, 1000 parameter replicates, 
10 statistical replicates, a genome length of 1000 and a mutation supply rate going as high as ~500
mutants per generation (with 2000 generations per replicate) and it took about 1 month to complete
while running as many as 760 concurrent jobs.  Also, at the end more than 1 Tb of data was produced.
However, while stress testing SHAPE I've run sets with 32 parameter combinations, a genome length of 100,
and mutation rates leading to ~50 mutations per generation and have completed these overnight on a server.

Once data generation is complete SHAPE will have created numerous raw output files detailing the complete
through time demograhpics of the experiment and the fitness of every genotype.  The raw demographics output
is costly in terms of disk space and so you can set that this be removed with the parameter:
results_removeSteps <- TRUE
Regardless, SHAPE will also have output some processed summary RData files relating to each of the demographic
outputs.  These can be found in the "processed_runData_from_XXX.RData" file which contain two R lists "runDemographics"
and "info_estLines".  The first gives summaries of the per generation demographics only lacking specifics 
concerning which genotypes had which births, deaths and mutations (these would be found in the raw output 
databases).  The second details information concerning genotypes which became established ( controlled with
the "size_ofEstablishment" parameter) and provides information concerning which genotypes are part of 
their evolutionary trajectory.  

Once all the jobs of your experiment are complete you should run the "SHAPE_plotting.v.#.r" script found in 
your experiment's directory.  That will process the results of your experiment, generate several plots I find 
usefull, but also produce several RData files that contain a wealthy summary of information regarding your 
experiment.  Using those files, and working with the "SHAPE_plotting.v.#.r" script, you should be able
to run most any analyis you'd like.  I have not built scripts to extract and visualise every possible analytic
conceivable but the files provided should give you a sense of how you could get there yourself provided you 
know how to do a bit of R scripting.  If there are particular plots or analytics you'd like, or think should be
added as standard output for this script please let me know and I'll try to work with you on it.

Regardless, I hope you'll find SHAPE usefull for supporting theoretical and empirical study of 
evolution in an asexual system.  I'll note that I've done my best to make SHAPE robust to errors
but it's not perfect.  If you try to input parameters that will break the program I have no doubt you can.
Yet, provided you're trying to work with the tool it should run smoothly.  If you find your runs are failing
I have included the "SHAPE_farmingReset.v.#.r" script, which will be copied to your experiment directory,
that you can run to check what jobs have completed sucessfully and to automaticall re-submit ones that 
have failed.  
	NOTE: At present there is a rare instance that causes SHAPE to crash, I have not had success (due to limited
			time) to successfully troubleshoot this, but seems to be related to parameter combinations which result
			in few mutants that lead from different parental genotypes to the same offspring genotype.	

Again, thanks for choosing to use SHAPE, please feel free to contact me if you have any questions or concerns:

Sincerely,

Jonathan Dench
Jdenc017@gmail.com


DEPENDENT R Libraries - You need these installed for SHAPE to run:
abind
colorRamps
compiler
DBI
evd
foreach
lattice
plot3D
RColorBrewer
RSQLite
sn
VGAM

AND either
doSNOW   - For Windows based machines
doMC	 - For all others

##############################################################################
############################## CITATION ######################################
##############################################################################
If you plan to cite SHAPE the related manuscript is pending peer review submission but can be found on BioRxiv:
Dench, Jonathan (2018) The SHAPE of logistic growth shows that timing does matter, bioRxiv, doi:10.1101/392944

BibTex Citation Script:
@article {Dench:2018aa,
	author = {Dench, Jonathan},
	title = {The SHAPE of logistic growth shows that timing does matter},
	year = {2018},
	doi = {10.1101/392944},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Experimental evolution is a powerful tool for studying evolutionary questions and the use of in silico systems is growing as they offer complete parameter control and more fine grained tracking of dynamics. However, existing software are limited by the models implemented, output obtainable, or their lack of interpretability within biological systems. Here I introduce SHAPE (Simulated Haploid Asexual Population Evolution) a forward time simulation tool closely replicating microbial experimental evolution that offers the flexibility to implement alternative models and study detailed output. I first validate SHAPE by replicating the evolutionary outcome expected by several theoretical works. I then extend theory by contrasting how serial passaging affects the evolutionary outcome for microbial communities undergoing exponential and logistic growth.},
	URL = {https://www.biorxiv.org/content/early/2018/08/16/392944},
	eprint = {https://www.biorxiv.org/content/early/2018/08/16/392944.full.pdf},
	journal = {bioRxiv}
}
