##################
# Update history #
##################

2016-06-02:
	Created	Cluster project.

2016-06-08:
	Removed	trimming process as BWA & GATK tools can handle quality based trimming and adapters.

2016-06-13:
	Created	base reference file to reduce data replication across multiple scripts.

2016-06-28:
	Added	contig definition to for automation.

2016-07-01:
	Added	FASTQ scan definitions.

2016-07-05:
	Added	SLURM based memory calculation to java arguments.

2016-07-07:
	Added	check for previously finished jobs.

2016-07-11:
	Added	global module version definitions.

2016-07-12:
	Added	duplicate metrics storage for sample.

2016-07-14:
	Added	parallelization option to HaplotypeCaller function as this has been resolved in GATK and it will reduce runtime.

2016-07-15:
	Added	SLURM based parallelization calculation to GATK arguments.
	Changed	HaplotypeCaller mem-per-cpu from 16G to 32G to test if extra memory decreases run-time

2016-07-18:
	Changed	HaplotypeCaller cpu-per-task from 1 to 2 to test if -nct parallelization has been resolved in this version of GATK.

2016-07-19:
	Created	CatPrintReads & Index scripts to replace MergeReads.
			Picard MergeSamFiles takes 6 to 10 hours to complete.
			SAMTools Cat takes 30 minutes.
	Added	Global temp directory definition to minimise network thrashing.

2016-07-20:
	Added	scriptFailed function to collect basic node data when a job fails for any reason.
	Fixed	Coverage & Gender Determination not collecting GL* PrintRead dependencies.
	Changed	HaplotypeCaller cpu-per-task & mem-per-cpu from 2x32 to 4x8 to boost parallelization but reduce overall memory requirement as this extra memory no longer affects run-time.
	Changed	SortSam mem-per-cpu from 32G to 16G as this has no effect on runtime.

2016-07-21:
	Changed	CatVarInputs to generate independent of main cycle so X&Y don't need to be resorted by CatVar.
	Changed	Global temp directory to be SLURM provided one so minimize chances of job fail due to temp directory failure.
	Changed	CatPrintReads walltime from 3 hours to 1. 6x to 2x 30m known run-time.
	Changed	CatPrintReadsIndex walltime from 1 hour to 1.5. 1.25x to 2x 45m known run-time.
	Changed	HaplotypeCaller walltime from 6 hours to 3 as parallelization with -nct has proven effective and cpus-per-task & mem-per-cpu from 4x8192 to 8x4096.
	Changed	CatVariants walltime to 3 hours as longest runtime is ~2 hours
	Changed	SortSam walltime from 3 hours to 1.5 as longest runtime is 30 minutes
	Changed	Mark Duplicates memory allocation from 32G to 16G as no different in runtime and walltime from 6 hours to 2 as longest runtime is just over 1 hour.
	Added	scriptFailed function to all non-0 exit points.
	Added	scontrol show job $SLURM_JOBID to scriptFailed function.
	Added	Reads file size to job name to size-vs-speed filtering.
	Created	Automatic sample progression. (ASP)

2016-07-22:
	Fixed	X sub contigs not being linked correctly in their folders .../X/X:1-2699520.g.vxf.gz, etc.
	Fixed	MT contig being included in HaplotypeCaller after removing input collection from primary loop.

2016-07-25
	Added	Check for job submission failures. Pipeline will exit 1 in that event.
	Added	Local spooler lauching sequencial jobs without input.

2016-7-27
	Fixed	ASP seeing a failed grep as a failed ssh connection. grep failure is ExitCode 1, ssh network failure exitcode is > 1

2016-7-28
	Added	Contig Count to baserefs.sh

2016-07-29
	Changed	ContigSplit function to be array job. Now submits 1 job array per block instead of 84. Array is dynamically created based on .done file existence.
	Merged	CatReads and CatReadsIndex jobs into 1. Increase walltime to accomodate both jobs. ~30m cat, 45m index.
	Changed	Log output header to 2 character combos. RS: Read Split, PA: Primary alignment, SS: SortSAM, etc...
	Changed	Log output to single line per post-merge contig.
	Changed	changelot.txt to README.md so changelog is visible.
	
#########
# To do #
#########

Investigate job arrays to better manage bulk job sumission. NeSI queue length is 5000 jobs. Each array job can contain 1000 sub jobs?

Find out why MarkDuplicates sometimes fails to find input. The dependencies are transmitted but the output doesn't exist when MD starts.

Continue to work out minimum requirement to obtain 1 hour max runtimme per segment within high partition.

Build multi-sample pipeline with identity comparison and merge function.

Add fingerprinting function post bam collection?

Build list of times for each job type and contig to set per-contig wall-times.

Add ability to re-try a job if it fails because of a node issues. Could trigger a chain of dependant jobs to take it as far as it can outwith the overall process. This will minimize delay on restarting the job later on.

Migrate input file collection to slurm script. This will allow failed jobs to possibly be run in time for the catvar or catreads functions to collect them. Still wont be 100% effective as job execution depends on cluster availability. Should work well for smaller contigs.

Re-write everything. Check-Blocks is a bloated mess.
