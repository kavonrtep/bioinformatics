# How to run genome assembly on Metacentrum
Metacentrum Documentation:
- pbs assembler : https://metavo.metacentrum.cz/pbsmon2/qsub_pbspro
- Metacentrum personal view: https://metavo.metacentrum.cz/pbsmon2/person
- Running simple job: https://docs.metacentrum.cz/computing/run-basic-job/
- Frontends - https://docs.metacentrum.cz/computing/infrastructure/frontends/
- Storages - https://docs.metacentrum.cz/computing/infrastructure/storages/
- Available software


# Running Genome Assembly on MetaCentrum (Velvet Example) 
When working with large-scale computations, such as genome assembly, your personal computer often lacks the necessary computational power, memory, or time to complete the task efficiently. MetaCentrum provides a cluster of powerful computing nodes managed by a job scheduler called PBS (Portable Batch System). Using PBS:

You submit jobs with specific resource requirements (CPUs, memory, walltime).
PBS distributes these jobs among the available compute nodes.


---

## Understanding the PBS Parameters 

The first lines of your PBS script define how the job is scheduled:
 
- `#PBS -N velvet_assembly`
**Job Name** : Sets the name of your job, making it easier to identify in the queue.
 
- `#PBS -l select=1:ncpus=4:mem=16gb:scratch_local=50gb`
**Resources Selection** : Requests 1 compute node (`select=1`), 4 CPU cores (`ncpus=4`), 16 GB of RAM (`mem=16gb`), and 50 GB of local scratch disk space (`scratch_local=50gb`).
 
- `#PBS -l walltime=12:00:00`
**Time Limit** : The maximum wall time for your job. If it exceeds 12 hours, the job will be terminated.
---


## Log Files (Standard Output and Standard Error) 
By default, when you submit a PBS job, two log files are created in the directory where you ran `qsub`: 
- `velvet_assembly.o<JOBID>`: This file contains the **standard output (stdout)**  of your job. It logs all the normal output messages produced by your scripts and commands.
 
- `velvet_assembly.e<JOBID>`: This file contains the **standard error (stderr)**  of your job. Any error messages or warnings are written here.

These logs are useful for:

- Debugging errors if the job fails.

- Reviewing the progress and performance of your job.

- Ensuring that the workflow ran as expected.


---
## Why Use a Scratch Directory (`$SCRATCHDIR`)? 
When your job runs on MetaCentrum, it executes on a compute node. Each compute node provides a local **scratch directory**  (`$SCRATCHDIR`) which is a temporary storage space unique to that node. Using `$SCRATCHDIR` offers several benefits: 
- **Faster I/O Performance:**  Reading and writing data directly on the compute node’s local storage is typically faster than accessing data over networked filesystems.
 
- **Reduced Network Load:**  Running I/O-heavy operations directly in `$SCRATCHDIR` helps prevent slowdowns that could occur if multiple jobs were simultaneously reading and writing from a shared network directory.
 
- **Efficiency and Stability:**  Your job is less likely to experience I/O bottlenecks, and other users’ workflows are not impacted by your job’s data transfers.
At the end of your job, you copy the final results back to your permanent home directory. Once the job finishes, `$SCRATCHDIR` is automatically cleaned.


## Step-by-Step Instructions 

### Step 1: Obtain Data Locally 
Use `wget` on your local machine to download the input FASTQ files:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_2.fastq.gz
```
After these commands, you have `SRR022852_1.fastq.gz` and `SRR022852_2.fastq.gz` on your local machine.
### Step 2: Create a Directory on the Remote Server Using FileZilla 
 
1. Start FileZilla and connect:
 
  - **Host** : `nympha.meta.zcu.cz`
 
  - **Protocol** : SFTP - SSH File Transfer Protocol
 
  - **Username** : Your MetaCentrum username

After logging in, navigate to your home directory on the remote server.
 
2. Right-click and select **"Create directory"** . Name the directory `ngs_assembly`.
Your directory path will be something like:


```arduino
/auto/plzen1/home/username/ngs_assembly
```
Replace `username` with your actual MetaCentrum username.

### Step 3: Transfer Data to the Remote Server 

Use FileZilla to transfer the two FASTQ files from your local machine to:


```bash
/auto/plzen1/home/username/ngs_assembly
```

### Step 4: Prepare the PBS Job Submission Script 
 
1. **Create the script locally**  using a text editor (e.g., `gedit`). Name it `run_assembly.sh`:

2. **Important** before saving script change `username` string to your real `username` 

```bash
#!/bin/bash
#PBS -N velvet_assembly
#PBS -l select=1:ncpus=2:mem=4gb:scratch_local=50gb
#PBS -l walltime=1:00:00

# Set data directory in your home directory
DATADIR=/auto/plzen1/home/username/ngs_assembly

# Load the velvet module
module add velvet

# Copy input data from home directory to scratch directory
cp $DATADIR/SRR022852_1.fastq.gz $SCRATCHDIR
cp $DATADIR/SRR022852_2.fastq.gz $SCRATCHDIR

# Move into the scratch directory
cd $SCRATCHDIR

# Run velvet: using local scratch for input and output
velveth run_25_paired 25 -fastq.gz -shortPaired -separate SRR022852_1.fastq.gz SRR022852_2.fastq.gz
velvetg run_25_paired -ins_length 350

# Copy results back to home directory
cp -r run_25_paired $DATADIR/
```
Replace `username` in `DATADIR` with your actual username.
 
2. **Transfer the script to the remote server** :
Use FileZilla to upload `run_assembly.sh` to `/auto/plzen1/home/username/ngs_assembly`.

### Step 5: Submit the Job 
 
1. **Log in via SSH** :

```bash
ssh username@nympha.metacentrum.cz
```
 
2. **Navigate to your working directory** :

```bash
# see current location
pwd
# go to directory with data and pbs script
cd /auto/plzen1/home/username/ngs_assembly
# or simply 
cd ngs_assembly # if you are already in /auto/plzen1/home/username
```
NOTE: Some network directories on metacentrum have alternative names, so for example `/auto/plzen1/home/username/` point to the same location as `/storage/plzen1/home/username/` 

3. **Submit the job** :

```bash
qsub run_assembly.sh
```
After execution you should see message like:
```bash
6724234234.pbs-m1.metacentrum.cz
```
This is the job ID. 
The job will be placed in the queue, and once it starts running, it will use the requested resources and `$SCRATCHDIR`.
 
4. **Check job status** :
On frontend shell you can use command:
```bash
qstat -u username
```
Status can be also monitored on https://metavo.metacentrum.cz/pbsmon2/person (this require login with EduID)

Once the job finishes, you’ll find a directory `run_25_paired` with `contigs.fa` and other assembly results in your `ngs_assembly` directory.
Working directory will also contain files:
```bash
velvet_assembly.o<JOBID>
velvet_assembly.e<JOBID>
```
This files contain standard output and standard error messages produced by script and are usefull when debugin pbs script.

### Step 6: Retrieve Results 
After the job is complete, use FileZilla again to connect to `nympha.meta.zcu.cz` and navigate to:

```arduino
/auto/plzen1/home/username/ngs_assembly/run_25_paired
```
Download `contigs.fa` (and any other files, for example complete directory with results) to your local machine.

---

**Summary:**  
1. **Local Setup** : Download FASTQ files with `wget`.
 
2. **Remote Setup** : Create a directory on the remote server with FileZilla.
 
3. **Data Transfer** : Upload input files to the remote directory using FileZilla.
 
4. **Job Script** : Create a PBS job script locally, then upload it to the remote server.
 
5. **Submit and Run** : Use `qsub` to submit the job. It uses `$SCRATCHDIR` for efficient I/O.
 
6. **Results and Logs** : Once finished, retrieve assembly results and check log files (`.o<JOBID>` and `.e<JOBID>`) to ensure successful completion.
 
7. **Download Output** : Use FileZilla to bring the results back to your local machine.

By following these steps, you efficiently run a genome assembly on MetaCentrum, taking advantage of scratch space for faster performance and properly managing logs and resources.
