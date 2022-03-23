#!/usr/bin/env python
import sys
import os
import time
import subprocess

def submit(command, runtime, cores, ram, directory='', modules='', group='vaccine',
    jobscript='jobscript', output='/dev/null', error='/dev/null'):
    """
    Function to submit a job to the Queueing System - with jobscript file
    Parameters are:
    command:   The command/program you want executed together with any parameters.
               Must use full path unless the directory is given and program is there.
    directory: Working directory - where should your program run, place of your data.
               If not specified, uses current directory.
    modules:   String of space separated modules needed for the run.
    runtime:   Time in minutes set aside for execution of the job.
    cores:     How many cores are used for the job.
    ram:       How much memory in GB is used for the job.
    group:     Accounting - which group pays for the compute.
    jobscript: Standard name for the jobscript that needs to be made.
               You should number your jobscripts if you submit more than one.
    output:    Output file of your job.
    error:     Error file of your job.
    """
    runtime = int(runtime)
    cores = int(cores)
    ram = int(ram)
    if cores > 10:
        print("Can't use more than 10 cores on a node")
        sys.exit(1)
    if ram > 120:
        print("Can't use more than 120 GB on a node")
        sys.exit(1)
    if runtime < 1:
        print("Must allocate at least 1 minute runtime")
        sys.exit(1)
    minutes = runtime % 60
    hours = int(runtime/60)
    walltime = "{:d}:{:02d}:00".format(hours, minutes)
    if directory == '':
        directory = os.getcwd()
    # Making a jobscript
    script = '#!/bin/sh\n'
    script += '#PBS -A ' + group + ' -W group_list=' + group + '\n'
    script += '#PBS -e ' + error + ' -o ' + output + '\n'
    script += '#PBS -d ' + directory + '\n'
    script += '#PBS -l nodes=1:ppn=' + str(cores) + ':thinnode,mem=' + str(ram) + 'GB' + '\n'
    script += '#PBS -l walltime=' + walltime + '\n'
    if modules != '':
        script += 'module load ' + modules + '\n'
    script += command + '\n'
    if not jobscript.startswith('/'):
        jobscript = directory + '/' + jobscript
    with open(jobscript, 'wt') as jobfile:
        jobfile.write(script)
    # The submit
    job = subprocess.run(['qsub', jobscript],stdout=subprocess.PIPE, universal_newlines=True)
    jobid = job.stdout.split('.')[0]
    return jobid


# wdir = "/home/mathias/Bachelor/Bachelor/code/specialkursus/tmp"    # local
wdir = "/home/projects/vaccine/people/matbor"                    #computerome

rand_pep_files = ["smaller_rand_20mers_aa","smaller_rand_20mers_ab","smaller_rand_20mers_ac","smaller_rand_20mers_ad","smaller_rand_20mers_ae"]
netMHCIIpan_path = "/home/projects/vaccine/people/morni/netMHCIIpan-4.1/netMHCIIpan"
allele_file = open(sys.argv[1], "r")

if not os.path.exists(wdir + "/jobscripts"):
    os.mkdir(wdir + "/jobscripts")
if not os.path.exists(wdir + "/results"):
    os.mkdir(wdir + "/results")

for allele in allele_file:
    for i, pep_file in enumerate(rand_pep_files):
        allele = allele.strip()
        pep_file_path = wdir + "/data/" + pep_file
        xls_name = wdir + "/results/netMHCIIpan_{}_{}.xls".format(allele, i)

        submit(netMHCIIpan_path + " -a " + allele + " -f " + pep_file_path + " -xls " + "-xlsfile " + xls_name,
               runtime = 480, cores=1, ram=2, directory=wdir + "/jobscripts", group='vaccine',
               jobscript='job_{}_{}'.format(allele, i), output='/dev/null', error='/dev/null')

        time.sleep(1)

allele_file.close()
