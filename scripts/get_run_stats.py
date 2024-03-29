#! /usr/bin/env python3

import argparse
import datetime
import os
import re
import subprocess
import time
import yaml

import pandas as pd

try:
    #pylint: disable=locally-disabled,unused-import
    from StringIO import StringIO
except ImportError:
    from io import StringIO


def get_elapsed_secs(timestr):
    gt1day = re.match("^([0-9]+)-.*", timestr)
    if gt1day:
        days = int(gt1day.group(1))
        [hours, minutes, seconds] = [int(x) for x in timestr.split(':')]
        timedelt = datetime.timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)
    else:
        [hours, minutes, seconds] = [int(x) for x in timestr.split(':')]
        timedelt = datetime.timedelta(hours=hours, minutes=minutes, seconds=seconds)
    return(timedelt.seconds)


def get_job_stats(jobid):
    seffrun = subprocess.run(["seff", str(jobid)], capture_output=True, check=True)
    res_raw = seffrun.stdout.decode('utf-8')
    res_noparenth = re.sub("\((.*?)\)", "", res_raw) # parentheses info not informative, maybe problematic, remove
    res_stripped = "\n".join([x.strip() for x in res_noparenth.split("\n")])
    seff_dict = yaml.load(res_noparenth, Loader=yaml.SafeLoader)
    if seff_dict.get('CPU Efficiency'):
        strtest = re.match("([0-9.]+)% of ([0-9:-]+) .*", seff_dict['CPU Efficiency'])
        seff_dict['CPU Efficiency'] = strtest.group(1)
        seff_dict['Core Time'] = strtest.group(2)
    if seff_dict.get('Memory Efficiency'):
        strtest = re.match("([0-9.]+)% of ([0-9:]+.*)", seff_dict['Memory Efficiency'])
        seff_dict['Memory Efficiency'] = strtest.group(1)
        seff_dict['Memory Total'] = strtest.group(2)
    if seff_dict.get('Job Wall-clock time'):
        ## TWS TODO: Follow up after gathering usage data for better handling in the future
        try:
            seff_dict['Elapsed Seconds'] = get_elapsed_secs(seff_dict['Job Wall-clock time'])
        except TypeError as e:
            problem_field = seff_dict['Job Wall-clock time']
            print(f"Problem field: {problem_field}")
            print(f"Problem field type: {type(problem_field)}")
            raise e

    return(seff_dict)

def get_job_info(jobid):
    sacctrun = subprocess.run(["sacct", "-P", "-o", "JobID,JobName,Partition,User,Account,State,ExitCode,ReqCPUS,ReqMem,ReqTres", "-j", str(jobid)], capture_output=True, check=True)
    res_raw = sacctrun.stdout.decode('utf-8')
    res_strIO = StringIO(res_raw) # Can load a file-like object into pandas
    res = pd.read_csv(res_strIO, sep = '|')
    keep = res[res['JobID'] == jobid].set_index('JobID')
    # Pull "billing" metric from ReqTRES results
    keep["TRES Billing"] = keep["ReqTRES"].apply(lambda x: int(re.sub("billing=(.*?),.*", "\\1", x)))
    # Convert to dict
    res_dict = keep.to_dict(orient='records')
    res_dict = res_dict[0] # It's a dict inside of a 1-length list
    return(res_dict)

def main(args):
    jobs = {}
    with open(args.infile) as infile:
        for line in infile:
            match = re.match('Submitted job (\d+) with external jobid \'(\d+)\'', line)
            if match:
                clust_jobid = match.group(2)
                stats = get_job_stats(clust_jobid)
                info = get_job_info(clust_jobid)
                jobs[clust_jobid] = stats
                jobs[clust_jobid].update(info)
                # Calculate cost
                cost = {'Cost': jobs[clust_jobid]['TRES Billing'] * jobs[clust_jobid]['Elapsed Seconds'] / 60}
                jobs[clust_jobid].update(cost)

    # Create outfile name
    if not args.outfile:
        run_timestamp = re.sub("\.snakemake\.log", "", args.infile)
        args.outfile = "{}.info.yaml".format(run_timestamp)

    with open(args.outfile, "w") as outfile:
        yaml.dump(jobs, outfile, default_flow_style=False, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='python piperun_info.py', description = 'Gets cluster job ID from snakemake log, calls seff/sacct to get job stats and information for each job in a run.')
    parser.add_argument('-i', '--infile', required=True, help='Input file to process')
    parser.add_argument('-o', '--outfile', help='Filename for output. Defaults to {timestamp}.info.yaml, where timestamp is taken from input file {timestamp}.snakemake.log. The output will be placed alongside the input snakemake log.')

    args = parser.parse_args()

    main(args)
