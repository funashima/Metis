#!/usr/bin/env python3
import subprocess

job_id_start = 4300
job_id_end = 4400

for job_id in range(job_id_start, job_id_end + 1):
    subprocess.call('qdel {}'.format(job_id), shell=True)
