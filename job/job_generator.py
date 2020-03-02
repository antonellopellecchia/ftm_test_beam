#!/usr/bin/python3

import os, sys, stat

sh_fin = open('models/ftm_test_beam.sh')
sh_text_model = sh_fin.read()
sh_fin.close()

job_fin = open('models/ftm_test_beam.job')
job_text_model = job_fin.read()
job_fin.close()

def save_job(job_index):
    dirname = 'jobs/ftm_test_beam_%04d/'%(job_index)
    
    try:
        os.mkdir(dirname)
    except FileExistsError: pass
    print('saving', dirname)

    sh_out = sh_text_model.format('%04d'%(job_index))
    with open(dirname+'job.sh', 'w') as fout:
        fout.write(sh_out)
    os.chmod(dirname+'job.sh', 0o755)

    job_out = job_text_model.format('%04d'%(job_index))
    with open(dirname+'job.job', 'w') as fout:
        fout.write(job_out)

    return dirname

def main(argv):
    with open('job_run.sh', 'w') as run_file:
        #run_file.write('#!/bin/sh/\n\n')
        for job_index in range(1, 11):
            dirname_job = save_job(job_index)
            run_file.write('condor_submit ' + dirname_job + 'job.job -name ettore\n')

if __name__=='__main__': main(sys.argv)
