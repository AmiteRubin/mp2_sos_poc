import os
import subprocess
import argparse
import numpy as np
import time


def main(input_filename, formula, nk, c_nk, scale_a, ba, ri, exp, code, aux_ri, pyscf_ver, remove_linear, hz_basis,
         pseudu, structure, path_to_struct, path_to_basis_data, cpus_per_task=1, slurm_time='1-00:00:00', mem_mb=20000,
         email=False, job_name=None):
    # we need to adjust towards this data hierarchy to automate the union of data exported
    this_files_path = os.path.realpath(__file__)
    this_dirs_path = os.path.dirname(this_files_path)
    os.chdir(this_dirs_path)
    name_sub_dir = formula
    name_sub_sub_dir = "nk_" + str(int(nk)) + str(int(nk)) + str(int(nk)) + "_sc_" + str("{:.4f}".format(scale_a))
    full_path_sub = this_dirs_path + "/" + name_sub_dir
    full_path_sub_sub = this_dirs_path + "/" + name_sub_dir + "/" + name_sub_sub_dir
    if os.path.isdir(full_path_sub_sub):
        os.chdir(full_path_sub_sub)
        print("True")
    elif os.path.isdir(full_path_sub):
        os.chdir(full_path_sub)
        os.mkdir(name_sub_sub_dir)
        os.chdir(full_path_sub_sub)
        print("True" + formula)
    else:
        os.mkdir(name_sub_dir)
        os.chdir(full_path_sub)
        os.mkdir(name_sub_sub_dir)
        os.chdir(full_path_sub_sub)

    this_files_path_sub = os.path.realpath(__file__)
    this_dirs_path_sub = os.path.dirname(this_files_path_sub)
    print("path", this_dirs_path_sub)
    if remove_linear == 'yes':
        rml = 'rmli'
    else:
        rml = 'nrmli'

    if hz_basis == 'yes':
        if ri == 'ri':
            job_name = job_name + "_" + structure + "_" + pyscf_ver + "_hz_" + ba + "_" + ri + "_" + aux_ri + "_" + rml + "_" + str(
                "{:02d}".format(int(exp * 10))) + "_" + name_sub_dir + "_" + name_sub_sub_dir + "_" + code

        else:
            job_name = job_name + "_" + structure + "_" + pyscf_ver + "_hz_" + ba + "_" + ri + "_" + rml + "_" + str(
                "{:02d}".format(int(exp * 10))) + "_" + name_sub_dir + "_" + name_sub_sub_dir + "_" + code

    else:
        if ri == 'ri':
            job_name = job_name + "_" + structure + "_" + pyscf_ver + "_" + ba + "_" + pseudu + "_" + ri + "_" + aux_ri + "_" + rml + "_" + str(
                "{:02d}".format(int(exp * 10))) + "_" + name_sub_dir + "_" + name_sub_sub_dir + "_" + code
        else:
            job_name = job_name + "_" + structure + "_" + pyscf_ver + "_" + ba + "_" + pseudu + "_" + ri + "_" + rml + "_" + str(
                "{:02d}".format(int(exp * 10))) + "_" + name_sub_dir + "_" + name_sub_sub_dir + "_" + code

    # The input file starts from line 2 (the second line), first line is any comments so
    # don't want the system to read that one!
    # Give the user an option to e-mail themselves when it finishes
    if email:
        email_txt = '''\
#SBATCH --mail-type=END
#SBATCH --mail-user=rubinam4@biu.ac.il
'''
    else:
        email_txt = ''

    # This is a large string of all of the text that will be created. You can add other
    # lines if you wish, or change around what gets formatted in and what is constant
    # between jobs/submission scripts.
    # Note that this assumes you have a file folder named "output" in the folder where
    # this is being submitted from.
    # Before starting the job, the SLURM script cd's into the directory that you SUBMIT
    # the script from, echoes valuable job info into the .out file.
    # The last line of the string actually calls what you want to submit, and in this case
    # runs those lines using python. However, you could just as easily replace python with
    # bash, etc. The command line voodoo is just to get only the portions of the input
    # file that you want, since I have just one set of instructions after another in my
    # python scripts for each of the individual jobs in the SBATCH array.
    slurm_txt = '''\
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={0:d} # maximum of 28, 4 is a good starting max for pyscf
#SBATCH --account=berkelbach
#SBATCH -t {1}
#SBATCH --mem={2:d}M # in MB, ?? GB is the max
#SBATCH -o {3}.out
#SBATCH -e {3}.err
#SBATCH -J {4}
#SBATCH --exclude=umi[3-4]
{5}

module load pyscf/pyscf-{16}
export PYSCF_MAX_MEMORY=$SLURM_MEM_PER_NODE
export PYSCF_TMPDIR=$SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

WORKDIR=$SLURM_SUBMIT_DIR
cd $WORKDIR
echo 'My job ID is' ${{SLURM_JOB_ID}}
echo 'My job name is' ${{SLURM_JOB_NAME}}
#module load openblas
python3 ../../{6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22}
'''.format(
        cpus_per_task,
        slurm_time,
        mem_mb,
        job_name,
        job_name,
        email_txt,
        input_filename, formula,
        nk, c_nk, scale_a, ba, ri, exp, code, aux_ri, pyscf_ver, remove_linear, hz_basis, pseudu, structure,
        path_to_struct, path_to_basis_data
    )

    # Create the SLURM file with the above text
    file_loc = this_dirs_path_sub + '/' + job_name + '.sh'
    # print("file_loc", file_loc)
    slurm_file = open(file_loc, 'w')
    slurm_file.write(slurm_txt)
    # print(slurm_txt)
    slurm_file.close()
    # Submit to the queue and remove the SLURM submission file
    subprocess.call(['sbatch', file_loc])


# if __name__ == '__main__':
# Fill in the description below with what SLRUM submission script this is sending to
# the queue to run.
parser = argparse.ArgumentParser(
    description="script for submitting pyscf calculations for mp2, scs_mp2, sos_mp2, reg_mp2"
)
parser.add_argument(
    "input_filename", type=str,
    help="mp2_scs_reg.py, name of python script running pyscf mp2_scs_reg"
)
parser.add_argument(
    "formula", type=str,
    help="formula of the atom in the unit cell"
)
parser.add_argument(
    "nk", type=float,
    help="number of k-points in isotropic mash [nkx, nky, nkz]"
)
parser.add_argument(
    "c_nk", type=float,
    help="center of k-points in isotropic mash [c-nk, c-nk, c-nk]"
)
parser.add_argument(
    "scale_a", type=float,
    help="rescaling the unitcell volume"
)
parser.add_argument(
    "ba", type=str,
    help="basis set to use"
)
parser.add_argument(
    "ri", type=str,
    help="ri to use ri in aux basis, nri not to use ri in aux basis"
)
parser.add_argument(
    "exp", type=float, default=0.0,
    help="number of CPUs assigned to this job"
)
parser.add_argument(
    "code", type=str,
    help="sos version code to use"
)
parser.add_argument(
    "aux_ri", type=str,
    help="auxiliary basis set-ri"
)
parser.add_argument(
    "pyscf_ver", type=str,
    help="pyscf version to use"
)
parser.add_argument(
    "remove_linear", type=str,
    help="to use remove linear option or not"
)
parser.add_argument(
    "hz_basis", type=str,
    help="to use HZ basis or not"
)
parser.add_argument(
    "pseudu", type=str,
    help="which pseudu basis set to use"
)
parser.add_argument(
    "structure", type=str,
    help="what structure to perform calculation on"
)
parser.add_argument(
    "path_to_struct", type=str,
    help="the path to structures to perform calculations on"
)
parser.add_argument(
    "path_to_basis_data", type=str,
    help="the path to basis sets to perform with calculations"
)
parser.add_argument(
    "-c", "--cpus_per_task", type=int, default=1,
    help="number of CPUs assigned to this job"
)
parser.add_argument(
    "-m", "--memory", type=int, default=20000,
    help="amount of memory (in MB) that this job will be assigned on the cluster, max is ???"
)
parser.add_argument(
    "-t", "--time", type=str, default='1-00:00:00',
    help="amount of time (day-hh:mm:ss) that will be assigned to the job"
)
parser.add_argument(
    "-e", "--email", action='store_true',
    help="whether or not to send an email when the job finishes"
)
parser.add_argument(
    "-j", "--jobname", type=str, default=None,
    help="rerunning only, which jobs to rerun if they failed the first time because of a bug/you fixed the input file part corresponding to this one"
)

# See argparse documentation, this collects the arguments that are given in the
# command line.
args = parser.parse_args()
main(
    args.input_filename, args.formula,
    args.nk, args.c_nk, args.scale_a, args.ba, args.ri, args.exp, args.code, args.aux_ri, args.pyscf_ver,
    args.remove_linear, args.hz_basis, args.pseudu, args.structure, args.path_to_struct, args.path_to_basis_data,
    cpus_per_task=args.cpus_per_task, slurm_time=args.time, mem_mb=args.memory, email=args.email, job_name=args.jobname)



