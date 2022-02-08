#!/bin/bash

function show_help {
    echo "Usage: ./sbatch_array.sh -f <cmdfile>"
    echo "  Submit commands to a slurm/sbatch system using an array job, where <cmdfile> is a file containing a list of commands to execute, one per line."
    echo "OPTIONS:"
    echo "  -f <val>: The name of the command file to submit. REQUIRED."
    echo "  -n <val>: Set the number of processes per command in the sbatch script. (Default: 1)"
    echo "  -w <val>: Set the walltime requested. (Default: 120, measured in hours)"
    echo "  -P <val>: Set the number of processes per node per command in the sbatch script. (Default: 1)"
    echo "  -p      : Pretend: do not actually run the commands, but generate the sbatch scripts for examination."
    echo "  -v      : verbose mode."
    echo "  -h      : This help."
}

cmdfile=""
walltime=120  #hours
nodes=1
ppnode=1
pretend=0
while getopts "h?vpP:n:w:f:" opt; do
    case "$opt" in
    f) cmdfile=$OPTARG
	;;
    w) walltime=$OPTARG
        ;;
    n) nodes=$OPTARG
        ;;
    P) ppnode=$OPTARG
        ;;
    p)  pretend=1
        ;;
    h|\?)
        show_help
	exit 0
        ;;
    esac
done

if [ $# -le 1 ]; then
    show_help
    exit 1
fi

slurmfile="$cmdfile.slm"
if [ ! -f $cmdfile ]; then
    echo "Command file $cmdfile does not exist!"
    exit 1
fi
ncmds=`cat $cmdfile | wc -l`
log=`mktemp -d $HOME/log/${cmdfile}.XXXXXX`
echo "Making sbatch script $slurmfile from $cmdfile with $ncmds commands"
echo "#!/bin/bash -login" > $slurmfile
echo "#SBATCH --time=$walltime:00:00"  >> $slurmfile
echo "#SBATCH --array=1-$ncmds" >> $slurmfile
echo "#SBATCH --nodes=$nodes" >> $slurmfile
echo "#SBATCH --cpus-per-task=$ppnode" >> $slurmfile
echo "#SBATCH --output=$log/output.txt" >> $slurmfile
## Main code
echo "cd \${SLURM_SUBMIT_DIR}" >> $slurmfile
echo "echo \"This is job number \$SLURM_JOB_ID in \$SLURM_SUBMIT_DIR running command number \$SLURM_ARRAY_TASK_ID on node \$SLURM_NODEID\" &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
echo "date  &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
echo "cmd=\`head -n \$SLURM_ARRAY_TASK_ID $cmdfile | tail -n 1\`">> $slurmfile
echo "echo \"Running: \$cmd\"  &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
echo "eval \$cmd  &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
echo "echo \"Completed...\"  &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
echo "date  &>> $log/\$SLURM_JOB_ID.out" >> $slurmfile
if [ "$pretend" == 0 ] ; then
    echo "Submitting sbatch script $slurmfile"
    sbatch $slurmfile
else
    echo "Pretend mode; not submitting sbatch script $slurmfile"
    echo "Command to submit:"
    echo "sbatch $slurmfile"
fi
