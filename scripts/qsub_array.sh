#!/bin/bash

function show_help {
    echo "Usage: ./qsub_array.sh -f <cmdfile>"
    echo "  Submit commands to a qsub system using an array job, where <cmdfile> is a file containing a list of commands to execute, one per line."
    echo "OPTIONS:"
    echo "  -f <val>: The name of the command file to submit. REQUIRED."
    echo "  -n <val>: Set the number of processes per command in the qsub script. (Default: 1)"
    echo "  -w <val>: Set the walltime requested. (Default: 120, measured in hours)"
    echo "  -P <val>: Set the number of processes per node per command in the qsub script. (Default: 1)"
    echo "  -p      : Pretend: do not actually run the commands, but generate the qsub scripts for examination."
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

pbsfile="$cmdfile.pbs"
if [ ! -f $cmdfile ]; then
    echo "Command file $cmdfile does not exist!"
    exit 1
fi
ncmds=`cat $cmdfile | wc -l`
log=`mktemp -d $HOME/log/${cmdfile}.XXXXXX`
echo "Making qsub script $pbsfile from $cmdfile with $ncmds commands"
echo "#!/bin/bash -login" > $pbsfile
echo "#PBS -l walltime=$walltime:00:00,nodes=$nodes:ppn=$ppnode"  >> $pbsfile
echo "#PBS -t 1-$ncmds" >> $pbsfile
## PBS output log is broken so instead of this, we workaround by redirecting all output
##echo "#PBS -o ~/log/\$PBS_JOBID/log.txt" >> $pbsfile 
##echo "#PBS -e ~/log/\$PBS_JOBID/error.txt" >> $pbsfile
##echo "mkdir -p ~/log/\$PBS_JOBID" >> $pbsfile
##
echo "cd \$PBS_O_WORKDIR" >> $pbsfile
echo "echo \"This is job number \$PBS_JOBID in \$PBS_O_WORKDIR running command number \$PBS_ARRAYID on node \$PBS_NODEFILE\" &>> $log/\$PBS_JOBID.out" >> $pbsfile
echo "date  &>> $log/\$PBS_JOBID.out" >> $pbsfile
echo "cmd=\`head -n \$PBS_ARRAYID $cmdfile | tail -n 1\`">> $pbsfile
echo "echo \"Running: \$cmd\"  &>> $log/\$PBS_JOBID.out" >> $pbsfile
echo "eval \$cmd  &>> $log/\$PBS_JOBID.out" >> $pbsfile
echo "echo \"Completed...\"  &>> $log/\$PBS_JOBID.out" >> $pbsfile
echo "date  &>> $log/\$PBS_JOBID.out" >> $pbsfile
if [ "$pretend" == 0 ] ; then
    echo "Submitting qsub script $pbsfile"
    qsub $pbsfile
else
    echo "Pretend mode; not submitting qsub script $pbsfile"
    echo "Command to submit:"
    echo "qsub $pbsfile"
fi
