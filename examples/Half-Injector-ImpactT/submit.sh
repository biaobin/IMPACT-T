#PBS -N half-injector
#PBS -l nodes=1:ppn=4
#PBS -j oe 
#PBS -l walltime=24:00:00 
#PBS -q  queue_46


cd $PBS_O_WORKDIR 

echo my job id is $PBS_JOBID | tee mpi.log
echo run nodes is following: | tee -a mpi.log
cat $PBS_NODEFILE | tee -a mpi.log

echo begin time is `date` | tee -a mpi.log
id=`echo $PBS_JOBID|awk -F. '{print $1}'`

NP=`cat $PBS_NODEFILE|wc -l` 
exe=/public/home/biaobin/gitproj/IMPACT-T/src/ImpactT.exe

mpirun -np $NP $exe

echo end time is `date` | tee -a mpi.log

