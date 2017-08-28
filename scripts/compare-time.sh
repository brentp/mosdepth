set -euo pipefail
export PATH=/scratch/ucgd/lustre/work/u6000771/Projects/src/samtools/:/scratch/ucgd/lustre/work/u6000771/Projects/src/bedtools2/bin/:$PATH
export LD_LIBRARY_PATH=/scratch/ucgd/lustre/work/u6000771/Projects/src/htslib/

export bam=/scratch/ucgd/lustre/work/u6000771/Data/NA12878.mem.bam
export bam=/scratch/ucgd/lustre/work/u6000771/Data/mosdepth/ERR1395576.30X.bam
#export bam=/scratch/ucgd/lustre/work/u6000771/Data/mosdepth/small.bam
export cram=/scratch/ucgd/lustre/work/u6000771/Data/mosdepth/ERR1395576.30X.cram
export fasta=/scratch/ucgd/lustre/work/u6000771/bcbio/genomes/Hsapiens/g1k_v37_decoy/seq/g1k_v37_decoy.fa
#export cram=/scratch/ucgd/lustre/work/u6000771/Data/NA12878.mem.cram
#export bam=test.bam

export TIMEFORMAT='%lU'

dotime() {
	# returns the max memory use in MiB and time in seconds for the given command.
	# usage:
	# sets _MEM, _TIME env vars.
	# so use can also be:
    # $ dotime "sleep 2" > /dev/null
	# $ echo  $_MEM $_TIME
	set -euo pipefail
	export cmd="$1"
	T="$(date +%s)"
	res="$(/usr/bin/time -v bash -c "$cmd 2>xx.stderr" 2>&1 1>/dev/null)" # | grep Maximum | awk '{ print $NF}' )"
	_MEM=$(echo "$res" | grep Maximum | awk '{ print $NF }')
	export _MEM=$(echo "scale=2; $_MEM / 1024" | bc -l)
	export _USER_TIME="$(echo "$res" | grep "User time" | awk '{ print $NF }')"
	export _WALL_TIME="$(echo "$res" | grep "Elapsed " | awk '{ print $NF }')"
}

report_time() {
	dotime "$1"
	echo -e "$1|$_USER_TIME|$_WALL_TIME|$_MEM"
}



report_time "sleep 1"
report_time "mosdepth -t 2 $cram"
report_time "mosdepth -t 1 $cram"
report_time "mosdepth -t 0 $cram"

report_time "mosdepth -t 2 $bam"
report_time "mosdepth -t 1 $bam"
report_time "mosdepth -t 0 $bam"

report_time "mosdepth -t 2 --by 500 $bam"
report_time "mosdepth -t 1 --by 500 $bam"
report_time "mosdepth -t 0 --by 500 $bam"

report_time "mosdepth -t 2 --by 500 $cram"
report_time "mosdepth -t 1 --by 500 $cram"
report_time "mosdepth -t 0 --by 500 $cram"
report_time	"sambamba_v0.6.6 depth base --fix-mate-overlaps $bam"

report_time "samtools depth -d 100000 $bam"
report_time "samtools depth -d 100000 $cram"
report_time "bedtools genomecov -ibam $bam"
