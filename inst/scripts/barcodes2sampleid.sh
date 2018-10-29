#!/usr/bin/env bash
#
# barcodes2sampleid.sh
#
# Gerhard Schöfl, 2016
#
# USAGE:
# barcodes2sampleid [targetdirectory] [conversion.csv]
#
# conversion.csv with columns <BC,SAMPLEID,LOCUS> in exactly that
# order and exactly those names!
#
# Barcodes in a single conversion file must be unique
#

## GLOBALS ###################################################################

## consts ##
declare -r TARGETDIR="$1"
declare -r CONVCSV="$2"
declare -r SEPARATOR="§"
declare -r BACBIO_DIR="${TARGETDIR}/pacbio"
declare -r NANOPORE_DIR="${TARGETDIR}/nanopore"
declare -r ILLUMINA_DIR="${TARGETDIR}/illumina"

if [ -z "$TARGETDIR" -o "$TARGETDIR" = "-h" ]; then
    echo ""
    echo "USAGE:"
    echo "barcodes2sampleid [targetdirectory] [conversion.csv]"
    echo ""
    echo "conversion.csv with columns <BC,SAMPLEID,LOCUS> in"
    echo "exactly that order and exactly those names!"
    echo ""
    exit 0
fi

## File regexps ##

## Expectation, e.g.: "lbc10--lbc10.fastq"
declare -r PACBIO_REGEX="^lbc[0-9]+--lbc[0-9]+.+fastq"
## Expectation, e.g.: "ONBC_009_ONT1_A_2D.fastq"
declare -r NANOPORE_REGEX="^ONBC_[0-9]+.+fastq"
## Expectation, e.g.: "9_S9_L001_R1_001.fastq" and "9_S9_L001_R2_001.fastq"
declare -r ILLUMINA_REGEX="^[0-9]+.+R[12].+fastq"

## arrays ##
declare -a BC
declare -a SAMPLEID
declare -a LOCUS
declare -a TABLE

## integers ##
declare -i INDEX=0

## chars ##
declare FASTQ_OLD
declare FASTQ_NEW

## FUNCS #####################################################################
function add_to_table() {
	# @param $1 barcode (KEY)
	# @param $2 sampleid (VAL1)
	# @param $3 locus (VAL2)
	TABLE[$INDEX]="$1"$SEPARATOR"$2"$SEPARATOR"$3"
	#printf "Added %s => %s => %s at [%s]\n" $1 $2 $3 $INDEX
	let "INDEX++"
}
declare -rf add_to_table

function lookup() {
	# @param $1 barcode (KEY)
	# @return $2 result (SAMPLEID)
	# @return $3 result (LOCUS)
	local -i __elem=0
	local __key
	local __val
	local __sampleid="$2"
	local __locus="$3"
	## search for key
	while [ $__elem -lt $INDEX ]; do
		__key="${TABLE[$__elem]}"
		__key="${__key%%$SEPARATOR*}"
		if [ "$__key" -eq "$1" ]; then
			break
		else
			let "__elem++"
		fi
	done
	## Extract value
	__val="${TABLE[$__elem]}"
	__val="${__val#*$SEPARATOR}"
	local __tmp=()
	__tmp=(${__val//$SEPARATOR/ })
	eval $__sampleid='"${__tmp[0]}"'
	eval $__locus='"${__tmp[1]}"'
}
declare -rf lookup

function barcode_from_bacbio_file() {
	# @param  $1 bacbio file e.g. "lbc10--lbc10.fastq"
	# @return $2 barcode "10"
	# @return $3 tail "lbc10.fastq"
	if [[ "$1" =~ $PACBIO_REGEX ]]; then
		local __ret1="$2"
		local __ret2="$3"
		local __tmp0="${1%%--lbc*}"
		local __tmp1="${__tmp0#lbc*}" ## barcode
		local __tmp2="${1#*--}"       ## tail
		eval $__ret1="'$__tmp1'"
		eval $__ret2="'$__tmp2'"
		return 0 ## all is well
	else
		printf "Couldn't recognise PacBio file format: %s\n" $1
		return 1 ## not all is well
	fi
}
declare -rf barcode_from_bacbio_file

function barcode_from_nanopore_file() {
	# @param  $1 nanopore file e.g. "ONBC_009_ONT1_A_2D.fastq"
	# @return $2 barcode "9"
	# @return $3 tail "ONT1_A_2D.fastq"
	if [[ "$1" =~ $NANOPORE_REGEX ]]; then
		local __ret1="$2"
		local __ret2="$3"
		local __tmp0="${1#*_}"
		local __tmp1="${__tmp0%%_*}"	## barcode
		#local __tmp2="${__tmp0#*_}"	## tail
		local __tmp2="$__tmp0"			## tail
		eval $__ret1="'$__tmp1'"
		eval $__ret2="'$__tmp2'"
		return 0 ## all is well
	else
		printf "Couldn't recognise Nanopore file format: %s\n" $1
		return 1 ## not all is well
	fi
}
declare -rf barcode_from_bacbio_file

function barcode_from_illumina_file() {
	# @param  $1 illumina file e.g. "9_S9_L001_R1_001.fastq"
	# @return $2 barcode "9"
	# @return $3 tail "S9_L001_R1_001.fastq"
	if [[ "$1" =~ $ILLUMINA_REGEX ]]; then
		local __ret1="$2"
		local __ret2="$3"
		local __tmp1="${1%%_*}"		## barcode
		local __tmp2="${1#*_}"		## tail
		eval $__ret1="'$__tmp1'"
		eval $__ret2="'$__tmp2'"
		return 0 ## all is well
	else
		printf "Couldn't recognise Illumina file format: %s\n" $1
		return 1 ## not all is well
	fi
}
declare -rf barcode_from_bacbio_file

function move_pacbio_files() {
	local __file
	local __bc
	local __tail
	local __sid
	local __loc
	local __old
	local __new
	for __file in $(ls $BACBIO_DIR); do
		__old="${BACBIO_DIR}/${__file}"
		if barcode_from_bacbio_file $__file __bc __tail; then
			lookup $__bc __sid __loc
			__new="${BACBIO_DIR}/${__sid}_${__loc}_${__tail}"
			printf "Renaming %s => %s\n" $(basename $__old) $(basename $__new)
			mv -n "$__old" "$__new"
		else
			printf "Cannot rename %s\n" $(basename $__old)
		fi
	done
}
declare -rf move_pacbio_files

function move_nanopore_files() {
	local __file
	local __bc
	local __tail
	local __sid
	local __loc
	local __old
	local __new
	for __file in $(ls $NANOPORE_DIR); do
		__old="${NANOPORE_DIR}/${__file}"
		if barcode_from_nanopore_file $__file __bc __tail; then
			lookup $__bc __sid __loc
			__new="${NANOPORE_DIR}/${__sid}_${__loc}_${__tail}"
			printf "Renaming %s => %s\n" $(basename $__old) $(basename $__new)
			mv -n "$__old" "$__new"
		else
			printf "Cannot rename %s\n" $(basename $__old)
		fi
	done
}
declare -rf move_nanopore_files

function move_illumina_files() {
	local __file
	local __bc
	local __tail
	local __sid
	local __loc
	local __old
	local __new
	for __file in $(ls $ILLUMINA_DIR); do
		__old="${ILLUMINA_DIR}/${__file}"
		if barcode_from_illumina_file $__file __bc __tail; then
			lookup $__bc __sid __loc
			__new="${ILLUMINA_DIR}/${__sid}_${__loc}_${__tail}"
			printf "Renaming %s => %s\n" $(basename $__old) $(basename $__new)
			mv -n "$__old" "$__new"
		else
			printf "Cannot rename %s\n" $(basename $__old)
		fi
	done
}
declare -rf move_illumina_files

## Main ######################################################################

## read CSV
OLDIFS=$IFS
IFS=","
[ ! -f $CONVCSV ] && { echo "$CONVCSV file not found"; exit 1; }
while read bc__ sid__ loc__ __; do
	BC+=("$bc__")
	SAMPLEID+=("$sid__")
	LOCUS+=("$loc__")
done < "$2"

## Sanity checks
echo "<${BC[0]}>"
if [ "${BC[0]}" != "BC" ]; then
	echo "The first column in the conversion file must contain barcodes and be named 'BC'"
	exit 1
fi

echo "<${SAMPLEID[0]}>"
if [ "${SAMPLEID[0]}" != "SAMPLEID" ]; then
	echo "The second column in the conversion file must contain sample ids and be named 'SAMPLEID'"
	exit 1
fi

echo "<${LOCUS[0]}>"
if [ "${LOCUS[0]}" != "LOCUS" ]; then
	echo "The third column in the conversion file must contain locus names and be named 'LOCUS'"
	exit 1
fi

unset bc__
unset sid__
unset loc__
IFS=$OLDIFS

## create lookup table:
##   BC => SAMPLEID, LOCUS
for ((i=1; i<${#BC[*]}; i++)); do
    add_to_table ${BC[i]} ${SAMPLEID[i]} ${LOCUS[i]}
done

## move bacbio files
if [ -d "$BACBIO_DIR" ]; then
	echo "Entering 'pacbio' directory"
	move_pacbio_files
else
	echo "Couldn't find 'pacbio' directory"
fi

## move nanopore files
if [ -d "$NANOPORE_DIR" ]; then
	echo "Entering 'nanopore' directory"
	move_nanopore_files
else
	echo "Couldn't find 'nanopore' directory"
fi

## move illumina files
if [ -d "$ILLUMINA_DIR" ]; then
	echo "Entering 'illumina' directory"
	move_illumina_files
else
	echo "Couldn't find 'illumina' directory"
fi

exit 0 # all is well
