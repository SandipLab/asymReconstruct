# Sample script for the sorting and angular reorientation of cyclic pseudosymmetric particles in Cryo-EM data processing
# Written by Zhuowen Li
# Sandip Basak Lab @ NISB, Nanyang Technological University
# 
# Usage: bash [script.sh] [input.star] [classes] [Ns] [configuration] [output.star]
#
# input.star          The run_data.star file from a RELION 3D classification.
# classes             An n-letter string where n is the number of classes in the input.star, the letter on the m-th position denotes the conformation indicator of _rlnClassNumber #m in the corresponding input.star, e.g. "aabbbb" stands for "a" conformation for class 1 and class 2 while "b" conformation assigned to class 3 to 6.
# Ns                  Number of subunits for each cyclically symmetric particle.
# configuration       An Ns-letter string standing for the target particle configuration of subunit conformations. This script is designated for sorting and aligning cyclically organized oligomeric single particles with its symmetry-applied reconstruction performed in RELION software. The letters used in this string should conform to those defined by the above "classes" argument. 
# output.star         Output file name to be specified.

# Print Usage
if [ $# -lt 5 ] || [ $1 = "-h" ] || [ $1 = "--help" ]; then 
	head -11 $0 | tail -7
	exit 1
fi

# Identify column number and names of input.star
# Required columns:
# 	_rlnMicrographName
#	_rlnCoordinateX
#	_rlnCoordinateY
#	_rlnAngleRot
#	_rlnClassNumber
#
col_id() {
	awk -v query=$1 '{if(index($1,query)>0){print substr($2,2); exit}}' $2
}
number_cols() {
	awk '{cnt[NF]+=1} END{for (num in cnt) {if(cnt[num]>maxCnt) {maxCnt=cnt[num]; maxCol=num}}; print maxCol}' $1
}

# Sort and Align
awk -v class_assignments=$2 -v Ns=$3 -v configuration=$4 -v col_Micrograph=$(col_id _rlnMicrographName $1) -v col_X=$(col_id _rlnCoordinateX $1) -v col_Y=$(col_id _rlnCoordinateY $1) -v col_Rot=$(col_id _rlnAngleRot $1) -v col_Class=$(col_id _rlnClassNumber $1) -v Ncol=$(number_cols $1) '
function as_string(as_list, as_ID, as_Ns){
	as_str=""
	for (as_i=0;as_i<as_Ns;as_i++){
		as_str=as_str as_list[as_ID " " as_i]
	}
	return as_str
}
function angle_to_rank(atr_ang, atr_Ns){
	if(atr_ang>=0) return int(atr_ang/(360/atr_Ns))
	else return int(atr_ang/(360/atr_Ns))+(atr_Ns-1)
}
function is_circularly_equal(ice_str1,ice_str2){
	for(ice_i=1;ice_i<=length(ice_str1);ice_i++){
		if(ice_str2 == substr(ice_str1,ice_i) substr(ice_str1,1,ice_i-1)) return ice_i-1
	}
	return -1
}
BEGIN{
	# ptl_ID               = $col_Micrograph " " $col_X " " $col_Y
	# subunit_ID           = $col_Micrograph " " $col_X " " $col_Y " " $col_Rot
	# subunit_conformation = substr(class_assignments, $col_Class,1)
	# subunit_rank         = angle_to_rank($col_Rot)
	# record               = $0
}
{
	if($Ncol==""){
		print
	}
	if($Ncol!=""){
		ptl_list[$col_Micrograph " " $col_X " " $col_Y]=1
		subunits[$col_Micrograph " " $col_X " " $col_Y " " $col_Rot] = $0

		ptl_configuration[$col_Micrograph " " $col_X " " $col_Y " " angle_to_rank($col_Rot, Ns)]=substr(class_assignments, $col_Class, 1)
		ptl_ang_at_rank[$col_Micrograph " " $col_X " " $col_Y " " angle_to_rank($col_Rot, Ns)]=$col_Rot
	}
}
END{
	for (ptl in ptl_list){
		if(is_circularly_equal(as_string(ptl_configuration, ptl, Ns), configuration)!=-1){
			print subunits[ptl " " ptl_ang_at_rank[ptl " " is_circularly_equal(as_string(ptl_configuration, ptl,Ns), configuration)]]
		}
	}
}' $1 > $5

# Exit state and stats
nHeader=$(awk -v maxColEntry=$(number_cols $1) '{if(NF<maxColEntry) cntHeader+=1} END{print cntHeader}' $1)
nOut=$(wc -l $5 | awk '{print $1}')
if [ $nOut -gt $nHeader ]; then
	echo "Alignment is done. Number of aligned particles with compositon $4 is $(($nOut - $nHeader))."
fi