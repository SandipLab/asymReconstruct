# Sample script for configuration statistics of cyclic pseudosymmetric proteins in Cryo-EM data processing
# Written by Zhuowen Li
# Sandip Basak Lab @ NISB, Nanyang Technological University
# 
# Usage: bash asymmetric_conformations_stats.sh [input.star] [classes] [Ns]
#
# input.star          The run_data.star file from a RELION 3D classification.
# classes             An n-letter string where n is the number of classes in the input.star, the letter on the m-th position denotes the conformation indicator of _rlnClassNumber #m in the corresponding input.star, e.g. "aabbbb" stands for "a" conformation for class 1 and class 2 while "b" conformation assigned to class 3 to 6.
# Ns                  Number of subunits for each cyclically symmetric particle.
#
#
#
#
#
#
# Print Usage
if [ $# -lt 3 ] || [ $1 = "-h" ] || [ $1 = "--help" ]; then
	head -9 $0 | tail -5
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
#
# Sort and Align
awk -v class_assignments=$2 -v Ns=$3 -v col_Micrograph=$(col_id _rlnMicrographName $1) -v col_X=$(col_id _rlnCoordinateX $1) -v col_Y=$(col_id _rlnCoordinateY $1) -v col_Rot=$(col_id _rlnAngleRot $1) -v col_Class=$(col_id _rlnClassNumber $1) -v Ncol=$(number_cols $1) '
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
	#if($Ncol==""){
	#	print
	#}
	if($Ncol!=""){
		ptl_list[$col_Micrograph " " $col_X " " $col_Y]=1
		subunits[$col_Micrograph " " $col_X " " $col_Y " " $col_Rot] = $0

		ptl_composition[$col_Micrograph " " $col_X " " $col_Y " " angle_to_rank($col_Rot, Ns)]=substr(class_assignments, $col_Class, 1)
		ptl_ang_at_rank[$col_Micrograph " " $col_X " " $col_Y " " angle_to_rank($col_Rot, Ns)]=$col_Rot
	}
}
END{
	num_compositions=0
	for (ptl in ptl_list){
		query_composition=as_string(ptl_composition, ptl, Ns)
		if(num_compositions==0){
			num_compositions+=1
			compositions[num_compositions]=query_composition
			compositions_cnt[num_compositions]+=1
		}
		else{
			flag_found=0
			for(qc_i=1;qc_i<=num_compositions;qc_i++){
				if(is_circularly_equal(query_composition, compositions[qc_i])!=-1){
					compositions_cnt[qc_i]+=1
					flag_found=1
					break
				}
			}
			if(flag_found==0){
				num_compositions+=1
				compositions[num_compositions]=query_composition
				compositions_cnt[num_compositions]+=1
			}
		}

	}
	# print out the statistical result
	for(pri=1;pri<=num_compositions;pri++){
		print compositions[pri], compositions_cnt[pri]
	}
}' $1
