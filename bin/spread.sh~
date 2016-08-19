for file in `ls /home/markov/oliveira/saint2_orig/scripts/5PTI_c_n_lib_5pti_10000_1000_t2.5_linear/0/decoy*`
do
	j=0
	for file2 in `ls /home/markov/oliveira/saint2_orig/scripts/5PTI_c_n_lib_5pti_10000_1000_t2.5_linear/0/decoy*`
	do
		/home/markov/oliveira/saint2_orig/3rdparty/TMscore $file $file2 | awk '/^TM-score    =/ { printf "%f\n", $3; }'
	done
done
