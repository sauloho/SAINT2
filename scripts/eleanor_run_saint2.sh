#!/bin/bash

#export SAINT2=$HOME/saint2_orig/
SCORING=false
OUTPUT=$1
METHODS=$2
cd $3 > /dev/null
export DATA_PATH=`pwd | tr -d '\n'`
cd - > /dev/null
SEGMENT_LENGTH=$4
LENGTH=`cat $DATA_PATH/$OUTPUT.length | tr -d '\n'`
INITIAL_LENGTH=9
GROWTH_MOVES=`bc <<< "10000*($LENGTH-$SEGMENT_LENGTH)^2/($LENGTH)^2"`
MOVES=`bc <<< "1000*($LENGTH-$SEGMENT_LENGTH)/$LENGTH"`
#SNAPSHOTS=$5
HOST=`hostname`

if [ $SCORING = true ]
then
	mkdir -p $HOST/$OUTPUT
fi

if [[ $METHODS =~ c ]]
then
	#COTRANS

	# set up directory for output
	OUTPATH=${OUTPUT}_c_n_${OUTPUT}.flib_${GROWTH_MOVES}_${MOVES}_t2.5_Seg${SEGMENT_LENGTH}_linear/$HOST
	mkdir -p $OUTPATH

	# make a cut version of the PDB for the non-segment part if it doesn't already exist
	CUTPDB=${DATA_PATH}/${OUTPUT}_inverseSeg${SEGMENT_LENGTH}_fwd.pdb
	if [ ! -a $CUTPDB ]
	then
		awk -v seglength="$SEGMENT_LENGTH" '$1 == "ATOM" && $6 > seglength {print $0}' $DATA_PATH/$OUTPUT.pdb > $CUTPDB
	fi

	### run saint2 and store the filename of the decoy generated on the variable "$FILE" ###
	FILE=`$SAINT2/scripts/eleanor_run_cotrans2 1 $OUTPUT n $DATA_PATH/$OUTPUT.fasta.txt $DATA_PATH/$OUTPUT.flib $GROWTH_MOVES $MOVES 2.5 segment $SEGMENT_LENGTH linear`
	awk -v seglength="$SEGMENT_LENGTH" '$1 == "ATOM" && $6 > seglength {print $0}' $OUTPATH/$FILE > $OUTPATH/$FILE.cut

	# get the last part of the filename after the last underscore (process ID)
	PID=${FILE##*_}

	### get the scores of the final decoy that has been generated ###
	#TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')

	if [ $SCORING = true ]
	then
		# get the TM score of the non-segment part against the non-segment part of the pdb
		TMCUT=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE.cut $CUTPDB | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
		TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
		#rm $OUTPATH/$FILE.cut
	fi
	
	ERRFILE=$OUTPATH/err_${HOST}_$PID
	paste <(grep -e "^MOVE LOCATION:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE FRAGMENT:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE FRAGLENG:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE ACCEPTED:" $ERRFILE | cut -f 2) \
		> $OUTPATH/moves_$PID
	grep -ve "^MOVE " $ERRFILE > $OUTPATH/err_${HOST}_$PID.temp
	if [ "$?" = "0" ]; then
		mv $OUTPATH/err_${HOST}_$PID.temp $OUTPATH/err_${HOST}_$PID
	fi

	rm -f $OUTPATH/scmatrix*
	# Package up the part and perc files into MODELs
	#for FRAC in perc part; do
	#for FRAC in part; do
	#	COUNT=1
	#	for i in `ls -v $OUTPATH/${FILE}_${FRAC}*`; do
	#		echo "MODEL       " $COUNT >> $OUTPATH/${FILE}_${FRAC}
	#		grep "^ATOM" $i >> $OUTPATH/${FILE}_${FRAC}
	#		rm $i
	#		echo ENDMDL >> $OUTPATH/${FILE}_${FRAC}
	#		((COUNT++))
	#	done
	#done

	if [ $SCORING = true ]
	then
		# score the decoy using the scoring or running config file
		#$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_scoring -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
		$SAINT2/bin/saint2 config_${OUTPUT}_c_n*_Seg${SEGMENT_LENGTH}_linear -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
		if [ "$?" = "0" ]; then
			SOLV=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Solvation =/ { print $NF; }'`
			ORIE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Orientation =/ { print $NF; }'`
			LJ=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Lennard-Jones =/ { print $NF; }'`
			RAPDF=`cat $HOST/$OUTPUT/temp_$$ | awk '/^RAPDF =/ { print $NF; }'`
			CORE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^CORE =/ { print $NF; }'`
			PREDSS=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredSS =/ { print $NF; }'`
			SAULO=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Saulo =/ { print $NF; }'`
			ELEANOR=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Eleanor =/ { print $NF; }'`
			COMB=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Combined score =/ { print $NF; }'`
			RG=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Radius of gyration =/ { print $NF; }'`
			DIAM=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Diameter =/ { print $NF; }'`
			TORSION=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredTor =/ { print $NF; }'`
			if [ -n "$TM" ]; then
				#echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $PREDSS $SAULO $ELEANOR $RG $DIAM $TORSION $COMB >> $HOST/$OUTPUT/scores_cotrans_$$.txt
				#echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $SAULO $ELEANOR $RG $DIAM $COMB >> $HOST/$OUTPUT/scores_cotrans_$PID.txt
				echo $TM $TMCUT $FILE $SOLV $ORIE $RAPDF $LJ $CORE $SAULO $ELEANOR $RG $DIAM $COMB >> $HOST/$OUTPUT/scores_cotrans_$PID.txt
			fi
		fi
		rm $HOST/$OUTPUT/temp_$$
	fi
fi

if [[ $METHODS =~ v ]]
then
	#REVERSE

	# set up directory for output
	OUTPATH=${OUTPUT}_c_n_${OUTPUT}.flib_${GROWTH_MOVES}_${MOVES}_t2.5_Seg${SEGMENT_LENGTH}_linear_rev/$HOST
	mkdir -p $OUTPATH

	# make a cut version of the PDB for the non-segment part if it doesn't already exist
	CUTPDB=${DATA_PATH}/${OUTPUT}_inverseSeg${SEGMENT_LENGTH}_rev.pdb
	if [ ! -a $CUTPDB ]
	then
		awk -v seglength="`expr $LENGTH - $SEGMENT_LENGTH`" '$1 == "ATOM" && $6 < seglength {print $0}' $DATA_PATH/$OUTPUT.pdb > $CUTPDB
	fi

	### run saint2 and store the filename of the decoy generated on the variable "$FILE" ###
	FILE=`$SAINT2/scripts/eleanor_run_cotrans2 1 $OUTPUT n $DATA_PATH/$OUTPUT.fasta.txt $DATA_PATH/$OUTPUT.flib $GROWTH_MOVES $MOVES 2.5 segment $SEGMENT_LENGTH linear rev`
	awk -v seglength="`expr $LENGTH - $SEGMENT_LENGTH`" '$1 == "ATOM" && $6 < seglength {print $0}' $OUTPATH/$FILE > $OUTPATH/$FILE.cut

	# get the last part of the filename after the last underscore (process ID)
	PID=${FILE##*_}

	### get the scores of the final decoy that has been generated ###
	#TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
	
	if [ $SCORING = true ]
	then
		# get the TM score of the non-segment part against the non-segment part of the pdb
		TMCUT=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE.cut $CUTPDB | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
		TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
		#rm $OUTPATH/$FILE.cut
	fi
	
	ERRFILE=$OUTPATH/err_${HOST}_$PID
	paste <(grep -e "^MOVE LOCATION:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE FRAGMENT:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE FRAGLENG:" $ERRFILE | cut -f 2) \
		<(grep -e "^MOVE ACCEPTED:" $ERRFILE | cut -f 2) \
		> $OUTPATH/moves_$PID
	grep -ve "^MOVE " $ERRFILE > $OUTPATH/err_${HOST}_$PID.temp
	if [ "$?" = "0" ]; then
		mv $OUTPATH/err_${HOST}_$PID.temp $OUTPATH/err_${HOST}_$PID
	fi

	rm -f $OUTPATH/scmatrix*
	# Package up the part and perc files into MODELs
	#for FRAC in perc part; do
	#for FRAC in part; do
	#	COUNT=1
	#	for i in `ls -v $OUTPATH/${FILE}_${FRAC}*`; do
	#		echo "MODEL       " $COUNT >> $OUTPATH/${FILE}_${FRAC}
	#		grep "^ATOM" $i >> $OUTPATH/${FILE}_${FRAC}
	#		rm $i
	#		echo ENDMDL >> $OUTPATH/${FILE}_${FRAC}
	#		((COUNT++))
	#	done
	#done

	if [ $SCORING = true ]
	then
		# score the decoy using the scoring or running config file
		#$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_scoring -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
		$SAINT2/bin/saint2 config_${OUTPUT}_c_n*_Seg${SEGMENT_LENGTH}_linear_rev -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
		if [ "$?" = "0" ]; then
			SOLV=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Solvation =/ { print $NF; }'`
			ORIE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Orientation =/ { print $NF; }'`
			LJ=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Lennard-Jones =/ { print $NF; }'`
			RAPDF=`cat $HOST/$OUTPUT/temp_$$ | awk '/^RAPDF =/ { print $NF; }'`
			CORE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^CORE =/ { print $NF; }'`
			PREDSS=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredSS =/ { print $NF; }'`
			SAULO=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Saulo =/ { print $NF; }'`
			ELEANOR=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Eleanor =/ { print $NF; }'`
			COMB=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Combined score =/ { print $NF; }'`
			RG=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Radius of gyration =/ { print $NF; }'`
			DIAM=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Diameter =/ { print $NF; }'`
			TORSION=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredTor =/ { print $NF; }'`
			if [ -n "$TM" ]; then
				#echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $PREDSS $SAULO $ELEANOR $RG $DIAM $TORSION $COMB >> $HOST/$OUTPUT/scores_cotrans_$$.txt
				#echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $SAULO $ELEANOR $RG $DIAM $COMB >> $HOST/$OUTPUT/scores_reverse_$PID.txt
				echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $SAULO $RG $DIAM $COMB >> $HOST/$OUTPUT/scores_reverse_$PID.txt
			fi
		fi
		rm $HOST/$OUTPUT/temp_$$
	fi
fi

if [[ $METHODS =~ r ]]
then
	#REVERSE

	# set up directory for output
	OUTPATH=${OUTPUT}_c_n_${OUTPUT}.flib_10000_1000_t2.5_linear_rev/$HOST
	mkdir -p $OUTPATH

	### run saint2 and store the FILENAME OF THE decoy generated on the variable "$FILE" ###
	FILE=`$SAINT2/scripts/eleanor_run_cotrans2 1 $OUTPUT n $DATA_PATH/$OUTPUT.fasta.txt $DATA_PATH/$OUTPUT.flib 10000 1000 2.5 linear rev`

	# get the last part of the filename after the last underscore (process ID)
	PID=${FILE##*_}

	### get the scores of the final decoy that has been generated ###
	TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
	
	# Package up the part and perc files into MODELs
	for FRAC in perc part; do
		COUNT=1
		for i in `ls -v $OUTPATH/${FILE}_${FRAC}*`; do
			echo "MODEL       " $COUNT >> $OUTPATH/${FILE}_${FRAC}
			grep "^ATOM" $i >> $OUTPATH/${FILE}_${FRAC}
			rm $i
			echo ENDMDL >> $OUTPATH/${FILE}_${FRAC}
			((COUNT++))
		done
	done

	# score the decoy using the scoring config file
	#$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_scoring -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
	$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_c_n*linear_rev -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
	if [ "$?" = "0" ]; then
		SOLV=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Solvation =/ { print $NF; }'`
		ORIE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Orientation =/ { print $NF; }'`
		LJ=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Lennard-Jones =/ { print $NF; }'`
		RAPDF=`cat $HOST/$OUTPUT/temp_$$ | awk '/^RAPDF =/ { print $NF; }'`
		CORE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^CORE =/ { print $NF; }'`
		PREDSS=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredSS =/ { print $NF; }'`
		SAULO=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Saulo =/ { print $NF; }'`
		COMB=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Combined score =/ { print $NF; }'`
		RG=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Radius of gyration =/ { print $NF; }'`
		DIAM=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Diameter =/ { print $NF; }'`
		TORSION=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredTor =/ { print $NF; }'`
		if [ -n "$TM" ]; then
			echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $PREDSS $SAULO $RG $DIAM $TORSION $COMB >> $HOST/$OUTPUT/scores_reverse_$$.txt
		fi
	fi
	rm $HOST/$OUTPUT/temp_$$
fi

if [[ $METHODS =~ i ]]
then
	#INVITRO

	# set up directory for output
	OUTPATH=${OUTPUT}_i_${OUTPUT}.flib_11000_t2.5/$HOST
	mkdir -p $OUTPATH

	### run saint2 and store the filename of the decoy generated on the variable "$FILE" ###
	FILE=`$SAINT2/scripts/eleanor_run_invitro2 1 $OUTPUT $DATA_PATH/$OUTPUT.fasta.txt $DATA_PATH/$OUTPUT.flib 11000 2.5` 

	# get the last part of the filename after the last underscore (process ID)
	PID=${FILE##*_}

	### get the scores of the final decoy that has been generated ###
	TM=$($SAINT2/3rdparty/TMalign $OUTPATH/$FILE $DATA_PATH/$OUTPUT.pdb | grep -m 1 TM-score= | awk '{ printf "%f",$2; }')
	
	# score the decoy using the scoring config file
	#$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_scoring -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
	$SAINT2/bin/saint2 $DATA_PATH/config_${OUTPUT}_i_* -- $OUTPATH/$FILE > $HOST/$OUTPUT/temp_$$
	if [ "$?" = "0" ]; then
		SOLV=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Solvation =/ { print $NF; }'`
		ORIE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Orientation =/ { print $NF; }'`
		LJ=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Lennard-Jones =/ { print $NF; }'`
		RAPDF=`cat $HOST/$OUTPUT/temp_$$ | awk '/^RAPDF =/ { print $NF; }'`
		CORE=`cat $HOST/$OUTPUT/temp_$$ | awk '/^CORE =/ { print $NF; }'`
		PREDSS=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredSS =/ { print $NF; }'`
		SAULO=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Saulo =/ { print $NF; }'`
		COMB=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Combined score =/ { print $NF; }'`
		RG=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Radius of gyration =/ { print $NF; }'`
		DIAM=`cat $HOST/$OUTPUT/temp_$$ | awk '/^Diameter =/ { print $NF; }'`
		TORSION=`cat $HOST/$OUTPUT/temp_$$ | awk '/^PredTor =/ { print $NF; }'`
		if [ -n "$TM" ]; then
			echo $TM $FILE $SOLV $ORIE $RAPDF $LJ $CORE $PREDSS $SAULO $RG $DIAM $TORSION $COMB >> $HOST/$OUTPUT/scores_invitro_$$.txt
		fi
	fi
	rm $HOST/$OUTPUT/temp_$$
fi

