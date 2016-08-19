if [ -z "$SAINT2" ]
then
	export SAINT2=/home/markov/oliveira/saint2_orig
fi

	$SAINT2/scripts/run_cotrans2 200 1QUU n /homes/oliveira/1QUU.fasta.txt /homes/oliveira/fragments/1QUU.flib 10000 1000 2.5 linear
