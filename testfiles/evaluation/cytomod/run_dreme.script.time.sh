for file in ../files/SIM_seqs/cytomod_format/FILENAME; do newname=$(basename $file).dreme;echo $file $newname; mytime="$(time ( ./run_dreme.sh $file) 2>&1 1>/dev/null )"; echo $mytime > $newname.log;done




