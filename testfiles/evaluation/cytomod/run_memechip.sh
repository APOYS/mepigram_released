file=$1
newname=$(basename $file)
echo $file
echo $newname
time meme-chip -spamo-skip -fimo-skip -oc $newname.memechip -time 300 -xalph alphabet.meme.txt -order 1 -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 3 -dreme-e 0.05 $file
