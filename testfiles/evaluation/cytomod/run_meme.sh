file=$1
newname=$(basename $file).meme
echo $file
echo $newname


#echo dreme -v 1 -oc $newname -alph alphabet.meme.txt -p $file -t 18000 -e 0.05
#dreme -v 1 -oc $newname -alph alphabet.meme.txt -p $file -t 18000 -e 0.05

echo meme $file -alph alphabet.meme.txt -oc $newname -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -revcomp

meme $file -alph alphabet.meme.txt -oc $newname -nostatus -time 18000 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -revcomp
