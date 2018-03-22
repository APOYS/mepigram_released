file=$1
newname=$(basename $file).dreme
echo $file
echo $newname

#time dreme -v 5 -p $file -alph alphabet.meme.txt -eps -o $file.dreme.outi

#mytime=$"(time (dreme -v 1 -oc $newname -alph alphabet.meme.txt -p $file -t 18000 -e 0.05) 2>&1 1>/dev/null)"
echo dreme -v 1 -oc $newname -alph alphabet.meme.txt -p $file -t 18000 -e 0.05
dreme -v 1 -oc $newname -alph alphabet.meme.txt -p $file -t 18000 -e 0.05
#echo "$mytime" > $newname.log  
