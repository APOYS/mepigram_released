file=$1

mytime="$(time ( ./run_dreme.sh $file) 2>&1 1>/dev/null )"
echo "$mytime" 
