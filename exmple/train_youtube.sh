set -x

# download youtube-links dataset
zipfile="youtube-links.txt.gz"
if test -e "$zipfile";
then
    echo 'zip file exists'
else
    wget http://socialnetworks.mpi-sws.mpg.de/data/youtube-links.txt.gz
fi

# generate the network
zcat youtube-links.txt.gz | awk -F '	' '{print $1" "$2" 1"}' > net.txt

# run the comment
../cli/cli -model HPE -train net.txt -save rep.txt -dimensions 64 -sample_times 5 -walk_steps 5 -negative_samples 5 -alpha 0.025 -threads 1
