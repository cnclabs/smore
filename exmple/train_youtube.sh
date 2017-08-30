set -x
Threads=1

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
../cli/deepwalk -train net.txt -save rep_dw.txt -undirected 1 -dimensions 64 -walk_times 1 -walk_steps 40 -window_size 5 -negative_samples 5 -alpha 0.025 -threads $Threads
../cli/walklets -train net.txt -save rep_wl.txt -undirected 1 -dimensions 64 -walk_times 1 -walk_steps 40 -window_min 2 -window_max 5 -negative_samples 5 -alpha 0.025 -threads $Threads
../cli/line -train net.txt -save rep_line1.txt -undirected 1 -order 1 -dimensions 64 -sample_times 10 -negative_samples 5 -alpha 0.025 -threads $Threads
../cli/line -train net.txt -save rep_line2.txt -undirected 1 -order 2 -dimensions 64 -sample_times 10 -negative_samples 5 -alpha 0.025 -threads $Threads
../cli/hpe -train net.txt -save rep_hpe.txt -undirected 1 -dimensions 64 -sample_times 10 -walk_steps 5 -negative_samples 5 -alpha 0.025 -threads $Threads
