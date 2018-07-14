# WFCP

make
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
#
c="data_03"
#
mkdir runs/${dir}
path=$(pwd)
cd ${path}'/src'
./wfcp -ft ${path}'/data/'${c}'.turb' -fc ${path}'/data/'${c}'.cbl' -C 10 -time_limit 50 -rins 5 -relax 1 -polishing_time 150 -seed 9 -gap 0.1 -CC 3 -names 1
mv plot '../runs/'${dir}
mkdir plot
