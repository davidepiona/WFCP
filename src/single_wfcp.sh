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
./wfcp -ft ${path}'/data/'${c}'.turb' -fc ${path}'/data/'${c}'.cbl' -C 3 -time_limit 150 -rins 10 -relax 1 -polishing_time 150 -CC 1

