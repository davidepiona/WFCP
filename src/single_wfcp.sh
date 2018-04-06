# WFCP

make
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
#
t="wf01"
c="wf01_cb01_capex.cbl"
#
mkdir runs/${dir}
path=$(pwd)
cd ${path}'/src'
./wfcp -ft ${path}'/data/'${t}'/'${t}'.turb' -fc ${path}'/data/'${t}'/'${c} -C 10 -time_limit 250 -rins 10 -relax 1 -polishing_time 150 > ../runs/${dir}/run_${c}.log

