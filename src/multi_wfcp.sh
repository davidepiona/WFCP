# WFCP
make
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir runs/${dir}
path=$(pwd)
cd ${path}'/data'
for t in `find . -type d -name 'w*' | cut -c 3-`
do
	echo ${t}	
	cd ${path}'/data/'${t}
	pwd
	for c in `find . -type f -name '*.cbl' | cut -c 3-`
	do
		echo "---"$c
		cd ${path}'/src'
		./wfcp -ft ${path}'/data/'${t}'/'${t}'.turb' -fc ${path}'/data/'${t}'/'${c} -C 10 -time_limit 250 -rins 10 -relax 1 -polishing_time 150 > ../runs/${dir}/run_${c}.log
	done
done
cd ${path}'/runs/'${dir}
grep "STAT" *.log |sort > results.csv