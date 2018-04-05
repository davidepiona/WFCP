# WFCP
cd "/home/davide/Scrivania/WFCP/data"
for t in `find . -type d -name 'wf02' | cut -c 3-`
do
	echo ${t}	
	cd "/home/davide/Scrivania/WFCP/data/"${t}
	pwd
	for c in `find . -type f -name '*.cbl' | cut -c 3-`
	do
		echo "---"$c
		cd "/home/davide/Scrivania/WFCP/src"
		./wfcp -ft '/home/davide/Scrivania/WFCP/data/'${t}'/'${t}'.turb' -fc '/home/davide/Scrivania/WFCP/data/'${t}'/'${c} -C 10 -time_limit 15 -rins 10 -relax 1 > ../runs/run_${c}.log
	done
done
cd "/home/davide/Scrivania/WFCP/runs"
grep "STAT" *.log |sort > results.csv
