# WFCP
make
rm -r plot
mkdir plot
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir runs/${dir}
path=$(pwd)
cd ${path}'/data'
for c in `find . -type f -name '*.turb' | cut -c 3-9 |sort`
do
	echo "---"$c
	cd ${path}'/src'
	./wfcp -ft ${path}'/data/'${c}'.turb' -fc ${path}'/data/'${c}'.cbl' -C 10 -rins 5 -relax 1 -gap 0.1 -seed 9 -CC 3 -time_limit 30 -names 1 > ../runs/${dir}/run_${c}.log
done
#mv *.png '../runs/'${dir}
mv plot '../runs/'${dir}
mkdir plot
cd ${path}'/runs/'${dir}
grep "STAT" *.log |sort > results.csv
