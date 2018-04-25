# WFCP
make
rm -r plot
mkdir plot
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir runs/${dir}
path=$(pwd)
settings="-C 10 -rins 5 -relax 1 -seed 9 -CC 3 -time_limit 30 -polishing_time 260 -names 1"
cd ${path}'/data'
for c in `find . -type f -name '*.turb' | cut -c 3-9 |sort`
do
	echo "---"$c
	cd ${path}'/src'
	./wfcp -ft ${path}'/data/'${c}'.turb' -fc ${path}'/data/'${c}'.cbl' ${settings} > ../runs/${dir}/run_${c}.log
done
mv *.png '../runs/'${dir}
mv plot '../runs/'${dir}
mkdir plot
cd ${path}'/runs/'${dir}
echo ${settings} > settings.txt
grep "STAT" *.log |sort > results.csv
