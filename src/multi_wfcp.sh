# WFCP
make
rm -r plot
mkdir plot
cd ..
dir=$(date +"%Y-%m-%d_%H-%M-%S")
mkdir runs/${dir}
path=$(pwd)
settings="-model 0 -CC 2 -time_limit 600 -time_loop 120 -relax 3 -names 1"
cd ${path}'/data'
count=0
for c in `find . -type f -name '*.turb' | cut -c 3-9 |sort`
do
	echo "---"$c
	if [ $count -le 5 ]; then
		cSub="-C 10"
	fi
	if [ $count -gt 5 ] && [ $count -le 13 ]; then
		cSub="-C 100"
	fi
	if [ $count -gt 13 ] && [ $count -le 17 ]; then
		cSub="-C 4"
	fi
	if [ $count -gt 17 ]; then
		cSub="-C 10"	
	fi
	cd ${path}'/src'
	./wfcp -ft ${path}'/data/'${c}'.turb' -fc ${path}'/data/'${c}'.cbl' ${cSub} ${settings} > ../runs/${dir}/run_${c}.log
done
mv *.png '../runs/'${dir}
mv plot '../runs/'${dir}
mkdir plot
cd ${path}'/runs/'${dir}
echo ${settings} > settings.txt
grep "STAT" *.log |sort > results.csv
