# WFCP
cd "/home/davide/Scrivania/WFCP/data/wf01"
for F in `find . -type f -name '*.cbl' | cut -c 3-`
do
echo ${F}
cd "/home/davide/Scrivania/WFCP/src"
./wfcp -ft '/home/davide/Scrivania/WFCP/data/wf01/wf01.turb' -fc '/home/davide/Scrivania/WFCP/data/wf01/'$F -C 10 -time_limit 100 -rins 10 -relax 1 > run_wf01_$F.log
done

grep "STAT" *.log |sort > results.csv
