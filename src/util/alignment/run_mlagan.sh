#echo $1
fn=`python split.py $1`
#echo $fn
python run_mlagan_once.py $1 > ./alignment/$fn

