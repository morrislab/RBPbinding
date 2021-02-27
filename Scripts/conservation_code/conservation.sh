list=$(ls ${1}*${2})

len=$(expr length $1)

for f in $list
do
echo $f
echo ${f::-6}.cons
sen=${f%${2}}
seq=${sen:${len}}
python score_conservation.py -s js_divergence -w 3 -g 0.6 -a $seq -o ${f::-6}.cons $f 
#python score_conservation.py -s js_divergence -w 3 -g 0.6 -o ${f::-6}.cons $f
done




