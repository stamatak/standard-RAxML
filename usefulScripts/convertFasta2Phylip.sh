#! /bin/sh

if [ $# != 1 ]; then
    echo "USAGE: ./script <fasta-file>"
    exit
fi

numSpec=$(grep -c  ">" $1)
tmp=$(cat $1 | sed "s/>[ ]*\(\w*\).*/;\1</"  | tr -d "\n" | tr -d ' '  | sed 's/^;//' | tr "<" " " )
length=$(($(echo $tmp | sed 's/[^ ]* \([^;]*\);.*/\1/'   | wc -m ) - 1))

echo "$numSpec $length"
echo  $tmp | tr ";" "\n"
