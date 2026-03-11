name=$1

rm *.bp -rf

./example_pcms $name 0 &> $name\0.log &
./example_pcms $name 1 &> $name\1.log &

wait

echo "=====rdv0====="
cat $name\0.log

echo "=====rdv1====="
cat $name\1.log
