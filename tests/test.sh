#!/usr/bin/env bash
#
# test.sh
#
# Tests for find_differential_primers.py

echo "Testing find_differential_primers.py"

# Test 1: simple, default operation
test1="find_differential_primers.py -i test.conf -l test1.log -o test1"

# Test 2: generate hybridisation probe oligos
test2="find_differential_primers.py -i test.conf -l test2.log -o test2 \
    --hybridprobe --numreturn 5"

# Test 3: filter on GC content at 3' end
test3="find_differential_primers.py -i test.conf -l test3.log -o test3 \
    --filtergc3prime"

# Test 4: screen against E.coli
test4="find_differential_primers.py -i test.conf -l test4.log -o test4 \
    --blastdb=sequences/e_coli_screen.fna"

cmdarray=( "${test1}" "${test2}" "${test3}" "${test4}" )
for t in "${cmdarray[@]}"
do
    echo ${t}
    eval time ${t}
done


