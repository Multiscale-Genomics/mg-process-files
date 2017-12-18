

pytest -s tests/test_bed_functions.py
tc=$?
rc=$(($rc + $tc))
./tidy_data.sh

pytest -s tests/test_gff3_functions.py
tc=$?
rc=$(($rc + $tc))
./tidy_data.sh

pytest -s tests/test_wig_functions.py
tc=$?
rc=$(($rc + $tc))
./tidy_data.sh

pytest -s tests/test_json3d_functions.py
tc=$?
rc=$(($rc + $tc))
./tidy_data.sh

if [[ $rc != 0 ]]; then exit $rc; fi