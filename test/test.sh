#!/bin/bash

set -e

(
    cd ..
    python rhoterm_caller.py --untreated test/untreated.gff --treated test/treated.gff
    diff Test_output_summary.txt test/expected_summary.txt
    diff Test_output_peaks.txt test/expected_peaks.txt
)
