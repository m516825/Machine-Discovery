#!/bin/bash
BIGRAM_FILE=./bigram.txt
ENCODE_FILE=./encode.txt
TEST_FILE=./test.txt
OUTPUT_FILE=./pred.txt

python hw1-1.py --bigram_file $BIGRAM_FILE --encode_file $ENCODE_FILE --output_file $OUTPUT_FILE --test_file $TEST_FILE