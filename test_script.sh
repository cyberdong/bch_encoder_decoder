#!/bin/bash
# BCH Encoder/Decoder

# gcc data_generator.c -o data.exe
# gcc bch_encoder.c -o bch_encoder.exe
# gcc error.c -o error.exe
# gcc bch_decoder.c -o bch_decoder.exe

make

# MM=13
# KK=4096
# TT=8
# PP=8
# EE=16

# MM -- BCH code over GF(2**mm), where mm is Galois fields
# 以下参数适用于A15的NAND控制器, Galois fields=15
MM=15
echo "MM="$MM

KK_bytes=2048
let KK_bits=$KK_bytes*8
echo "KK_bytes="$KK_bytes
echo "KK_bits="$KK_bits

TT=16
echo "TT="$TT
#TT -- Number of errors that can be corrected
#TT is 16 bits per 2048 Bytes, 会产生30 Bytes的BCH ECC;可纠正2048字节中的16bit的错误; 注入超过16位的错误，就会有不能成功纠错的位

# PP -- Number of substreams to calculate in parallel
PP=8

# EE -- the number of error bits injected
EE=16
echo "EE="$EE

./data_gen -n $KK_bytes > data_in.txt
./bch_encoder -m $MM -k $KK_bits -t $TT -p $PP < data_in.txt > data_codeword.txt
./error -e $EE < data_codeword.txt > data_error.txt
./bch_decoder -m $MM -k $KK_bits -t $TT -p $PP < data_error.txt > data_out.txt



