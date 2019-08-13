# objects = data_generator.o bch_encoder.o error.o bch_decoder.o

CC = gcc

all: data bch_encoder error bch_decoder

data: data_generator.o
	$(CC) -o data_gen data_generator.o -lm

bch_encoder: bch_encoder.o
	$(CC) -o bch_encoder bch_encoder.o -lm

error: error.o
	$(CC) -o error error.o -lm

bch_decoder: bch_decoder.o
	$(CC) -o bch_decoder bch_decoder.o -lm

.PHONY : clean
clean :
	-rm -f data_gen bch_encoder error bch_decoder *.o

