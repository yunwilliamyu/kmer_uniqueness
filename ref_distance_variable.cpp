/* Copyright (c) 2015 Yun William Yu */
/* ref_distance_variable.cpp */

#include <stdio.h>
#include <string>
#include <cstdio>
#include <vector>
#include <iostream>
#include <parallel/algorithm>
#include <algorithm>
#include <omp.h>
#include <cassert>

#ifdef GCCINT128
typedef __uint128_t readseq;
const int ksize_max = 64;
#else
typedef uint64_t readseq;
const int ksize_max = 32;
#endif
/*
 * A = 00
 * C = 01
 * G = 10
 * T = 11
*/ 
int ksize = ksize_max; // default k-mer size, but *always* gets overwritten

std::vector<readseq> encode_read_vector_full(FILE *input)
{
	std::vector<readseq> anslist;
	anslist.reserve(4ul*1024*1024*1024);
	readseq ans=0;
	int i=0;
    int ch;
    ch = getc(input);
	while ( ch != EOF )
	{
		switch (ch) {
			case 'A': case 'a': // store as 00
				ans >>= 2;
                ++i;
				break;
			case 'C': case 'c': // store as 01
				ans = (ans>>2)+(1ul<<(ksize*2 - 2));
                ++i;
				break;
			case 'G': case 'g': // store as 10
				ans = (ans>>2)+(2ul<<(ksize*2 - 2));
                ++i;
				break;
			case 'T': case 't': // store as 11
				ans = (ans>>2)+(3ul<<(ksize*2 - 2));
                ++i;
				break;
			case '\n':
                i = ( i < (ksize-1) ) ? i : (ksize-1) ;
                break;
            // if improper input, restart count
            case 'N':
			default:
                i = 0;
                ans = 0;
		}
        assert((ans >> (ksize*2-2)) ==0);
		if (i>=ksize)
			anslist.push_back(ans);
        ch = getc(input);
	}
	return anslist;
}

std::string decode_read(readseq e)
{
	char a[ksize+1];
	a[ksize]=0;
	for (int i=31; i>=0; --i)
	{
		switch ((e & (3ul<<(ksize*2 - 2)))>>(ksize*2 - 2)) {
			case 0:
				a[i]='A';
				break;
			case 1:
				a[i]='C';
				break;
			case 2:
				a[i]='G';
				break;
			case 3:
				a[i]='T';
				break;
		}
		e <<= 2;
	}
	std::string ans = a;
	return ans;
}

// All the potential Hamming neighbors
inline readseq mask(int i) {return ~(3UL << (2*i));}
inline readseq add(int i,readseq j) {return j << (2*i); }
inline readseq testers(int i, int j, readseq work) {return (work & mask(i)) + add(i,j);}

// Returns 2 if has a hamming neighbor in dict
int has_hamming_neighbor(const std::vector<readseq> & dict, const readseq x) {
    readseq work = x;
    readseq test;
    for (int i=0; i<ksize; ++i) {
        for (int j=0; j<4; ++j) {
            test = testers(i,j,work);
            if((test!=work)&& \
                std::binary_search(dict.begin(),dict.end(),test)) {
                return 2;
                break;
            }
        }
    }
    return 0;
}

int main( int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " kmer_size input_file" << std::endl;
		std::cerr << "\tExamines the uniqueness of the reference with respect to exact duplicates and Hamming neighbors." << std::endl;
		std::cerr << "\tkmer_size should be between 0 and " << (ksize_max) << "." << std::endl;
		exit(0);
    }

    ksize = atoi(argv[1]);
    // Must fit in integer size.
    if ((ksize <= 0)||(ksize > ksize_max)) {
		std::cerr << "\tkmer_size should be between 0 and " << (ksize_max) << "." << std::endl;
        exit(0);
    }

    FILE *input = fopen(argv[2], "r");

    std::vector<readseq> database;
    std::cerr << "Reading reference to database" << std::endl;
    database = encode_read_vector_full(input);
    std::cout << "\tTotal loci:\t" << database.size() << std::endl;

    std::cerr << "Sorting database" << std::endl;
	__gnu_parallel::sort(database.begin(),database.end());

    std::cerr << "Deduplicating" << std::endl;
    std::vector<readseq> dict; // unique'd version of database
    dict.reserve(database.size());
    std::vector<char> key; // 0 = unique; 1 = exact match; 2 = no exact match but hamming neighbor
    key.resize(database.size());
    long i = 0;
    readseq last = -1;
    for (auto it = database.begin(); it != database.end(); ++it) {
        if (last != *it) {
            dict.push_back(*it);
            ++i;
        } else {
            key[i-1] = 1;
        }
        last = *it;
    }
    // Clear database
    std::vector<readseq>().swap(database);
    std::cerr << "\tUnique k-mers:\t" << dict.size() << std::endl;
    long int exact_match_count = 0;
    for (unsigned long i = 0; i != dict.size(); ++i ) {
        if (key[i]==1) {
            ++exact_match_count;
        }
    }
    std::cerr << "\tNon-unique loci:\t" << exact_match_count << std::endl;
    std::cerr << "\tUnique loci:\t" << dict.size() - exact_match_count << std::endl;

    std::cerr << "Notating all locations with Hamming neighbors: " << std::endl;
    long dict_size = dict.size();
#pragma omp parallel for
    for (long k = 0; k <= dict_size; ++k ) {
        if ((key[k]==0)&&(has_hamming_neighbor(dict, dict[k])==2)) {
            key[k] = 2;
        }
        /*if (k%1000 == 0) {
            std::cerr << k << "\r";
        }*/
    }
    std::cerr << std::endl;
    long int hamming_match_count = 0;
    for (unsigned long i = 0; i != dict.size(); ++i ) {
        if (key[i]==2) {
            ++hamming_match_count;
        }
    }
    std::cerr << "\tHamming-ed locations: " << hamming_match_count << std::endl;
    std::cerr << "\tUnambiguous locations: " << dict.size() - exact_match_count - hamming_match_count << std::endl;


/*
    for (auto it = database.begin(); it != database.end(); ++it) {
        std::cout << decode_read(*it) << std::endl;
    }
    std::cout << "Pause" << std::endl;
    for (unsigned long i = 0; i != dict.size(); ++i) {
        std::cout << decode_read(dict[i]) << ": " << int(key[i]) << std::endl;
    }
  */  

}
