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
typedef unsigned __int128 readseq;
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

const readseq big_one = 1;
const readseq big_two = 2;
const readseq big_three = 3;

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
				ans = (ans>>2)+(big_one<<(ksize*2 - 2));
                ++i;
				break;
			case 'G': case 'g': // store as 10
				ans = (ans>>2)+(big_two<<(ksize*2 - 2));
                ++i;
				break;
			case 'T': case 't': // store as 11
				ans = (ans>>2)+(big_three<<(ksize*2 - 2));
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
        assert((ans >> (ksize*2)) ==0);
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
		switch ((e & (big_three<<(ksize*2 - 2)))>>(ksize*2 - 2)) {
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
    /*
    unsigned __int128 a = 0;
    unsigned __int128 b = 6 ;
    unsigned __int128 c = b << 96;

    if (a == b) {
        std::cerr << "A = B" << std::endl;
    } else {
        std::cerr << "NOPENOPE A = B" << std::endl;
    }
    if (a == c) {
        std::cerr << "A = B" << std::endl;
    } else {
        std::cerr << "NOPENOPE A = B" << std::endl;
    }


    std::cerr << sizeof(short) << std::endl;
    std::cerr << sizeof(int) << std::endl;
    std::cerr << sizeof(long) << std::endl;
    std::cerr << sizeof(long long) << std::endl;
    std::cerr << sizeof(unsigned __int128) << std::endl;
    */
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
    unsigned int totalloci = database.size();
    std::cout << "\tTotal loci:\t" << totalloci << std::endl;

    std::cerr << "Sorting database" << std::endl;
	//__gnu_parallel::sort(database.begin(),database.end());
    std::sort(database.begin(),database.end());

    std::cerr << "Deduplicating" << std::endl;
    // We no longer create another copy of database called dict, but do everything in place to save RAM
    //std::vector<readseq> dict; // unique'd version of database
    //dict.reserve(database.size());
    std::vector<int> exact_matches; // # of copies including itself
    exact_matches.resize(database.size(),1);
    long i = 0;
    for (unsigned int j=0; j < database.size(); ++j) {
        if (database[i] == database[j]) {
            ++(exact_matches[i]);
        } else {
            database[++i] = database[j];
        }
    }
    database.resize(i+1);
    exact_matches.resize(database.size());
    std::cerr << "\tUnique k-mers:\t" << database.size() << std::endl;
    long int exact_match_count = 0;
    for (unsigned long i = 0; i != database.size(); ++i ) {
        if (exact_matches[i]>1) {
            ++exact_match_count;
        }
    }
    std::cerr << "\tNon-unique loci:\t" << exact_match_count << std::endl;
    std::cerr << "\tUnique loci:\t" << database.size() - exact_match_count << std::endl;

    std::cerr << "Notating all locations with Hamming neighbors: " << std::endl;
    long database_size = database.size();
    std::vector<char> hamming_matches; // # of copies including itself
    hamming_matches.resize(database.size(),1);
#pragma omp parallel for
    for (long k = 0; k <= database_size; ++k ) {
        if (has_hamming_neighbor(database, database[k])==2) {
            hamming_matches[k] = 2;
        }
        /*if (k%1000 == 0) {
            std::cerr << k << "\r";
        }*/
    }
    std::cerr << std::endl;
    long int hamming_match_count = 0;
    long int unambig_count = 0;
    for (unsigned long i = 0; i != database.size(); ++i ) {
        if (hamming_matches[i]==2) {
            ++hamming_match_count;
        } else if (exact_matches[i]>1) {
            // deliberately left empty
        } else {
            ++unambig_count;
        }
    }
    std::cerr << "\tHamming-ed locations (incl. locations that already have exact matches): " << hamming_match_count << std::endl;
    //std::cerr << "\tUnambiguous locations: " << database.size() - exact_match_count - hamming_match_count << std::endl;
    std::cerr << "\tUnambiguous locations: " << unambig_count << std::endl;


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
