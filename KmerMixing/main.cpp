#include <zlib.h>
#include <cstring>
#include <string>
#include <iostream>
#include <regex>
#include <omp.h>
#include <fstream>
#include <bitset> 
#include "MurmurHash3.h"
#include "robin_hood.h"
#include "kseq.h"
#include "minimizer.h"

KSEQ_INIT(gzFile, gzread);  
uint32_t K=31;
uint32_t seed = 101;
uint64_t MIN_HASH_VAL = 0x07FFFFFFFFFFFFFF;
uint64_t mask = 0x3FFFFFFFFFFFFFFF;

// trim from start (in place)
inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

inline void trim_str(std::string& str)
{
    rtrim(str);
    ltrim(str);
}
std::vector<std::string> split_into_fields(std::string const & s, std::string const & tokens) {
  char *c_s = const_cast<char*>(s.c_str());
  char *c_tokens = const_cast<char*>(tokens.c_str());
  char *c_split = strtok(c_s, c_tokens); // split string into delimited cstrings
  std::vector<std::string> fields;
  while (c_split != NULL) {
    fields.push_back(c_split);
    c_split = strtok(NULL, c_tokens);
  }
  return fields;
}

bool passes_hash(mzr::Kmer kmer) {
    uint64_t * dat = &kmer.seq;
    uint64_t ho[2];
    MurmurHash3_x64_128(dat, 8, seed, ho);
    return ho[0] <= MIN_HASH_VAL;
}

void build_kmer_set_from_fq(std::string const & r1, std::string const & r2, uint32_t k, robin_hood::unordered_set<uint64_t> & rset) {


    gzFile r1_fin = gzopen(r1.c_str(), "r");
    gzFile r2_fin = gzopen(r2.c_str(), "r");
    kseq_t *fq_r1 = kseq_init(r1_fin);
    kseq_t *fq_r2 = kseq_init(r2_fin);
    uint64_t hashed = 0, total = 0;
    while (kseq_read(fq_r1) >= 0) {
        mzr::KmerGenerator kmer_gen(fq_r1->seq.s, k, mask, true);
        while (!kmer_gen.empty()) {
            mzr::Kmer kmer = kmer_gen.get_kmer();
            if (passes_hash(kmer)) {
                rset.insert(kmer.seq);
                hashed++;
            }
            total++;
        }
    }
    while (kseq_read(fq_r2) >= 0) {
        mzr::KmerGenerator kmer_gen(fq_r2->seq.s, k, mask, true);
        while (!kmer_gen.empty()) {
            mzr::Kmer kmer = kmer_gen.get_kmer();
            if (passes_hash(kmer)) {
                rset.insert(kmer.seq);
                hashed++;
            }
            total++;
        }
    }

    // clean up
    std::cout << "For sample " << r1 << std::endl;
    std::cout << "Fraction of kmers hashed: " << hashed / (double) total << std::endl;
    std::cout << "Hashed : total: " << hashed  << " : " <<  total << std::endl;
    kseq_destroy(fq_r1);
    kseq_destroy(fq_r2);
    gzclose(r1_fin);
    gzclose(r2_fin);
}

void add_sample_to_counter(robin_hood::unordered_map<uint64_t, uint8_t> & counter, robin_hood::unordered_set<uint64_t> const & sample_set) {
    for (uint64_t kmer : sample_set) {
        counter[kmer]++;
    }
}

robin_hood::unordered_map<uint64_t, uint8_t> get_kmer_frequency_by_sample(std::vector<std::string> const & sample_list, uint32_t k, int threads) {
    std::vector<robin_hood::unordered_set<uint64_t>> sets(threads, robin_hood::unordered_set<uint64_t>());
    robin_hood::unordered_map<uint64_t, uint8_t> freq;
    std::cout << "extracting kmers from reads using " << threads << " threads..." << std::endl;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (auto s: sample_list) {
            build_kmer_set_from_fq(s + "_R1.fastq.gz", s + "_R2.fastq.gz", k, sets[omp_get_thread_num()]);
        }
    }

    std::cout << "counting sample occurence of read kmers..." << std::endl;
    for (auto sample_set : sets) {
        add_sample_to_counter(freq, sample_set);
    }
    return freq;
}

robin_hood::unordered_map<uint64_t, uint16_t> get_kmer_frequency_by_node(std::string const & node_fasta, uint32_t k, bool is_copan) {

    robin_hood::unordered_map<uint64_t, uint16_t> freq;
    robin_hood::unordered_set<uint64_t> node_set;
    
    std::string current_node_name = "";
    std::ifstream contig_file(node_fasta);
    std::string line;
    std::cout << node_fasta << std::endl;
    uint64_t hashed = 0, total = 0;
    while (contig_file){
        std::getline(contig_file, line);
        trim_str(line);
        if (line.size() > 0 && line[0] == '>')
        {

            // get next element
            std::string hdr = line;
            std::string seq;
            while(contig_file.peek() != '>' && !contig_file.eof())                                      
            {
                std::getline(contig_file, line);
                trim_str(line);
                seq.append(line);
            }


            std::string next_node_name;
            if (is_copan) {
                auto fields = split_into_fields(hdr, ":");
                next_node_name = fields[0];
            }
            else {
                std::regex re("(>NODE_[0-9]+)_");
                std::smatch m;
                std::regex_search(hdr, m, re);
                next_node_name = m[1].str();
            }

            // we've reached a new node...
            if (next_node_name != current_node_name) {
                // add kmers to the counter
                for (uint64_t kmer : node_set) {
                    freq[kmer]++;
                }
                // reset the node_set and name
                node_set.clear();
                current_node_name = next_node_name;
            }

            // add sequence to node set
            mzr::KmerGenerator kmer_gen(seq, k, mask, true);
            while (!kmer_gen.empty()) {
                mzr::Kmer kmer = kmer_gen.get_kmer();
                if (passes_hash(kmer)) {
                    hashed++; 
                    node_set.insert(kmer.seq);
                }
                total++;
            }
        }
    }

    // add the last node's data
    for (uint64_t kmer : node_set) {
        freq[kmer]++;
    }
    std::cout << "Fraction of kmers hashed: " << hashed / (double) total << " for tool " << node_fasta << std::endl;
    return freq;
}

std::vector<std::string> get_sample_list(std::string sample_list) {
    std::vector<std::string> slist;
    std::ifstream fin(sample_list);
    std::string line;
    while (fin) {
        std::getline(fin, line);
        trim_str(line);
        if (line.empty()) break;
        slist.emplace_back(line);
    }
    return slist;
}

int main(int argc, char* argv[]) { 

    if (argc != 5) {
        std::cout << "Usage: <exe> <list_of_fq_prefixes> <datadir> <outfile_prefix> <max_threads>" << std::endl;
        std::cout << argc << std::endl;
        exit(-1);
    }

    std::string sample_list = argv[1];
    std::string datadir = argv[2];
    std::string out_file = argv[3];
    int threads = std::stoi(argv[4]);
    std::vector<std::string> slist = get_sample_list(sample_list);
    if (slist.size() < threads) {
        threads = slist.size();
    }
    omp_set_num_threads(threads);


    std::string cpg_ms_1000_dt_02 = datadir + "coasm_50_rep_2_ms_1000_dt_0.02.fasta";
    std::string cpg_ms_1000_dt_05 = datadir + "coasm_50_rep_2_ms_1000_dt_0.05.fasta";
    std::string cpg_ms_1000_dt_1 = datadir + "coasm_50_rep_2_ms_1000_dt_0.1.fasta";
    std::string cpg_ms_100_dt_02 = datadir + "coasm_50_rep_2_ms_100.fasta";
    std::string k141 = datadir + "coasm_50_rep_2/k141.fastg";
    std::string k59 =  datadir + "coasm_50_rep_2/k59.fastg";
    std::cout << cpg_ms_1000_dt_02 << std::endl;
    std::cout << cpg_ms_1000_dt_05 << std::endl;
    std::cout << cpg_ms_100_dt_02 << std::endl;
    std::cout << k141 << std::endl;
    std::cout << k59 << std::endl;

    auto kmer_to_sample_count = get_kmer_frequency_by_sample(slist, K, threads);
    std::cout << "total kmers across all reads: " << kmer_to_sample_count.size() << std::endl;

    #pragma omp parallel sections
    {
	#pragma omp section
	{
    		auto A = get_kmer_frequency_by_node(cpg_ms_1000_dt_02, K, true);
    		std::cout << "total kmers in copangraph: " << A.size() << std::endl;
    		std::ofstream fout_A(out_file + "copan_ms_1000_dt_02_kmix.csv");
    		fout_A << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : A) {
    		    fout_A << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
    		fout_A.close();
    		A.clear();
	}
		
        #pragma omp section
	{
    		auto B = get_kmer_frequency_by_node(cpg_ms_1000_dt_05, K, true);
    		std::cout << "total kmers in copangraph: " << B.size() << std::endl;
		std::ofstream fout_B(out_file + "copan_ms_1000_dt_05_kmix.csv");
    		fout_B << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : B) {
    		    fout_B << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
    		fout_B.close();
    		B.clear();
	}
	
	# pragma omp section
	{
    		auto C = get_kmer_frequency_by_node(cpg_ms_1000_dt_1, K, true);
    		std::cout << "total kmers in copangraph: " << C.size() << std::endl;
		std::ofstream fout_C(out_file + "copan_ms_1000_dt_1_kmix.csv");
    		fout_C << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : C) {
    		    fout_C << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
    		fout_C.close();
    		C.clear();
	}	
	
	# pragma omp section
    	{
    		auto D = get_kmer_frequency_by_node(cpg_ms_100_dt_02, K, true);
    		std::cout << "total kmers in copangraph: " << D.size() << std::endl;
		std::ofstream fout_D(out_file + "copan_ms_100_dt_02_kmix.csv");
    		fout_D << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : D) {
    		    fout_D << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
    		fout_D.close();
    		D.clear();
	}

	#pragma omp section
	{
    		auto E= get_kmer_frequency_by_node(k141, K, false);
    		std::cout << "total kmers in coasm: " << E.size() << std::endl;
		std::ofstream fout_E(out_file + "k141_kmix.csv");
    		fout_E << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : E) {
    		    fout_E << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
    		fout_E.close();
    		E.clear();
	}
	
	#pragma omp section
	{
    		auto F = get_kmer_frequency_by_node(k59, K, false);
    		std::cout << "total kmers in coasm: " << F.size() << std::endl;
		std::ofstream fout_F(out_file + "k59_kmix.csv");
    		fout_F << "node_count,sample_count" << std::endl;
    		for (auto const & kmer : F) {
    		    fout_F << kmer.second << ',' << unsigned(kmer_to_sample_count[kmer.first]) << '\n';
    		}
		fout_F.close(); 
		F.clear();
	}
    }

}
