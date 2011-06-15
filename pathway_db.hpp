#ifndef _PATHWAY_DB_H_
#define _PATHWAY_DB_H_
#include "aligator.hpp"
using namespace boost;
namespace po = boost::program_options;
#include <boost/random/variate_generator.hpp>
#include <boost/random/geometric_distribution.hpp>
#include <boost/random.hpp>

struct my_pathway_data{
    double catepvalue;
    double pvalue;
    string name;

    int id;
};

bool sort_pred(const my_pathway_data& left, const my_pathway_data& right);


class gene_info{
public:
    int chr_id;
    int start;
    int end;
    gene_info(int chr_id,int start,int end);
    gene_info();
};
bool sort_gene_info(const gene_info &left, const gene_info &right);
class pathway_db{
    vector<string> snp_ids_;
    vector<double> snp_pvalues_; // corresponding pvalues of SNPs_
    map<string,int> gene_name2id_;
    vector<vector<int> > snp2genes_;
    vector<vector<int> > gene2pathways_;
    vector<vector<int> > pathway2genes_;
    vector<int > removed_pathway_;
    vector<vector<double> > pathway_correlation_;
    vector<vector<int> > snp2pathways_;
    vector<gene_info> geneinfo_;
    map<string,int> snp_idx_;
    vector<string> pathway_names_;
    map<string,int> pathway_name2id_;
    double sig_pvalue_;
    double sig_gene_fraction_;
    string ofile;
    string lfile;
    FILE * lfd;

    int max_genes_;
    int min_genes_;
    double max_correlation_;
    vector<int> total_num_gene_pathway;
    int permutation_number;
    int simulate_permutation_number;
    double ld_threshold_;
    vector<int> snp_contributed_pathway_;

    int seed_;
public:
    pathway_db(int argc, char *argv[]);
    int study();

    int permutation_prepare(int permutation_number, int N, vector<vector<unsigned short> > *sim_num_gene_pathway);
    int permutation_study(vector<vector<unsigned short> > *sim_matrix,int N,vector<unsigned short> * num_gene_pathway, vector<double> *catepvalue_pathway);
    int calc_correlation_pathway();
    int remove_correlated_pathway( vector<double> &catepavlue);
    int count_significant_pathway(vector<unsigned short> *sgenes, vector<unsigned short> *spathways);
    int count_significant_pathway_by_snp(vector<int> *ssnps, vector<unsigned short> *spathways);
    int count_significant_genes_by_ld(vector<int> *genelist);
    int foreach_count_significant_genes_by_ld(vector<unsigned short> *genenum, vector<vector<int> > *genelists);
};

#endif /* _PATHWAY_DB_H_ */


