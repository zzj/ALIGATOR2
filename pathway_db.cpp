#include "pathway_db.hpp"
#define MAX 1000000

typedef std::pair<double, int> my_pair;

bool sort_pred(const my_pathway_data& left, const my_pathway_data& right)
{
    if (abs(left.pvalue-right.pvalue)<0.0000001){
        if (abs(left.catepvalue-right.catepvalue)>0.0000001)
            return left.catepvalue < right.catepvalue;
        else return left.name<right.name;
    }
    return left.pvalue < right.pvalue;
}
bool sort_pair_by_first(const pair<double,int> &a, const pair<double,int> &b){
    return a.first<b.first;
}
bool sort_gene_info(const gene_info& left, const gene_info& right)
{
    if (left.chr_id!=right.chr_id){
        return left.chr_id<right.chr_id;
    }
    double lmean=(left.start+left.end)/2;
    double rmean=(right.start+right.end)/2;
    return lmean<rmean;
}

gene_info::gene_info(int chrid,int start,int end){
    this->chr_id=chrid;
    this->start=start;
    this->end=end;
}
gene_info::gene_info(){
    gene_info(0,0,0);
}
pathway_db::pathway_db(int argc, char * argv[]){
    string sfile, gfile, pfile;
    string cfile;
    string ginfofile;
    po::options_description desc("This is a C++ implementation of aligator. \n Allowed options are: ");
     
    desc.add_options()
        ("help", "produce help message")
        ("snps_file", po::value<string>(&sfile),

	   
         "The snp file with pvalue ")

        ("gene_file", po::value<string>(&gfile),

         "The snp file with corresponding genes ")

        ("gene_info_file", po::value<string>(&ginfofile),
         "The file for gene position information")
                    
        ("min_genes", po::value<int>(&min_genes_)->default_value(0),
         "The minimum number of genes in each pathway (excluded)")
          
        ("max_genes", po::value<int>(&max_genes_)->default_value(MAX),
         "The maximum number of genes in each pathway (excluded)")

        ("max_correlation", po::value<double>(&max_correlation_)->default_value(1),
         "The maximum correlation")
          
        ("ld_threshold", po::value<double>(&ld_threshold_),
         "The minimum gene distance in a pathway, the unit is Mbps, not bps.")
		  
        ("pathway_file", po::value<string>(&pfile),
         "The gene file with corresponding pathways. ")

        ("config", po::value<string>(&cfile)->default_value("NOCONFIG"),
         "Configuration file ")
	  
        ("log_file", po::value<string>(&lfile)->default_value("NOCONFIG"),
         "Log file ")
          
        ("output_file", po::value<string>(&ofile)->default_value("output1.txt"),
         "outputfile ")

        ("significant_pvalue", po::value<double>(&sig_pvalue_)->default_value(0.00197607),
         "The SNP will be counted as a significant SNP if it is <=significant_pvalue")

        ("significant_gene_fraction", po::value<double> (&sig_gene_fraction_)->default_value(1),
         "The significant gene fraction (an alternative to --significant_pvalue. ")
	  
        ("permutation_number", po::value<int>(&permutation_number)->default_value(10000),
         "The permutations needed to get categorical pvalue")
	  
        ("simulate_permutation_number", po::value<int> (&simulate_permutation_number)->default_value(1000),
         "The number of simulations for multiple tests correction"	       )

        ("seed", po::value<int>(&seed_)->default_value(time(NULL)))
        //int permutation_number=10000;
        //int simulate_permutation_number=1000;
     
        ;
    po::variables_map vm;
    
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (cfile!="NOCONFIG"){
        std::ifstream in(cfile.c_str() ); 
        po::store(po::parse_config_file(in, desc), vm);
        po::notify(vm);
    }
    if (vm.count("help")) {
        cout << desc << "\n";
        exit(0);
    }
    ld_threshold_=ld_threshold_*1000000;
    if (vm.count("gene_file")<1||vm.count("snps_file")<1 || vm.count("pathway_file")<1){
	  
        cout << "gene_file or snps_file or pathway_file is not specified\n";
        exit(0);
    }

    vector<string> geneline;
    vector<string> snp_list;
    vector<vector<int> > snp2genes;
    map<string,double> snp_set;
    set<string> snp_names;
    int i,j;
    printf("Loading snp file ...\n");

    //snp_db.set_seperator(' ');
    FILE *file_snp=fopen(sfile.c_str(),"r");
    if (file_snp==NULL) {
        fprintf(stderr, "Cannot open file %s!", sfile.c_str());
    }
	 
    double snp_pvalue;
    char snp_name[1000];
    snp_ids_.clear();
    snp_pvalues_.clear();
    while(fscanf(file_snp,"%s%lf", snp_name, &snp_pvalue)==2){
        snp_ids_.push_back(snp_name);
        snp_names.insert(snp_name);
        snp_pvalues_.push_back(snp_pvalue);
    }
     
    //snp_db.cvs_load(sfile.c_str());
    //snp_ids_=snp_db.cvs_col_str(0);
    //snp_pvalues_=snp_db.cvs_col_double(1);
    if (snp_ids_.size()!=snp_pvalues_.size()) {
        fprintf(stderr, "ERROR: error at %s:%d\n",
                __FILE__, __LINE__);
        exit(0);
    }
     
    printf("Loading gene file ...\n");
    FILE *fd=open_file((const char *)gfile.c_str(),"r");
    if (fd==NULL) {
        fprintf(stderr, "Cannot open file %s!", gfile.c_str());
    }
    char genename[100];
    char snpname[20];
    char pathwayname[20];
    int num_gene;
    while(fscanf(fd,"%s%d",snpname, &num_gene)==2){
        snp_list.push_back(snpname);
        if (snp_names.find(snpname)==snp_names.end()){
            for (i=0;i<num_gene;i++){
                fscanf(fd,"%s", genename);
            }
            vector<int> temp;
            snp2genes.push_back(temp);
            continue;
        }
        vector<int> temp;
        for (i=0;i<num_gene;i++){
            fscanf(fd,"%s", genename);
            if (gene_name2id_.find(genename)==gene_name2id_.end()){
                int id=gene_name2id_.size();
                gene_name2id_[genename]=id;
                //printf("%s\n",temp[j].c_str());
            }
            temp.push_back(gene_name2id_[genename]);
        }
        snp2genes.push_back(temp);
    }

    gene2pathways_.resize(gene_name2id_.size());
    printf("There are %d genes in the gene file.\n", (int)gene_name2id_.size());
    printf("There are %d SNPs in the gene file.\n",(int)snp_list.size());

    printf("Loading gene information %s...", ginfofile.c_str());


    FILE* fginfo=open_file((char *)ginfofile.c_str(),"r");
    int chrid,genestart,geneend;
    if (fginfo==NULL){
        fprintf(stderr, " Cannot open file %s!", ginfofile.c_str());
    }
    geneinfo_.resize(gene_name2id_.size());
	 
    while(fscanf(fginfo,"%d%d%d%s",&chrid,&genestart,&geneend,genename)==4){
        map<string,int>::iterator ret=gene_name2id_.find(genename);
        if(ret!=gene_name2id_.end()){
            geneinfo_[ret->second]=gene_info(chrid,genestart,geneend);
        }
    }


    for (size_t i=0;i<geneinfo_.size();i++){
        //          printf("%d\t%d\t%d\n", geneinfo_[i].chr_id, geneinfo_[i].start, geneinfo_[i].end);
        ;
    }
    printf("Buiding snp index ...\n");
     
    for (i=0;
         i<(int)snp_ids_.size();
         i++)
        {
            snp_set[snp_ids_[i]]=snp_pvalues_[i];
        }
    snp_ids_.clear();
    snp_pvalues_.clear();
    for (int i=0;i<(int)snp_list.size();i++){
        if (snp_set.find(snp_list[i])!=snp_set.end()){
            snp_idx_[snp_list[i]]=snp_ids_.size();
            snp_ids_.push_back(snp_list[i]);
            snp_pvalues_.push_back(snp_set[snp_list[i]]);
            snp2genes_.push_back(snp2genes[i]);
        }
    }
    printf("%d SNPs was built into index (some SNPs do not have related genes or p-values)...\n",(int)snp_ids_.size());
    printf("Loading pathway file ...\n");

    // read it first to count the number of genes in each pathway
    FILE *fdp=open_file((char *)pfile.c_str(),"r");
    if (fdp==NULL) {
        fprintf(stderr, "Cannot open file %s!", pfile.c_str());
    }
    vector<string> temp_pathway_name;
     
    vector<string> temp;
    int gene_num=0;
    while(fscanf(fdp, "%s%s", genename, pathwayname)==2){
	  
        if (gene_name2id_.find(genename)!=gene_name2id_.end())
            ;
        //     printf("$%s$\n",temp[0].c_str());
        else continue;
	       
        if (pathway_name2id_.find(pathwayname)==pathway_name2id_.end()){
            pathway_name2id_[pathwayname]=pathway_names_.size();
            pathway_names_.push_back(pathwayname);
            total_num_gene_pathway.push_back(0);
        }
        total_num_gene_pathway[pathway_name2id_[pathwayname]]++;
    }
    fclose(fdp);

     



    // read it again to filter the pathways have proper number of genes
     
    fdp=open_file((char *)pfile.c_str(),"r");
    if (fdp==NULL) {
        fprintf(stderr, "Cannot open file %s!", pfile.c_str());
    }
    pathway2genes_.clear();
    pathway2genes_.resize(total_num_gene_pathway.size());
     
    gene_num=0;
    while(fscanf(fdp, "%s%s", genename, pathwayname)==2){
	  
        if (gene_name2id_.find(genename)!=gene_name2id_.end())
            ;
        //     printf("$%s$\n",temp[0].c_str());
        else continue;
	       
        //printf("here");
        if (pathway_name2id_.find(pathwayname)==pathway_name2id_.end()){
            fprintf(stderr, "ERROR: pathway name does not exist when read the pathway file again.  at %s:%d\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        pathway2genes_[pathway_name2id_[pathwayname]].push_back(gene_name2id_[genename]);
        int num_pathways=total_num_gene_pathway[pathway_name2id_[pathwayname]];
        if (num_pathways>min_genes_ && num_pathways<max_genes_){
            int id=gene_name2id_[genename];
            gene2pathways_[id].push_back(pathway_name2id_[pathwayname]);
            gene_num++;
        }
          
    }
    for (int i=0;i<pathway2genes_.size();i++){
        sort(pathway2genes_[i].begin(),pathway2genes_[i].end());
    }
    
    
    fclose(fdp);


     

    printf("%d genes(maybe same) added to the pathway index.\n", gene_num);
    printf("%d pathways added to the pathway index.\n", (int)pathway_names_.size());
    snp2pathways_.resize(snp_pvalues_.size());

    snp_contributed_pathway_.clear();
    int id=0;
    for (i=0;i<(int)snp2genes_.size();i++){
        set<int> temp;
        for (j=0;j<(int)snp2genes_[i].size();j++){
            for (int t=0;t<(int)gene2pathways_[snp2genes_[i][j]].size();t++)
                temp.insert(gene2pathways_[snp2genes_[i][j]][t]);
        }
        for (set<int>::iterator k=temp.begin();k!=temp.end();k++){
            snp2pathways_[i].push_back(id);
        }
    }
    if (sig_gene_fraction_!=1){
        vector<double> gene_pvalue;
        gene_pvalue.resize(gene_name2id_.size(),100);
        for (size_t i=0;i<snp2genes_.size();i++){
            for (size_t j=0;j<snp2genes_[i].size();j++){
                if(gene_pvalue[snp2genes_[i][j]]>snp_pvalues_[i]){
                    gene_pvalue[snp2genes_[i][j]]=snp_pvalues_[i];
                }
            }
	       
        }
        sort(gene_pvalue.begin(), gene_pvalue.end());
        int id=(int)(sig_gene_fraction_*gene_pvalue.size());
        sig_pvalue_=gene_pvalue[id];
        printf("significant_gene_franction is enabled, the significant_pvalue is adjusted to %lf.\n", sig_pvalue_);
    }
     
    if (lfile!="NOCONFIG"){
        lfd=open_file(lfile.c_str(),"w+");
        fprintf(lfd,"This study includes %lu GO categories or pathways, %d genes and %lu snps.\n", pathway_names_.size(), gene_num, snp_ids_.size());
        fprintf(lfd,"The GO categories (pathway) information file is at %s\n", pfile.c_str());
        fprintf(lfd,"The gene file is at %s\n", gfile.c_str());
        fprintf(lfd,"The snp file is at %s\n", sfile.c_str());
        fprintf(lfd,"The significant SNP threshold is %lf\n", sig_pvalue_);
        fclose(lfd);
    }
     
}


int pathway_db::study(){
    // Obtain list of N significant genes from significant SNPs
    // Counting each gene only once

    size_t i,j,k;
    vector<int> significant_gene_table, significant_gene;
    vector<int> significant_snp;
    int N;
    vector<int> num_gene_pathway;
    vector<vector<int> > gene_list_pathway;
    boost::mt19937 rng;                 // produces randomness out of thin air
    // set random seed
    rng.seed((boost::mt19937::result_type)seed_);
    num_gene_pathway.clear();
    num_gene_pathway.resize(pathway_names_.size(),0);
    gene_list_pathway.resize(pathway_names_.size());
    // find significant gene based on pvalue threshold. 
    significant_gene_table.resize(gene_name2id_.size());
    for (i=0;i<snp_pvalues_.size();i++){
        if (snp_pvalues_[i]<sig_pvalue_+0.000000001){
            set<int> check_list;
            for (j=0;j<snp2genes_[i].size();j++){
                if (significant_gene_table[snp2genes_[i][j]]==0){
                    significant_gene_table[snp2genes_[i][j]]=1;
                    significant_gene.push_back(snp2genes_[i][j]);
                    int temp=snp2genes_[i][j];
                    for (int k=0;k<gene2pathways_[temp].size();k++){
                        if (check_list.find(gene2pathways_[temp][k])==check_list.end()){
                            if (ld_threshold_==0) {
                                num_gene_pathway[gene2pathways_[temp][k]]++;
                            }
                            else {
                                gene_list_pathway[gene2pathways_[temp][k]].push_back(temp);
                            }
                            check_list.insert(gene2pathways_[temp][k]);
                        }
                    }
                }
            }
            significant_snp.push_back(i);
        }
    }

    if (ld_threshold_!=0){
        foreach_count_significant_genes_by_ld(&num_gene_pathway, &gene_list_pathway);
    }
	 
    N=accumulate(significant_gene_table.begin(),significant_gene_table.end(),0);
    if (N<=2) {
        fprintf(stderr, "ERROR: Not enough significant genes found, stop! at %s:%d\n",
                __FILE__, __LINE__);
        exit(0);
    }
    // Count number of significant genes in each pathway

    vector<int> sim_num_gene_pathway;
    vector<double> catepvalue_pathway;
    vector<my_pathway_data > pvalue_pathway;
     
    catepvalue_pathway.resize(pathway_names_.size());
    pvalue_pathway.resize(pathway_names_.size());
    for (i=0;i<pvalue_pathway.size();i++){
        pvalue_pathway[i].catepvalue=0;
        pvalue_pathway[i].pvalue=0;
        pvalue_pathway[i].id=i;
    }
    j=0;
    for (i=0;i<num_gene_pathway.size();i++){
	  
        if (num_gene_pathway[i]>0){
            //	       printf("%d\n", num_gene_pathway[i]);
            j++;
        }
    }
    printf("%d significant genes and %lu pathways found\n", N,j);
    boost::uniform_int<> gene_rand(0,N-1);      // distribution that maps to 0..N-1
     
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
        generate_gene_id(rng, gene_rand);             
    boost::uniform_int<> sample_rand(0,permutation_number-1);      // distribution that maps to 0..permutation_number-1
     
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
        generate_sample_id(rng, sample_rand);             
     
    vector<vector<int> > sim_matrix;
    printf("prepare permutation ...\n");

    //prepare random_gene_list, and get the sim_matrix for future use.
    //sim_matrix[i] stores the number of significant gene of each pathway for random_gene_list[i]
     
    permutation_prepare(permutation_number, N, &sim_matrix);
    // use the sim_matrix to get category specific pvalue. 
    permutation_study(&sim_matrix, N, &num_gene_pathway, &catepvalue_pathway);
    calc_correlation_pathway();
    remove_correlated_pathway( catepvalue_pathway);
    for (i=0;i<catepvalue_pathway.size();i++){
        if (num_gene_pathway[i]<2) catepvalue_pathway[i]=1;
    }
    // randomly select one replicate gene list as "observe data"
    // randomly sample (with replacement)
    printf("Starting correcting errors ...\n");
    vector<double> overall_pvalues;
    vector<double > pvalue_tables;
    vector<int> sim_gene;
    sim_gene.resize(N);
    printf("\n");
    double p1=0.05;
    double p2=0.01;
    double p3=0.001;
    vector<int> v1, v2, v3;
    v1.resize(simulate_permutation_number);
    v2.resize(simulate_permutation_number);
    v3.resize(simulate_permutation_number);
    overall_pvalues.resize(simulate_permutation_number);
    pvalue_tables.resize(simulate_permutation_number*pathway_names_.size());
    int sample_id,n1,n2,n3;
    vector<double> temppvalue_pathway;
    double min;
    vector<int> rsamples;

#pragma omp parallel for default (shared) private(i,sample_id,temppvalue_pathway,k,j,min,n1,n2,n3,sim_num_gene_pathway,rng)
    for (i=0;i<(int)simulate_permutation_number;i++){
        boost::mt19937 rng;                 // produces randomness out of thin air
        // set random seed
        rng.seed((boost::mt19937::result_type)(seed_+13475983*i));
        boost::uniform_int<> sample_rand(0,permutation_number-1);      // distribution that maps to 0..permutation_number-1

        boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
            generate_sample_id(rng, sample_rand);             
        
        // randomly select a gene list
        sample_id=generate_sample_id();
        temppvalue_pathway.clear();		  temppvalue_pathway.resize(pathway_names_.size());
        sim_num_gene_pathway=sim_matrix.at(sample_id);

        // use the other gene list to get significant gene value. 
        for (j=0;j<permutation_number;j++){
            int p=generate_sample_id();
            if (p==sample_id){j--; continue;}
            int t;
            for (t=0;t<sim_num_gene_pathway.size();t++){
                if (sim_matrix.at(p).at(t)>=sim_num_gene_pathway.at(t))
                    temppvalue_pathway.at(t)+=double(1)/double(permutation_number);
            }
        }
        min=100;;
        n1=0, n2=0, n3=0;
        remove_correlated_pathway(temppvalue_pathway);
        for (j=0;j<temppvalue_pathway.size();j++){
            // ignore the pathway has one or zero significant pvalue.
            if (sim_num_gene_pathway[j]<2) temppvalue_pathway[j]=1;
            if (temppvalue_pathway[j]<=p1) {
                n1++;
                if (temppvalue_pathway[j]<=p2){
                    n2++;
                    if (temppvalue_pathway[j]<=p3)
                        n3++;
                }
            }
            if (min>temppvalue_pathway[j])
                {
					min=temppvalue_pathway[j];
                }
        }
        //v1.push_back(n1);v2.push_back(n2);v3.push_back(n3);
        v1[i]=n1;
        v2[i]=n2;
        v3[i]=n3;
        // get the minimal pvalue from this simulation.
          
        overall_pvalues[i]=min;
        copy( temppvalue_pathway.begin(), temppvalue_pathway.end(), pvalue_tables.begin()+i*pathway_names_.size());

    }

#pragma omp barrier

    sort(pvalue_tables.begin(),pvalue_tables.end());
    sort(v1.begin(),v1.end());
    sort(v2.begin(),v2.end());
    sort(v3.begin(),v3.end());
    sort(overall_pvalues.begin(), overall_pvalues.end());
     
    vector<double> :: iterator down;
    for (j=0;j<catepvalue_pathway.size();j++){
        // binary search to find P2
        down=upper_bound(overall_pvalues.begin(), overall_pvalues.end(), catepvalue_pathway[j]);
        // correct pvalue by simulation pvalue. 
        pvalue_pathway[j].pvalue=( down-overall_pvalues.begin())/(double) overall_pvalues.size();

    }

    for (i=0;i<pvalue_pathway.size();i++){
        if (num_gene_pathway[i]<=1) pvalue_pathway[i].catepvalue=1;
        else pvalue_pathway[i].catepvalue=catepvalue_pathway[i];
        pvalue_pathway[i].name=pathway_names_[i];
        pvalue_pathway[i].id=i;
    }
    sort(pvalue_pathway.begin(),pvalue_pathway.end(), sort_pred);
     
    FILE *out=open_file(ofile.c_str(), "w+");
    if (out==NULL){
        fprintf(stderr, "ERROR: can not open file at %s:%d\n",
                __FILE__, __LINE__);
        exit(0);
    }
    n1=0,n2=0,n3=0;
    for (i=0;i<pvalue_pathway.size();i++){
        if (num_gene_pathway[pvalue_pathway[i].id]>1) {
            double expected_genes_on_pathway=0;
            double expected_hits_per_study=0;
            for (j=0;j<permutation_number;j++){
                expected_genes_on_pathway+=sim_matrix[j][pvalue_pathway[i].id];
            }
            if (pvalue_pathway[i].catepvalue<=p1) {
                n1++;
                if (pvalue_pathway[i].catepvalue<=p2){
                    n2++;
                    if (pvalue_pathway[i].catepvalue<=p3)
                        n3++;
                }
            }
            expected_genes_on_pathway/=permutation_number;
            down=upper_bound(pvalue_tables.begin(), pvalue_tables.end(), pvalue_pathway[i].catepvalue+0.000000001);
            expected_hits_per_study=(down-pvalue_tables.begin());
            expected_hits_per_study/=simulate_permutation_number;
            if (pvalue_pathway[i].catepvalue<2) 
                fprintf(out, "%s\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t\n",pvalue_pathway[i].name.c_str(), total_num_gene_pathway[pvalue_pathway[i].id], num_gene_pathway[pvalue_pathway[i].id], expected_genes_on_pathway, pvalue_pathway[i].catepvalue, expected_hits_per_study, pvalue_pathway[i].pvalue);
        }
    }
    vector<int>::iterator pos;
    int m1, m2, m3;
    pos=upper_bound(v1.begin(), v1.end(), n1);
    m1=(pos-v1.begin());
    pos=upper_bound(v2.begin(), v2.end(), n2);
    m2=(pos-v2.begin());
    pos=upper_bound(v3.begin(), v3.end(), n3);
    m3=(pos-v3.begin());

    fprintf(out, "p=%lf\t number of pathways < p =%d\t significant =%lf\n", p1, n1, (double)(simulate_permutation_number-m1)/simulate_permutation_number);
    fprintf(out, "p=%lf\t number of pathways < p =%d\t significant =%lf\n", p2, n2, (double)(simulate_permutation_number-m2)/simulate_permutation_number);
    fprintf(out, "p=%lf\t number of pathways < p =%d\t significant =%lf\n", p3, n3, (double)(simulate_permutation_number-m3)/simulate_permutation_number);
    exit(0);
    return 0;
}




/*
  N is the number of significant genes
*/
int pathway_db::permutation_prepare(int permutation_number, int N, vector<vector<int> > *sim_matrix){

    boost::mt19937 rng;                 // produces randomness out of thin air
    // see pseudo-random number generators
    boost::uniform_int<> six(0,snp_ids_.size()-1);      
    rng.seed((boost::mt19937::result_type)(seed_+1564245));
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
        generate_snp_id(rng, six);             
    sim_matrix->resize(permutation_number);
    int i,j;
    vector<int> random_numbers;
    random_numbers.resize(permutation_number*N*2);
    // randomly select SNPs
    // add corresponding gene(s) to gene list.
    // stop when list contains N genes.
    // count number of simulated significant genes in each category.
    // save the result input sim_matrix for reuse in the furture
    for ( i=0;i<random_numbers.size();i++){
        random_numbers[i]=generate_snp_id();
    }
     
    vector<int> sim_gene;
    vector<int> sim_gene_id;
    vector<int> sim_snp_id;
    vector<int> is_used_snp;
    vector<int> spathways;
    vector<vector<int> > gene_list_spathways;
#pragma omp parallel for default (shared) private(i,j,sim_gene_id,sim_gene,sim_snp_id,is_used_snp,spathways,gene_list_spathways)    
    for (i=0;i<permutation_number;i++){
        int  passes=0;
        int n=0;
        int k;
        is_used_snp.clear();
        is_used_snp.resize(snp2pathways_.size(),0);
        spathways.clear();
        spathways.resize(pathway_names_.size());
        gene_list_spathways.clear();
        gene_list_spathways.resize(pathway_names_.size());
        sim_gene.clear();
        sim_gene.resize(gene_name2id_.size());
        sim_snp_id.clear(); sim_gene_id.clear();
        while(true){
            int id=random_numbers[i*2*N+passes];
            passes++;
            if (id>=snp2genes_.size())  fprintf(stderr, "ERROR: error at %s:%d\n",
                                                __FILE__, __LINE__);
            sim_gene_id.clear();
            set<int> checklist;
            for (j=0;j<snp2genes_[id].size();j++){
                if (sim_gene[snp2genes_[id][j]]==0){
                    sim_gene[snp2genes_[id][j]]=1;
                    sim_gene_id.push_back(snp2genes_[id][j]);
                    n++;
                    if (n==N) break;
                }
            }
            if (is_used_snp[id]==0) {
                is_used_snp[id]=1;
                sim_snp_id.push_back(id);
                for ( j=0;j<sim_gene_id.size();j++){
                    for (k=0;k<gene2pathways_[sim_gene_id[j]].size();k++){
                        if (checklist.find(gene2pathways_[sim_gene_id[j]][k])==checklist.end()){
                            checklist.insert(gene2pathways_[sim_gene_id[j]][k]);
                            if (ld_threshold_==0)
                                spathways[gene2pathways_[sim_gene_id[j]][k]]++;
                            else {
                                gene_list_spathways[gene2pathways_[sim_gene_id[j]][k]].push_back(sim_gene_id[j]);
                            }
                        }
                    }
                }
            }
            if (n==N) break;
        }
        if (ld_threshold_!=0){
            foreach_count_significant_genes_by_ld(&spathways, &gene_list_spathways);
        }
		  
        sim_matrix->at(i)=spathways;
    }
#pragma omp barrier     
    return 0;
}


/* sim_matrix[i][j] stores the number of significant genes in jth pathway in ith permutation */
int pathway_db::permutation_study(vector<vector< int> > *sim_matrix,int N,vector<int> * num_gene_pathway, vector<double> *catepvalue_pathway){


    int i,j;
    vector<int> sim_num_gene_pathway;
    sim_num_gene_pathway.resize(num_gene_pathway->size());
    catepvalue_pathway->clear();
    catepvalue_pathway->resize(num_gene_pathway->size());
    for (i=0;i<sim_matrix->size();i++){
        for (j=0;j<sim_num_gene_pathway.size();j++){
            if (sim_matrix->at(i)[j]>=num_gene_pathway->at(j))
                catepvalue_pathway->at(j)+=double(1)/double(sim_matrix->size());
        }
    }
    return 0;
}
int pathway_db::count_significant_pathway(vector<int> *sgenes, vector<int> *spathways){
    int i,j;
    spathways->clear();
    spathways->resize(pathway_names_.size());
    for (i=0;i<sgenes->size();i++){
        for (j=0;j<gene2pathways_[sgenes->at(i)].size();j++){
            spathways->at(gene2pathways_[sgenes->at(i)][j])++;
        }
    }
    return 0;
}

int pathway_db::count_significant_pathway_by_snp(vector<int> *ssnps, vector<int> *spathways){
    int i,j;
    spathways->clear();
    spathways->resize(pathway_names_.size());
    for (i=0;i<ssnps->size();i++){
        for (j=0;j<snp2pathways_[ssnps->at(i)].size();j++){
            spathways->at(snp2pathways_[ssnps->at(i)][j])++;
        }
	 
    }
    return 0;
}
/* This function will count the number of significant genes in a gene list,
   all genes positions less than ld_threshold differences are count as 1 */
int pathway_db::count_significant_genes_by_ld(vector<int> *sgene){

    int ret=0;
    int lastchrid=0;
    double lastpos=0;
    vector<gene_info> geneinfos;
    for (int i=0;i<sgene->size();i++){
        geneinfos.push_back(geneinfo_[sgene->at(i)]);
    }
    sort(geneinfos.begin(),geneinfos.end(),sort_gene_info);
    for (int i=0;i<geneinfos.size();i++){
        if (geneinfos[i].chr_id==0){
            ret++;
            continue;
        }
        if (lastchrid!=geneinfos[i].chr_id) {
            ret++;
            lastchrid=geneinfos[i].chr_id;
            lastpos=(geneinfos[i].start+geneinfos[i].end)/2;
        }
        else {
            double currpos=(geneinfos[i].start+geneinfos[i].end)/2;
            if (currpos-lastpos>ld_threshold_){
                ret++;
                lastpos=currpos;
            }else {
                /* else keep last position info */
            }
        }
    }
    return ret;
}

int pathway_db::foreach_count_significant_genes_by_ld(vector<int> * num_gene_pathway, vector<vector<int> > *gene_list_pathway){


    if (num_gene_pathway->size()!=gene_list_pathway->size()){
        fprintf(stderr, "ERROR: size does not match at %s:%d\n",
                __FILE__, __LINE__);
        exit(0);
    }
    for (int i=0;i<num_gene_pathway->size();i++){
        num_gene_pathway->at(i)=count_significant_genes_by_ld(&gene_list_pathway->at(i));
    }
    return 0;
}

int pathway_db::calc_correlation_pathway(){
    pathway_correlation_.clear();
    vector<int> inter(1000000);
    int pathway_size=(int)pathway2genes_.size();
    pathway_correlation_.resize(pathway_size);
    int i,j;
    int joint, total;
    vector<int>::iterator it;
    for (i=0;i<pathway_size;i++){
        pathway_correlation_[i].resize(pathway_size);
    }
    for (i=0;i<pathway_size;i++){
        for (j=0;j<pathway_size;j++){
            inter.clear();
            it=set_intersection(pathway2genes_[i].begin(),pathway2genes_[i].end(),
                                pathway2genes_[j].begin(),pathway2genes_[j].end(),
                                inter.begin());
            joint=it-inter.begin();
            total=pathway2genes_[i].size()+pathway2genes_[j].size()-joint;
            if (total==0){
                pathway_correlation_[i][j]=0;
            }
            else {
                //                    printf("here\n");
                pathway_correlation_[i][j]=double(joint)/double(total);
            }
            //if (pathway_correlation_[i][j]>0) printf("%d\n%d\n%lf\n", i,j,pathway_correlation_[i][j]);
        }
    }
    //printf("%lf\n",min);
    return 0;
     
}
int pathway_db::remove_correlated_pathway( vector<double> &catepvalue){


    vector< pair<double,int> > pvalues_idx;
    pvalues_idx.resize(catepvalue.size());
    int i,j;
    for (i=0;i<(int)catepvalue.size();i++){
        pvalues_idx[i].first=catepvalue[i];
        pvalues_idx[i].second=i;
    }
    sort(pvalues_idx.begin(),pvalues_idx.end(),sort_pair_by_first);
    for (i=0;i<(int)pvalues_idx.size();i++){
        int iid=pvalues_idx[i].second;
        if (pvalues_idx[i].first>=1) continue;
        if (catepvalue[iid]==1) continue;
        for (j=i+1;j<(int)pvalues_idx.size();j++){
            int jid=pvalues_idx[j].second;
            if (catepvalue[jid]==1) continue;
            if (pathway_correlation_[iid][jid]>max_correlation_){
                catepvalue[jid]=10;
            }
        }
    }
    return 0;

}
