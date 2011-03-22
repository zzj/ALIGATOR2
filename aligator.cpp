#include "aligator.hpp"

int main(int argc, char * argv[])
{
     srand(123);
         srand48(time(NULL));

     pathway_db pd(argc,argv);
     pd.study();
}

