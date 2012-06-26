
#include "sample_arg.h"
#include "sequences.h"

using namespace arghmm;


int main(int argc, char **argv)
{
    
    char *fasta_file = argv[1];
    Sequences *sequences = read_fasta(fasta_file);
    
    int ntimes = 20;
    double times[] = {0.0, 46.238712364141541, 113.85760993922122, 212.74261506784549, 357.35077328916088, 568.82388184551303, 878.07943279474557, 1330.3307684175352, 1991.6972982819345, 2958.8711955630379, 4373.2538492690082, 6441.6288299503976, 9466.3937685609144, 13889.771066826961, 20358.461070818965, 29818.190039484831, 43651.975876399447, 63882.326155516348, 93466.929910452483, 136731.07349970445, 200000.0};
    double popsizes[ntimes];
    fill(popsizes, popsizes + ntimes, 1e4);
    double rho = 1.5e-9 * 20;
    double mu = 2.5e-9 * 20;
    
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees *trees = new LocalTrees();

    sample_arg_seq(&model, sequences, trees);
    printf("length %d\n", sequences->length());
    printf("ntrees %d\n", trees->get_num_trees());

    delete trees;
    delete sequences;

}
