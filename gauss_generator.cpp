/*
 * writed by over
 * ABSOLUTELY NO WARRANTY!
 * omsk 2013, omsu FCS sbs001
 */
#include "gauss_generator.h"
#include <stdio.h>
#define LINEAR_START 2

/*
 * Линейный конгруэнтный генератор
 * Xk+1 = (a*Xk + c) mod m
 */
unsigned int linear_generator(int seed, unsigned int a, unsigned int c, unsigned int mod){
    unsigned val = (a*seed +c) % mod;
    return val;
    //return (val<<1)>>17;
}
/*
 * За основу берется линейная конгруэнтная последовательность. Затем man ЦПТ утверждает что сумма n членов равномерно
 * распределенной последовательности стремится к мат ожиданию последовательности с нормальным распределением
 * @param fname - имя выходного файла
 * @param base_len - кол-во членов равномерно распределенной последовательности, необходимых для генерации одного члена
 * равномерно распределенной. В литературе фигурирует число 12 как достаточное
 * @param seq_len - кол-во членов выходной последовательности
 * @param linear_a - а для генератора линейной конгруэнтной последовательности
 * @param linear_b - b -/--/-
 * @param linear_mod - mod -/--/-
 *
 */
void generate_gauss_seq(const char* fname, int base_len, int seq_len, unsigned int linear_a, unsigned int linear_c, unsigned int linear_mod){
    int i,j;
    unsigned int seed;
    double nextval;
    std::ofstream out(fname);
    seed = linear_generator(LINEAR_START, linear_a, linear_c, linear_mod);
    for(i=0;i<seq_len;i++){
	nextval = 0;
	for(j=0;j<base_len;j++){
	    seed = linear_generator(seed, linear_a, linear_c, linear_mod);
	    #ifdef DEBUG
		printf("seed=%u\n", seed);
	    #endif
	    nextval += seed;
	}
	nextval /= base_len;
	out<<nextval<<" ";
    }
    out.close();
}