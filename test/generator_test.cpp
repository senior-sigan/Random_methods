#include "../gauss_generator.h"
#include <math.h>
#include <stdio.h>
int main(){
    printf("%u",((unsigned int)-1) );
    generate_gauss_seq("test.txt", 10, 10000, 69069, 5,((unsigned int)-1));
}