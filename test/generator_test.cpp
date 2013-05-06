#include "../gauss_generator.h"
#include <math.h>
int main(){
    generate_gauss_seq("test.txt", 10, 10000, 6364136223846793005L, 1442695040888963407L,((unsigned long)1)<<63);
}