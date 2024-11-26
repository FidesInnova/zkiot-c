#include <iostream>
#include "polynomial.h"

using namespace std;

int main () {
    int64_t p = 1588861;
    cout << "test" << endl;
    vector<int64_t> s_x = { 115, 3, 0, 0, 20, 1, 0, 17, 101, 0, 5 };
    int64_t eval_s_x = Polynomial::evaluatePolynomial(s_x, 0, p);
    int64_t alpha = Polynomial::hashAndExtractLower4Bytes(eval_s_x, p);

    cout << "eval_s_x: " << eval_s_x << endl;
    cout << "alpha: " << alpha << endl;

    return 0;
}