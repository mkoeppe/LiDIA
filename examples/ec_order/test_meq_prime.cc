#include <LiDIA/bigint.h>
#include <LiDIA/Fp_polynomial.h>
#include <LiDIA/meq_prime.h>
#include <LiDIA/bigmod.h>

using namespace std;
using namespace LiDIA;

int main(int argc, const char * argv[]) {
  meq_prime meq;
  cout << "success? " << meq.set_prime((unsigned long) 7) << endl;

  bigmod::set_modulus(bigint(7));
  bigmod c;
  c.assign(4);
  Fp_polynomial f;
  f.set_modulus(bigint(7));
  meq.build_poly_in_Y(f, c);
  cout << f << endl;
}
