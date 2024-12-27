#include <iostream>
#include <gmpxx.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <chrono>
#include <map>
#include <fstream>
#include <thread>
#include <mutex>
#include <memory>
#include <vector>
#include "gmp-impl.h"
#include "longlong.h"
using namespace std;

class MyTimer {
    std::chrono::time_point<std::chrono::system_clock> start;

public:
    void startCounter() {
        start = std::chrono::system_clock::now();        
    }

    int64_t getCounterNs() {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count();        
    }

    int64_t getCounterMs() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start).count();        
    }

    double getCounterMsPrecise() {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start).count()
                / 1000000.0;        
    }
};



void
my_mpz_fib_ui (mpz_ptr fn, unsigned long n)
{
    MyTimer timer;
  mp_ptr         fp, xp, yp;
  mp_size_t      size, xalloc;
  unsigned long  n2;
  mp_limb_t      c;
  TMP_DECL;
  timer.startCounter();

  if (n <= FIB_TABLE_LIMIT)
    {
      MPZ_NEWALLOC (fn, 1)[0] = FIB_TABLE (n);
      SIZ(fn) = (n != 0);      /* F[0]==0, others are !=0 */
      return;
    }

  n2 = n/2;
  xalloc = MPN_FIB2_SIZE (n2) + 1;
  fp = MPZ_NEWALLOC (fn, 2 * xalloc);

  TMP_MARK;
  TMP_ALLOC_LIMBS_2 (xp,xalloc, yp,xalloc);  
  size = mpn_fib2_ui (xp, yp, n2);
  // cout << "GMP phase 1 cost = " << timer.getCounterMsPrecise() << "\n";

  if (n & 1)
    {
      /* F[2k+1] = (2F[k]+F[k-1])*(2F[k]-F[k-1]) + 2*(-1)^k  */
      mp_size_t  xsize, ysize;

#if HAVE_NATIVE_mpn_add_n_sub_n
      xp[size] = mpn_lshift (xp, xp, size, 1);
      yp[size] = 0;
      ASSERT_NOCARRY (mpn_add_n_sub_n (xp, yp, xp, yp, size+1));
      xsize = size + (xp[size] != 0);
      ASSERT (yp[size] <= 1);
      ysize = size + yp[size];
#else
      mp_limb_t  c2;

      c2 = mpn_lshift (fp, xp, size, 1);
      c = c2 + mpn_add_n (xp, fp, yp, size);
      xp[size] = c;
      xsize = size + (c != 0);
      c2 -= mpn_sub_n (yp, fp, yp, size);
      yp[size] = c2;
      ASSERT (c2 <= 1);
      ysize = size + c2;
#endif

      size = xsize + ysize;
      c = mpn_mul (fp, xp, xsize, yp, ysize);

#if GMP_NUMB_BITS >= BITS_PER_ULONG
      /* no overflow, see comments above */
      ASSERT (n & 2 ? fp[0] >= 2 : fp[0] <= GMP_NUMB_MAX-2);
      fp[0] += (n & 2 ? -CNST_LIMB(2) : CNST_LIMB(2));
#else
      if (n & 2)
	{
	  ASSERT (fp[0] >= 2);
	  fp[0] -= 2;
	}
      else
	{
	  ASSERT (c != GMP_NUMB_MAX); /* because it's the high of a mul */
	  c += mpn_add_1 (fp, fp, size-1, CNST_LIMB(2));
	  fp[size-1] = c;
	}
#endif
    }
  else
    {
      /* F[2k] = F[k]*(F[k]+2F[k-1]) */

      mp_size_t  xsize, ysize;
#if HAVE_NATIVE_mpn_addlsh1_n
      c = mpn_addlsh1_n (yp, xp, yp, size);
#else
      c = mpn_lshift (yp, yp, size, 1);
      c += mpn_add_n (yp, yp, xp, size);
#endif
      yp[size] = c;
      xsize = size;
      ysize = size + (c != 0);
      size += ysize;
      c = mpn_mul (fp, yp, ysize, xp, xsize);
    }

  /* one or two high zeros */
  size -= (c == 0);
  size -= (fp[size-1] == 0);
  SIZ(fn) = size;

  TMP_FREE;
}

mpz_class gmp_fibo(int n) {
    mpz_class fib_n, fib_n1;
    my_mpz_fib_ui(fib_n.get_mpz_t(), n);
    return fib_n;
}

unordered_map<int, mpz_class> dp;
mpz_class& F(int n) {
    if (n <= 1) return dp[n];
    auto it = dp.find(n);
    if (it != dp.end()) return it->second;

    int k = n / 2;
    auto Fk = F(k);
    auto Fk1 = F(k - 1);
    if (n % 2 == 0) {           
        return dp[n] = Fk * Fk + Fk1 * Fk1;
    }

    return dp[n] = Fk * (Fk1 + Fk1 + Fk);
}

mpz_class dp_fibo(int n) {
    return F(n - 1);
}

void list_dependency(map<int,int>& mp, int n) {    
    if (n <= 1 || mp[n]) return;
    mp[n] = 1;
    list_dependency(mp, n / 2);
    list_dependency(mp, n / 2 - 1);
}

void mpz_class_reserve_bits(mpz_class& x, mp_bitcnt_t bitcount)
{
    // mpz_realloc2 is a C function, so we must call it on x.get_mpz_t().
    mpz_realloc2(x.get_mpz_t(), bitcount);
}

mpz_class best_fibo(int n) {
    if (n <= 100) return gmp_fibo(n);
    //cout << "AAA" << std::endl;
    map<int,int> mp;
    list_dependency(mp, n);
    
    mpz_class f[3], dummy;
    bool started = false;
    bool flag = 0;
    map<int, mpz_class*> temps;

    for (int i = 0; i < 3; i++) mpz_class_reserve_bits(f[i], n + 64);
    mpz_class_reserve_bits(dummy, 2 * n + 64);

    cout << "n = " << n << std::endl;
    MyTimer timer;
    timer.startCounter();
    // for (auto [key, value] : mp) cout << key << "\n";
    // cout << std::endl;

    for (auto &[key, value] : mp)
        if (key >= 20 && mp.count(key - 1) && !mp.count(key - 2))
    {
        int N = key;
        //cout << "key = " << N << std::endl;

        if (!started) {
            f[0] = gmp_fibo(N - 1);
            f[1] = gmp_fibo(N);
            f[2] = f[0] + f[1];
            temps[N - 1] = &f[0];
            temps[N] = &f[1];
            temps[N + 1] = &f[2];
            started = true;
            continue;
        }

        // edge cases: 160, 170
        // 2 3 4 5, 8 9 10, 18 19 20, 38 39 40, 79 80, 160
        // 2 3 4 5, 8 9 10, 19 20 21, 41 42, 84 85, 170

        // 13 14 15, 29 30 31, 61 62 63, 125 126 252
        // 13 14 15, 29 30 31, 61 62 63, 125 126 253
        // 13 14 15, 29 30 31, 60 61 62, 123 124 125, 248 249 250, 499 500, 1000
        // F[2k] = F[k]*(F[k]+2F[k-1])
        // F[2k+1] = (2F[k]+F[k-1])*(2F[k]-F[k-1]) + 2*(-1)^k
        // OR
        // F[2k+1] = 4*F[k]^2 - F[k-1]^2 + 2*(-1)^k
        // F[2k-1] =   F[k]^2 + F[k-1]^2
        // F[2k] = F[2k+1] - F[2k-1]

        int k = N / 2;
        auto &Fk = *temps[k];
        auto &Fk1 = *temps[k - 1];                
        //cout << "k = " << k << "\n" << Fk << "\n" << Fk1 << "\n-------\n" << std::endl;
        int sign = (k % 2 == 0) ? 1 : -1;
        //cout << N << " " << k << " " <<  Fk << " " << Fk1 << std::endl;

        if (N % 2 == 1) {
            // in this case, previous F[k+1] is unused. We that to store temporary result
            auto& Fkb = *temps[k + 1];            
            // Use f[k + 1] to store F[n - 1], f[k] = F[n], F[k - 1] = F[n + 1]
            // Fkb = Fk * (Fk + 2 * Fk1); // Fkb = F(2k) = F(N - 1)
            // Fk = (2 * Fk + Fk1) * (2 * Fk - Fk1) + 2 * sign;

            // this one is correct            
            // Fkb = 4 * (Fk * Fk) - Fk1 * Fk1 + 2 * sign; // Fkb = F[2 * k + 1] = F[N]
            // Fk = Fk * (Fk + 2 * Fk1); // Fk = F[2 * K] = F[n - 1]            

            Fk *= Fk;
            Fk1 *= Fk1;
            Fkb = 4 * Fk - Fk1 + 2 * sign; // Fkb = F[2 * k + 1] = F[N]
            Fk1 += Fk; // Fk1 = F[2 * k - 1] = F[n - 2]
            Fk = Fkb - Fk1; // F[k] = F[2 * k] = F[n - 1]

            if (mp.count(N + 1)) Fk1 = Fkb + Fk;
            temps.clear();
            temps[N - 1] = &Fk;
            temps[N] = &Fkb;
            temps[N + 1] = &Fk1;
        } else {
            // in this case, F[k - 2] is unused. Use it to store F[n - 1]
            auto& Fk2 = *temps[k - 2];            
            // Fk2 = (2 * Fk1 + Fk2) * (2 * Fk1 - Fk2) - 2 * sign;            
            // Fk1 = Fk * (Fk + (Fk1 + Fk1));

            // F[2k+1] = 4*F[k]^2 - F[k-1]^2 + 2*(-1)^k
            // F[2k-1] =   F[k]^2 + F[k-1]^2
            // F[2k] = F[2k+1] - F[2k-1]
            Fk *= Fk;
            Fk1 *= Fk1;
            Fk2 = 4 * Fk - Fk1 + 2 * sign; // F2k = F[2k + 1]
            Fk1 += Fk;
            Fk = Fk2 - Fk1;
            temps[N - 1] = &Fk1;            
            temps[N] = &Fk;
            temps[N + 1] = &Fk2;


            //Fk1 = 2 * Fk1 + Fk; // F(2k) = F(k) * (F(k) + 2 * F(k - 1)) = (Fk * Fk) + (2 * Fk1 * Fk1)
            // if (mp.count(N + 1)) Fk = Fk2 + Fk1;
            // temps.clear();
            // temps[N - 1] = &Fk2;
            // temps[N] = &Fk1;
            // temps[N + 1] = &Fk;
        }
    }

    //cout << "Finish phase 1" << std::endl;
    // cout << "phase 1 cost = " << timer.getCounterMsPrecise() << "\n";
    int k = n / 2;
    auto& Fk = *temps[k];
    auto& Fk1 = *temps[k - 1];
    //cout << "final k = " << k << " " << Fk << " " << Fk1 << "\n";

    // cout << "size = " << mpz_size(Fk.get_mpz_t()) << "\n";
    // auto fibo_n = gmp_fibo(n);
    // cout << "fibo size = " << mpz_size(fibo_n.get_mpz_t()) << "\n";

    if (n % 2 == 0) return Fk * (Fk + 2 * Fk1);
    int sign = (k % 2 == 0) ? 1 : -1;
    return (2 * Fk + Fk1) * (2 * Fk - Fk1) + 2 * sign;
}

mpz_class binet_fibo(int n)
{
    // Increase default precision so we don't lose accuracy for large n.
    // A rough rule of thumb is ~ (log2(phi) * n) + a margin; here we do ~ 2*n + 32 bits.
    mp_bitcnt_t bitcount = 2 * n + 32;
    
    mpf_set_default_prec(n + 64);

    // Use mpf_class for floating-point operations
    mpf_class sqrt5(5);
    sqrt5 = sqrt(sqrt5); // sqrt(5)

    mpf_class phi = (mpf_class(1) + sqrt5) / 2; // (1 + sqrt(5)) / 2

    // power = phi^n
    mpf_class power;//(0, 2 * n + 32);
    mpf_pow_ui(power.get_mpf_t(), phi.get_mpf_t(), n);

    // result_float = power / sqrt(5) + 0.5
    mpf_class result_float = power / sqrt5;
    result_float += 0.5;

    // Convert the floating-point approximation to an integer (mpz_class)
    mpz_class result = mpz_class(result_float);

    return result;
}

bool test(int L, int R)
{
    for (int n = L; n <= R; n++) {
        auto res1 = gmp_fibo(n);
        auto res2 = best_fibo(n);
        string s1 = res1.get_str();
        string s2 = res2.get_str();
        
        if (s1.length() != s2.length()) cout << "Wrong length\n";

        if (s1 != s2) {
            cout << s1 << " " << s2 << "\n";
            cout << "Fail at n = " << n << "\n";
            return false;
        }        
    }

    cout << "Pass all\n";
    return true;
}

bool test(int n) {
    MyTimer timer;    

    timer.startCounter();
    auto res2 = dp_fibo(n);
    double cost2 = timer.getCounterMsPrecise();
    cout << "dp cost = " << cost2 << std::endl;

    timer.startCounter();
    auto res1 = gmp_fibo(n);
    double cost1 = timer.getCounterMsPrecise();    
    cout << "mpz_fib_ui cost = " << cost1 << std::endl;


    timer.startCounter();
    auto res3 = best_fibo(n);
    double cost3 = timer.getCounterMsPrecise();    
    cout << "dp no-recursion cost = " << cost3 << std::endl;    

    timer.startCounter();
    auto res4 = binet_fibo(n);
    double cost4 = timer.getCounterMsPrecise();
    cout << "binet cost = " << cost3 << std::endl;

    string s1 = res1.get_str();
    string s2 = res2.get_str();
    string s3 = res3.get_str();
    string s4 = res4.get_str();
    
    bool ok = true;
    if (s2 != s1) {cout << "DP wrong answer\n"; ok = false;};
    if (s3 != s1) {cout << "Non-recursive DP wrong answer\n"; ok = false;}
    if (s4 != s1) {cout << "Binet wrong answer\n"; ok = false;}
    
    timer.startCounter();
    ofstream fo("output.txt");
    fo << s1 << "\n";
    fo.close();
    cout << "Output string cost = " << timer.getCounterMsPrecise() << "\n";

    return ok;
}

int main(int argc, char* argv[])
{
    std::ios_base::sync_with_stdio(false);
    cin.tie(0);
    
    int L = (argc > 1) ? atoi(argv[1]) : 10000000;
    int R = (argc > 2) ? atoi(argv[2]) : L;
    R = max(L, R);

    dp[0] = 1;
    dp[1] = 1;
    dp[2] = 2;
    dp[3] = 3;

    // mpz_init(temp1);
    // mpz_init(temp2);
    // mpz_init(temp3);
    // const int bits = 100000000;
    // mpz_realloc2(temp1, bits);
    // mpz_realloc2(temp2, bits);
    // mpz_realloc2(temp3, bits);


    // for (int i = 1; i <= 10; i++) cout << F(i) << "\n";
    // return 0;

    bool result;
    if (L == R) result = test(L);
    else result = test(L, R);    

    if (result) cout << "Correct\n";
    else cout << "Wrong\n";
    return 0;
}
 
