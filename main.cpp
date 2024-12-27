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

mpz_class gmp_fibo(int n) {
    mpz_class fib_n, fib_n1;
    mpz_fib_ui(fib_n.get_mpz_t(), n);
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

mpz_class best_fibo(int n) {
    if (n <= 100) return gmp_fibo(n);
    //cout << "AAA" << std::endl;
    map<int,int> mp;
    list_dependency(mp, n);
    
    mpz_class f[3][2];
    bool started = false;
    bool flag = 0;
    map<int, mpz_class*> temps;

    // cout << "n = " << n << std::endl;
    //for (auto [key, value] : mp) cout << key << "\n";
    // cout << std::endl;

    for (auto &[key, value] : mp)
        if (key >= 20 && mp.count(key - 1) && !mp.count(key - 2))
    {
        int N = key;
        //cout << "key = " << N << std::endl;

        if (!started) {
            f[0][0] = gmp_fibo(N - 1);
            f[1][0] = gmp_fibo(N);
            f[2][0] = f[0][0] + f[1][0];
            temps[N - 1] = &f[0][0];
            temps[N] = &f[1][0];
            temps[N + 1] = &f[2][0];
            flag ^= 1;
            started = true;
            continue;
        }

        // edge cases: 160, 170
        // 13 14 15, 29 30 31, 61 62 63, 125 126 252
        // 13 14 15, 29 30 31, 61 62 63, 125 126 253
        // 13 14 15, 29 30 31, 60 61 62, 123 124 125, 248 249 250, 499 500, 1000
        // F[2k] = F[k]*(F[k]+2F[k-1])
        // F[2k+1] = (2F[k]+F[k-1])*(2F[k]-F[k-1]) + 2*(-1)^k
        int k = N / 2;
        auto &Fk = *temps[k];
        auto &Fk1 = *temps[k - 1];
        //cout << "k = " << k << "\n" << Fk << "\n" << Fk1 << "\n-------\n" << std::endl;
        int sign = (k % 2 == 0) ? 1 : -1;
        if (N % 2 == 1) {
            f[0][flag] = Fk * (Fk + 2 * Fk1);
            f[1][flag] = (2 * Fk + Fk1) * (2 * Fk - Fk1) + 2 * sign;            
        } else {
            auto& Fk2 = *temps[k - 2];
            f[0][flag] = (2 * Fk1 + Fk2) * (2 * Fk1 - Fk2) - 2 * sign;
            f[1][flag] = Fk * (Fk + 2 * Fk1);
        }
        f[2][flag] = f[0][flag] + f[1][flag];
        //cout << "N = " << N << ": " << f[0][flag] << " " << f[1][flag] << " " << f[2][flag] << std::endl;

        temps.clear();
        temps[N - 1] = &f[0][flag];
        temps[N] = &f[1][flag];
        temps[N + 1] = &f[2][flag];
        flag ^= 1;
    }

    //cout << "Finish phase 1" << std::endl;
    
    int k = n / 2;
    auto& Fk = *temps[k];
    auto& Fk1 = *temps[k - 1];

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
    auto res1 = gmp_fibo(n);
    double cost1 = timer.getCounterMsPrecise();
    cout << "mpz_fib_ui cost = " << cost1 << std::endl;

    timer.startCounter();
    auto res2 = dp_fibo(n);
    double cost2 = timer.getCounterMsPrecise();
    cout << "dp cost = " << cost2 << std::endl;

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
    
    if (s2 != s1) cout << "DP wrong answer\n";
    if (s3 != s1) cout << "Non-recursive DP wrong answer\n";
    if (s4 != s1) cout << "Binet wrong answer\n";
    
    timer.startCounter();
    ofstream fo("output.txt");
    fo << s1 << "\n";
    fo.close();
    cout << "Output string cost = " << timer.getCounterMsPrecise() << "\n";

    return true;
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
