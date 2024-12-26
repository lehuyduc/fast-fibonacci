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
    mpz_class fib_n;
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
        //return dp[n] = F(k) * F(k) + F(k - 1) * F(k - 1);
        if (n <= 1'000'000) return dp[n] = Fk * Fk + Fk1 * Fk1;

        mpz_class Fk1_sqr;
        vector<thread> threads;
        threads.push_back(thread([&](){
            Fk1_sqr = Fk1 * Fk1;
        }));
        dp[n] = Fk * Fk;
        threads[0].join();        
                
        return dp[n] += Fk1_sqr;
    }
    
    //return dp[n] = F(k) * (F(k + 1) + F(k - 1));
    return dp[n] = Fk * (Fk1 + Fk1 + Fk);
}

mpz_class fibo(int n)
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
        auto res2 = F(n - 1);
        string s1 = res1.get_str();
        string s2 = res2.get_str();
        
        if (s1.length() != s2.length()) cout << "Wrong length\n";

        if (s1 != s2) {
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
    auto res2 = F(n - 1);
    double cost2 = timer.getCounterMsPrecise();
    cout << "dp cost = " << cost2 << std::endl;

    timer.startCounter();
    auto res3 = fibo(n);
    double cost3 = timer.getCounterMsPrecise();    
    cout << "binet cost = " << cost3 << std::endl;    

    string s1 = res1.get_str();
    string s2 = res2.get_str();
    string s3 = res3.get_str();
    
    if (s2 != s1) cout << "DP wrong answer\n";
    if (s3 != s1) cout << "Binet wrong answer\n";
    
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

    bool result;
    if (L == R) result = test(L);
    else result = test(L, R);

    if (result) cout << "Correct\n";
    else cout << "Wrong\n";
    return 0;
}
