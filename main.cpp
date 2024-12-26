#include <iostream>
#include <gmpxx.h>
#include <unordered_map>
#include <string>
#include <sstream>
#include <chrono>
#include <fstream>
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

unordered_map<int, mpz_class> dp;

mpz_class F(int n) {
    if (n <= 1) return 1;
    auto it = dp.find(n);
    if (it != dp.end()) return it->second;

    int k = n / 2;
    if (n % 2 == 0) return dp[n] = F(k) * F(k) + F(k - 1) * F(k - 1);
    return dp[n] = F(k) * (F(k + 1) + F(k - 1));
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
        auto res1 = fibo(n);
        auto res2 = F(n - 1);
        string s1 = res1.get_str();
        string s2 = res2.get_str();

        //cout << s1.length() << " " << s2.length() << "\n";
        if (s1.length() != s2.length()) cout << "Wrong lenght\n";
        if (res1 != res2) cout << "Wrong result\n";

        int dem = 0;
        for (int i = 0; i < s1.length(); i++) if (s1[i] != s2[i]) {
            cout << i << " " << s1[i] << " " << s2[i] << "\n";
            dem++;
            if (dem == 100) break;
        }

        if (dem > 0) {
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
    auto res1 = fibo(n);
    double cost1 = timer.getCounterMsPrecise();

    timer.startCounter();
    auto res2 = F(n - 1);
    double cost2 = timer.getCounterMsPrecise();

    string s1 = res1.get_str();
    string s2 = res2.get_str();
    int cnt = 0;
    for (int i = 0; i < s1.length(); i++) if (s1[i] != s2[i]) {
        cout << i << " " << s1[i] << " " << s2[i] << "\n";
        cnt++;
        if (cnt == 100) return 0;
    }

    if (cnt == 0) {
        cout << "binet cost = " << cost1 << "\n";
        cout << "dp cost = " << cost2 << "\n";

        timer.startCounter();
        ofstream fo("output.txt");
        fo << s1 << "\n";
        fo.close();
        cout << "Output string cost = " << timer.getCounterMsPrecise() << "\n";
    }

    
    return cnt == 0;
}


int main(int argc, char* argv[])
{
    std::ios_base::sync_with_stdio(false);
    cin.tie(0);    
    
    int L = (argc > 1) ? atoi(argv[1]) : 10000000;
    int R = (argc > 2) ? atoi(argv[2]) : L;
    R = max(L, R);

    bool result;
    if (L == R) result = test(L);
    else result = test(L, R);

    if (result) cout << "Correct\n";
    else cout << "Wrong\n";
    return 0;
}
