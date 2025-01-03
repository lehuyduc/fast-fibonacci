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

// Important functions:
// gmp_fibo
// F(int n), dp_fibo
// fast_doubling_fibo
// binet_fibo
// matrix_fibo

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

bool ktt = false;

mpz_class gmp_fibo(int n) {
    mpz_class fib_n, fib_n1;
    mpz_fib_ui(fib_n.get_mpz_t(), n);
    return fib_n;
}

unordered_map<int, mpz_class> dp;
mpz_class& F(int n) {
    if (n <= 2) return dp[n];
    auto it = dp.find(n);
    if (it != dp.end()) return it->second;

    if (dp.count(n - 1) && dp.count(n - 2)) return dp[n] = dp[n-1] + dp[n-2];
    if (dp.count(n - 1) && dp.count(n + 1)) return dp[n] = dp[n + 1] - dp[n - 1];
    if (dp.count(n + 1) && dp.count(n + 2)) return dp[n] = dp[n + 2] - dp[n + 1];

    int k = n / 2;
    auto Fk = F(k);
    auto Fk1 = F(k - 1);    

    if (n % 2 == 0) {           
        return dp[n] = Fk * (Fk + 2 * Fk1);
    }

    int sign = (k % 2 == 0) ? 1 : - 1;
    return dp[n] = (2 * Fk + Fk1) * (2 * Fk - Fk1) + 2 * sign;
}

mpz_class dp_fibo(int n) {
    return F(n);
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

void mpz_rsblsh2(mpz_class& c, mpz_class& a, mpz_class& b) {
    // Access the internal representation of a, b, and c
    auto a_mpz = a.get_mpz_t();
    auto b_mpz = b.get_mpz_t();
    auto c_mpz = c.get_mpz_t();

    // Get the limb data and sizes    
    mp_size_t a_size = abs(a_mpz->_mp_size); // Number of limbs in a
    mp_size_t b_size = abs(b_mpz->_mp_size); // Number of limbs in b

    // Compute the maximum size needed
    // Idk why it has to realloc every time. If I alloc one big block at the start, its size is reset to a small size anyway
    mp_size_t max_size = std::max(a_size, b_size) + 2; // +2 for carry/borrow    
    if (c_mpz->_mp_alloc < max_size) {
        //cout << "c_mpz resize " << std::endl;
        mpz_realloc2(c_mpz, (max_size) * GMP_LIMB_BITS);
    }
    if (a_mpz->_mp_alloc < max_size) {
        //cout << "a_mpz resize " << std::endl;
        mpz_realloc2(a_mpz, (max_size) * GMP_LIMB_BITS); 
    }
    if (b_mpz->_mp_alloc < max_size) {
        //cout << "b_mpz resize " << std::endl;
        mpz_realloc2(b_mpz, (max_size) * GMP_LIMB_BITS);
    }

    mp_ptr a_limbs = a_mpz->_mp_d;
    mp_ptr b_limbs = b_mpz->_mp_d;
    for (mp_size_t i = a_size; i < max_size; ++i) {
        a_mpz->_mp_d[i] = 0;
    }

    for (mp_size_t i = b_size; i < max_size; ++i) {
        b_mpz->_mp_d[i] = 0;
    }

    // Ensure c's limb data is clean
    std::fill(c_mpz->_mp_d + a_size, c_mpz->_mp_d + max_size, 0);
    
    mp_limb_t carry = mpn_rsblsh2_n(c_mpz->_mp_d, b_limbs, a_limbs, a_size);

    // Determine the number of significant limbs
    mp_size_t result_size = max_size;
    while (result_size > 0 && c_mpz->_mp_d[result_size - 1] == 0) {
        result_size--; // Trim trailing zeros
    }

    // Handle carry propagation correctly
    if (carry > 0) {
        c_mpz->_mp_d[result_size] = carry; // Add carry as a new highest limb
        result_size++;
    }

    // Update the size of c
    c_mpz->_mp_size = result_size;

    // Final sanity check for size mismatches
    if (result_size == 0 || result_size > max_size) {
        throw std::logic_error("Unexpected result size mismatch");
    }
}

void mpz_square_using_mpn(mpz_class &dest, const mpz_class &src) {
    // Get the raw limb pointer and size of `src`
    mp_size_t src_size = mpz_size(src.get_mpz_t());
    const mp_limb_t *src_limbs = mpz_limbs_read(src.get_mpz_t());

    // Calculate the required size for the result
    mp_size_t required_size = 2 * src_size;

    // Check if `dest` has enough space, and resize if necessary    
    if (dest.get_mpz_t()->_mp_alloc < required_size) {
        //cout << "mpz_square resize " << dest.get_mpz_t()->_mp_alloc << " " << required_size << "\n";
        mpz_realloc2(dest.get_mpz_t(), required_size * GMP_LIMB_BITS);
    }

    // Get a writable pointer for the limbs of `dest`
    mp_limb_t *dest_limbs = mpz_limbs_write(dest.get_mpz_t(), required_size);

    // Perform the square using `mpn_sqr`
    mpn_sqr(dest_limbs, src_limbs, src_size);

    // Set the size of `dest` to reflect the result
    mpz_limbs_finish(dest.get_mpz_t(), required_size);
}

const int THREAD_THRESHOLD = 500'000'000;

mpz_class fast_doubling_fibo(int n) {
    if (n <= 100) return gmp_fibo(n);
    map<int,int> mp;
    list_dependency(mp, n);
    
    mpz_class f[3], dummy;
    bool started = false;
    bool flag = 0;
    map<int, mpz_class*> temps;
    
    dummy = 0;
    mpz_class_reserve_bits(dummy, n + 64);    

    MyTimer timer;
    timer.startCounter();
    // for (auto [key, value] : mp) cout << key << "\n";
    // cout << std::endl;
    vector<thread> threads;

    for (auto &[key, value] : mp)
        if (key >= 20 && mp.count(key - 1) && !mp.count(key - 2))
    {
        int N = key;
        // cout << "key = " << N << std::endl;

        if (!started) {
            f[0] = gmp_fibo(N - 1);
            f[1] = gmp_fibo(N);
            f[2] = f[0] + f[1];
            mpz_class_reserve_bits(f[0], n + 64);
            mpz_class_reserve_bits(f[1], n + 64);
            mpz_class_reserve_bits(f[2], n + 64);
            //cout << "Initial size = " << (f[0].get_mpz_t())->_mp_alloc << "\n";
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
        int sign = (k % 2 == 0) ? 1 : -1;

        if (N % 2 == 1) {
            // in this case, previous F[k+1] is unused. We that to store temporary result
            auto& a = *temps[k - 1];
            auto& b = *temps[k];
            auto& c = *temps[k + 1];            
            // Use f[k + 1] to store F[n - 1], f[k] = F[n], F[k - 1] = F[n + 1]

            // if (n >= THREAD_THRESHOLD) {
            //     threads.clear();
            //     threads.emplace_back([&]() {
            //         Fk *= Fk;
            //         // mpz_square_using_mpn(Fkb, Fk);
            //         // Fk = std::move(Fkb);
            //     });
            //     Fk1 *= Fk1;
            //     // mpz_square_using_mpn(dummy, Fk1);
            //     // Fk1 = std::move(dummy);
            //     threads[0].join();
            // } else {
            //     // mpz_square_using_mpn(dummy, Fk);
            //     // Fk = std::move(dummy);
            //     // mpz_square_using_mpn(dummy, Fk1);
            //     // Fk1 = std::move(dummy);
            //     Fk *= Fk;
            //     Fk1 *= Fk1;
            // }

            mpz_square_using_mpn(c, b); // c = F[k]^2
            mpz_square_using_mpn(b, a); // b = F[k-1]^2
            mpz_rsblsh2(a, c, b); // a = 4 * F[k]^2 - F[k-1]^2
            if (sign > 0) mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 2); // a = F[n]
            else mpz_sub_ui(a.get_mpz_t(), a.get_mpz_t(), 2); // a = F[n]
                        
            mpz_add(b.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t()); // b = F[2k - 1] = F[n - 2]                        
            mpz_sub(c.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
            if (mp.count(N + 1)) mpz_add(b.get_mpz_t(), a.get_mpz_t(), c.get_mpz_t());

            temps[N - 1] = &c;
            temps[N] = &a;
            temps[N + 1] = &b;
        } else {
            auto &a = *temps[k - 1];
            auto &b = *temps[k];
            // in this case, F[k - 2] is unused. Use it to store F[n - 1]
            auto &c = *temps[k - 2];            

            // if (n >= THREAD_THRESHOLD) {
            //     threads.clear();
            //     threads.emplace_back([&]() {
            //         Fk *= Fk;
            //     });
            //     Fk1 *= Fk1;     
            //     threads[0].join();
            // } else {
            //     // mpz_square_using_mpn(dummy, Fk);
            //     // Fk = std::move(dummy);
            //     // mpz_square_using_mpn(dummy, Fk1);
            //     // Fk1 = std::move(dummy);
            //     Fk *= Fk;
            //     Fk1 *= Fk1;
            // }

            mpz_square_using_mpn(c, b); // c = F[k]^2
            mpz_square_using_mpn(b, a); // b = F[k-1]^2
            mpz_rsblsh2(a, c, b); // a = 4 * F[k]^2 - F[k-1]^2
            if (sign > 0) mpz_add_ui(a.get_mpz_t(), a.get_mpz_t(), 2); // a = F[2k + 1] = F[n + 1]
            else mpz_sub_ui(a.get_mpz_t(), a.get_mpz_t(), 2); // a = F[2k + 1] = F[n + 1]
            
            mpz_add(b.get_mpz_t(), b.get_mpz_t(), c.get_mpz_t()); // b = F[2k - 1] = F[n - 1]
            mpz_sub(c.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

            temps[N - 1] = &b;
            temps[N] = &c;
            temps[N + 1] = &a;
        }
    }
    //cout << "fast_doubling_fibo phase 1 = " << timer.getCounterMsPrecise() << "\n";
    int k = n / 2;
    auto& Fk = *temps[k];
    auto& Fk1 = *temps[k - 1];
    int sign = (k % 2 == 0) ? 1 : -1;

    if (n % 2 == 0) {
        //return Fk * (Fk + 2 * Fk1);               
        auto& tmp = *temps[k + 1];
        Fk1 *= 2;
        tmp = (Fk + Fk1);
        Fk1 = Fk * tmp;
        return std::move(Fk1);
        // Fk1 *= 2;
        // Fk *= (Fk + Fk1);        
        // return std::move(Fk);
    }

    //return (2 * Fk + Fk1) * (2 * Fk - Fk1) + 2 * sign;
    auto& tmp = dummy;
    Fk *= 2;
    tmp = Fk + Fk1;
    Fk -= Fk1;
    Fk1 = tmp * Fk + 2 * sign;
    return std::move(Fk1);
    // Fk *= 2;
    // Fk1 = (Fk + Fk1) * (Fk - Fk1) + 2 * sign;        
    // return std::move(Fk1);
}

unordered_map<int, mpz_class> lucas_dp;
mpz_class& lucas_dfs(int n, map<int,int> &mp) {
    if (n <= 3) return lucas_dp[n];
    auto it = lucas_dp.find(n);
    if (it != lucas_dp.end()) return it->second;

    if (lucas_dp.count(n - 1) && lucas_dp.count(n - 2)) return lucas_dp[n] = lucas_dp[n-1] + lucas_dp[n-2];
    if (lucas_dp.count(n - 1) && lucas_dp.count(n + 1)) return lucas_dp[n] = lucas_dp[n+1] - lucas_dp[n-1];
    if (lucas_dp.count(n + 1) && lucas_dp.count(n + 2)) return lucas_dp[n] = lucas_dp[n+2] - lucas_dp[n+1];

    int m = n / 2;
    auto& tmp1 = lucas_dfs(m, mp);
    auto& tmp2 = lucas_dfs(m + 1, mp);
    int q = (m % 2 == 0) ? 2 : -2;
    mpz_class u, v;
    
    if (n >= THREAD_THRESHOLD) {
        vector<thread> threads;
        threads.clear();
        threads.emplace_back([&]() {
            mpz_square_using_mpn(v, tmp2);
            v += q;
        });
        mpz_square_using_mpn(u, tmp1);
        u -= q;        
        threads[0].join();
    } else {
        mpz_square_using_mpn(u, tmp1);
        mpz_square_using_mpn(v, tmp2);
        u -= q;
        v += q;
    }

    if (n % 2 == 1) {
        // return v - u, v        
        if (mp.count(n + 1) && !lucas_dp.count(n + 1)) {
            lucas_dp[n + 1] = std::move(v);
        }
        u = lucas_dp[n + 1] - u;
        lucas_dp[n] = std::move(u);        
    } else {
        //  return u, v - u
        if (!lucas_dp.count(n + 1)) {
            v -= u;
            lucas_dp[n + 1] = std::move(v);
        }
        lucas_dp[n] = std::move(u);
    }

    if (mp.count(n + 2) && !lucas_dp.count(n + 2)) lucas_dp[n + 2] = lucas_dp[n] + lucas_dp[n + 1];
    return lucas_dp[n];
}

void list_dependency_lucas(map<int,int> &mp, int n) {
    if (mp.count(n)) return;
    mp[n] = 1;
    if (n <= 1) return;
    list_dependency_lucas(mp, n / 2);
    list_dependency_lucas(mp, n / 2 + 1); 
}

mpz_class lucas_fibo(int n) {
    if (n == 0) return 0;
    if (n <= 2) return 1;
    map<int,int> mp;
    list_dependency_lucas(mp, n);
    
    MyTimer timer;
    timer.startCounter();
    int m = n / 2;
    mpz_class u = lucas_dfs(m, mp);
    mpz_class v = lucas_dfs(m + 1, mp);
    //cout << "Lucas phase 1 cost = " << timer.getCounterMsPrecise() << "\n";

    mpz_class f = (2 * v - u) / 5;
    if (n % 2 == 1) {
        int q = (n & 2) - 1;
        return v * f - q;
    }
    return u * f;
}


mpz_class binet_fibo(int n)
{
    // Increase default precision so we don't lose accuracy for large n.    
    
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

struct Matrix {
    int n;
    vector<mpz_class> data;

    Matrix(int n) {
        this->n = n;
        data.resize(n * n);
    }

    Matrix(int n, const vector<int> &a) {
        this->n = n;
        data.resize(n * n);
        for (int i = 0; i < n * n; i++) data[i] = a[i];        
    }

    mpz_class& get(int row, int col) {
        return data[row * n + col];
    }

    const mpz_class& get(int row, int col) const {
        return data[row * n + col];
    }


    Matrix& operator*=(const Matrix& other) {
        Matrix dummy(n);

        for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            dummy.get(i, j) = 0;
            for (int k = 0; k < n; k++)
                dummy.get(i, j) += get(i, k) * other.get(k, j);
        }

        for (int i = 0; i < n * n; i++) data[i] = std::move(dummy.data[i]);        
        return *this;
    }
};

mpz_class matrix_fibo(int n) {
    if (n <= 2) return 1;
    n--;
    Matrix pow = Matrix(2, {1, 1, 1, 0});
    Matrix res = Matrix(2, {1, 0, 0, 1});
    
    while (n > 0) {
        if (n & 1) res *= pow;
        pow *= pow;
        n >>= 1;
    }
    return res.data[0];
}

bool stress_test(int times)
{
    MyTimer timer;
    double cost1 = 0, cost2 = 0;
    for (int t = 1; t <= times; t++) {
        int n = 1'000'000 + rand() % 300'000'000;
        timer.startCounter();
        auto res1 = gmp_fibo(n);
        cost1 += timer.getCounterMsPrecise();

        timer.startCounter();
        auto res2 = fast_doubling_fibo(n);
        cost2 += timer.getCounterMsPrecise();

        if (res1 != res2) {
            cout << "stress_test wrong result\n";
            return false;
        }

        if (t == 2) {
            cost1 = 0;
            cost2 = 0;
        }
        if (t >= 3)
            cout << n << " " << cost1 << " " << cost2 << " " << (cost1 / cost2) << std::endl;
    }

    cout << "stress_test success\n";
    return true;
}

bool test(int L, int R)
{
    for (int n = L; n <= R; n++) {
        cout << "n = " << n << "\n";
        auto res1 = gmp_fibo(n);
        auto res2 = fast_doubling_fibo(n);
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

    // Run this first to warm-up GMP and CPU
    timer.startCounter();
    auto res4 = dp_fibo(n);
    double cost4 = timer.getCounterMsPrecise();
    cout << "resursive doubling cost = " << cost4 << std::endl;

    timer.startCounter();    
    auto res1 = gmp_fibo(n);
    double cost1 = timer.getCounterMsPrecise();    
    cout << "mpz_fib_ui cost = " << cost1 << std::endl;

    timer.startCounter();    
    auto res2 = lucas_fibo(n);
    double cost2 = timer.getCounterMsPrecise();
    cout << "lucas_fibo cost = " << cost2 << std::endl;

    timer.startCounter();
    auto res3 = fast_doubling_fibo(n);
    double cost3 = timer.getCounterMsPrecise();    
    cout << "iterative fast doubling cost = " << cost3 << std::endl;

    //----------
    timer.startCounter();
    string s1 = res1.get_str();
    cout << "cost to convert number to base10 string = " << timer.getCounterMsPrecise() << std::endl;

    timer.startCounter();
    ofstream fo("output.txt");
    fo << s1 << "\n";
    fo.close();
    cout << "Output string cost = " << timer.getCounterMsPrecise() << std::endl;

    //----------    
    timer.startCounter();
    auto res5 = matrix_fibo(n);
    double cost5 = timer.getCounterMsPrecise();
    cout << "matrix fast expo cost = " << cost5 << std::endl;

    timer.startCounter();    
    auto res6 = binet_fibo(n);
    double cost6 = timer.getCounterMsPrecise();
    cout << "binet cost = " << cost6 << std::endl;
    
    bool ok = true;
    if (res2 != res1) {cout << "lucas_fibo WRONG ANSWER\n"; ok = false;};
    if (res3 != res1) {cout << "iterative fast doubling WRONG ANSWER\n"; ok = false;}
    if (res4 != res1) {cout << "recursive fast doubling WRONG ANSWER\n"; ok = false;}
    if (res5 != res1) {cout << "binet WRONG ANSWER\n"; ok = false;}
    if (res6 != res1) {cout << "matrix WRONG ANSWER\n"; ok = false;}    

    return ok;
}

int main(int argc, char* argv[])
{
    std::ios_base::sync_with_stdio(false);
    cin.tie(0);
    
    int L = (argc > 1) ? atoi(argv[1]) : 10000000;
    int R = (argc > 2) ? atoi(argv[2]) : L;
    R = max(L, R);

    dp[0] = 0;
    dp[1] = 1;
    dp[2] = 1;
    dp[3] = 2;

    lucas_dp[0] = 2;
    lucas_dp[1] = 1;
    lucas_dp[2] = 3;
    lucas_dp[3] = 4;
    for (int i = 4; i <= 10; i++) lucas_dp[i] = lucas_dp[i-1] + lucas_dp[i-2];

    stress_test(100);
    return 0;

    // for (int i = 0; i <= 30; i++) cout << i << " " << lucas_fibo(i) << "\n";
    // return 0;

    bool result;
    if (L == R) result = test(L);
    else result = test(L, R);    

    if (result) cout << "Correct\n";
    else cout << "WRONG\n";
    return 0;
}
 
