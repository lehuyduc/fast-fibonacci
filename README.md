On Ubuntu, just type
```
chmod +x full_run.sh
./full_run.sh
chmod +x run.sh

# Output result is in output.txt
# To find a specific number, use that number as an argument: ./run.sh {N}
# Default value is N = 10^6
./run.sh 1000000

# To calculate all numbers in a range and validate the result, use: ./run.sh {L} {R}
./run 1000 100000

# To find sum of even-value Fibonacci number <= Fibo(n), use: ./run.sh -{n}
./run.sh -1000
# The above print sum of Fibonacci number with even value <= F(1000)
```

Tested on Ubuntu 18.04 and Ryzen 5900X, `best_fibo` can compute the `4*10^8`-th Fibonacci number in 0.99 second

`best_fibo` uses the same idea/algorithm as GMP's `mpz_fib_ui`, which is calculating Fibonacci numbers in group of 3.
When using 1 thread, `best_fibo` and `gmp_fibo` has the exact same speed.
But this code uses 2 threads when calculating larger `Fib(n)`

The biggest bottleneck is actually converting the result to string. It's so bad that it makes any algorithmic improvement after fast doubling irrelavent.

N = 10^6
```
dp cost = 2.91069
mpz_fib_ui cost = 2.01109
dp no-recursion cost = 2.0593
binet cost = 43.5382
matrix cost = 14.6154
cost to convert number to base10 string = 5.30343
```

N = 10^7
```
dp cost = 27.7669
mpz_fib_ui cost = 16.284
dp no-recursion cost = 16.274
binet cost = 657.771
matrix cost = 232.092
cost to convert number to base10 string = 105.358
```

N = 10^8
```
dp cost = 397.618
mpz_fib_ui cost = 268.026
dp no-recursion cost = 239.268
binet cost = 12968.6
matrix cost = 3609.47
cost to convert number to base10 string = 1904.29
Output string cost = 124.072
```

N = 4 * 10^8
```
dp cost = 1883.05
mpz_fib_ui cost = 1296.01
dp no-recursion cost = 995.325
binet cost = 64367.3
matrix cost = 18170.9
cost to convert number to base10 string = 10385.2
```

N = 10^9
```
dp cost = 5390.61
mpz_fib_ui cost = 3801.59
dp no-recursion cost = 2980.57
binet cost = 183846
matrix cost = 50054.3
cost to convert number to base10 string = 30636.7
```
