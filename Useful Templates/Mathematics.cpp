#include "macros_debug.hpp"

/*COMPUTES SUM AND COUNT OF DIVISORS OF N*/
/*Time Complexity: O(sqrt(n))*/
pair<ll, ll> sumAndCountOfDivisors(ll n)
{
    ll cnt = 0;
    ll sum = 0;
    for (ll i = 1; i * i <= n; i++)
    {
        if (n % i == 0)
        {
            cnt++;
            sum += i;
            if (n / i != i)
            {
                cnt++;
                sum += n / i;
            }
        }
    }
    return make_pair(cnt,sum);
}

// For even better complexity use the Sum and Count formulae

/*PRIME FACTORIZATION OF N*/
/*TC=O(sqrt(n))*/
map<ll,ll> primeFactorization(ll n)
{
    map<ll,ll> pf;
    for (ll i = 2; i * i <= n; i++) while (n % i == 0) pf[i]++,n /= i;
    if (n > 1) pf[n]++;
    return pf;
}

/*SIEVE-OF-ERATOSTHENES*/
// It gives all the prime number from [1,N]
// Works well in CP for N upto 1e7
vector<bool> is_prime;
vector<ll> hp; // Highest Prime Divisor
vector<ll> lp;// Lowest Prime Divisor

void primeNumbersUptoN(ll n)
{
    is_prime.assign(n+1,true);
    hp.assign(n+1,0);
    lp.assign(n+1,0);
    is_prime[0] = is_prime[1] = 0;
    for (ll i = 2; i <= n; i++)
    {
        if (is_prime[i])
        {
            hp[i] = lp[i] = i;
            for (ll j = i * 2; j <= n; j += i)
            {
                is_prime[j] = 0;
                hp[j] = i;

                if (lp[j] == 0)
                    lp[j] = i;
            }
        }
    }
}

/*OPTIMIZED PRIME FACTORIZATION*/
/*Can be used only when we have computed lowestPrime for every n*/
/*This means this algorithm is applicable if for array elements are <=1e7 because we need to execute primeNumbersUptoN function before this*/
/*Time Complexity: O(8*23)=O(1)*/
map<ll,ll> primeFactorizationOpt(ll n)
{
    map<ll,ll> mp;
    while(n>1){
        ll x=lp[n],c=0;
        while(n%x==0) n/=x,c++;
        mp[x]=c;
    }
    return mp;
}

// TC=O(NlogN), Auxillary Space = O(1)
// Computes divisors of each number from [1,M] excluding 1 and the number itself
const ll M = 1e5 + 10;
vector<vector<ll>> divisors;
void divisorsUptoN(ll n)
{
    divisors.assign(n+1,vector<ll>(0));
    for (ll i = 2; i <= n; i++)
    {
        for (ll j = 2 * i; j <= n; j += i)
        {
            divisors[j].push_back(i);
        }
    }
}

// Calculating Number of divisors of N in O(logN)
ll countDivisors(ll x)
{
    ll ans = 1;
    while (x != 1)
    {
        ll h = hp[x];
        ll count = 0;
        while (x % h == 0)
        {
            x /= h;
            count++;
        }

        ans = ans * (count + 1);
    }
    return ans;
}

vector<ll> findAllDivisors(ll x)
{
    vector<ll> ans;
    for (ll i = 1; i * i <= x; i++)
    {
        if (x % i == 0)
        {
            ans.push_back(i);

            if (x / i != i)
                ans.push_back(x / i);
        }
    }
    return ans;
}

// Function to count number of divisors of each number till N in O(NlogN)
vector<ll> factors(ll(1e6 + 2), 1);
void precompute()
{
    ll n = factors.size();
    factors[0] = 0;
    for (ll i = 2; i < n; i++)
    {
        for (ll j = i; j < n; j += i)
        {
            factors[j]++;
        }
    }
}

/*EULER-TOTIENT FUNCTION*/
/*Input: Integer n
Output: phi(n)*/
/*Time Complexity: O(sqrt(n))*/
ll phiFunction(ll n) {
    ll result = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if (n > 1)
        result -= result / n;
    return result;
}

/*COMPUTES EULER-TOTIENT FUNCTION FROM 1 TO N*/
/*Time Complexity: O(nloglogn)*/
vector<ll> phi;
void phiUptoN(ll n) {
    phi.assign(n+1,0);
    for (int i = 0; i <= n; i++)
        phi[i] = i;

    for (int i = 2; i <= n; i++) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;
        }
    }
}

/*COMPUTES EULER-TOTIENT FUNCTION FROM 1 TO N (USING DIVISOR SUM PROPERTY)
phi(n)=n-summation(phi(d)) where d are proper divisors of n including 1*/
/*Time complexity: O(nlogn)*/
void phiUptoNv2(ll n) {
    phi[0] = 0, phi[1] = 1;
    for (int i = 2; i <= n; i++)
        phi[i] = i - 1;

    /*Note: This trick is quite common and stnadard. Obtaining results of a number using it's divisors or multiples*/
    for (int i = 2; i <= n; i++)
        for (int j = 2 * i; j <= n; j += i)
              phi[j] -= phi[i];
}


ll gcd(ll a, ll b, ll& x, ll& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    ll x1, y1;
    ll d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

bool find_any_solution(ll a, ll b, ll c, ll &x0, ll &y0, ll &g) {
    g = gcd(abs(a), abs(b), x0, y0);
    if (c % g) {
        return false;
    }

    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}

void shift_solution(ll & x, ll & y, ll a, ll b, ll cnt) {
    x += cnt * b;
    y -= cnt * a;
}

/*FINDING NUMBER OF SOLUTIONS OF ax+by=c IN THE INTERVAL x belonging to [minx,maxx] and y belonging to [miny,maxy]*/
/*If no solution exists, returns -1
NOTE: Edge cases are not handled*/
ll find_all_solutions(ll a, ll b, ll c, ll minx, ll maxx, ll miny, ll maxy) {
    ll x, y, g;
    if (!find_any_solution(a, b, c, x, y, g))
        return 0;
    a /= g;
    b /= g;

    ll sign_a = a > 0 ? +1 : -1;
    ll sign_b = b > 0 ? +1 : -1;

    shift_solution(x, y, a, b, (minx - x) / b);
    if (x < minx)
        shift_solution(x, y, a, b, sign_b);
    if (x > maxx)
        return 0;
    ll lx1 = x;

    shift_solution(x, y, a, b, (maxx - x) / b);
    if (x > maxx)
        shift_solution(x, y, a, b, -sign_b);
    ll rx1 = x;

    shift_solution(x, y, a, b, -(miny - y) / a);
    if (y < miny)
        shift_solution(x, y, a, b, -sign_a);
    if (y > maxy)
        return 0;
    ll lx2 = x;

    shift_solution(x, y, a, b, -(maxy - y) / a);
    if (y > maxy)
        shift_solution(x, y, a, b, sign_a);
    ll rx2 = x;

    if (lx2 > rx2)
        swap(lx2, rx2);
    ll lx = max(lx1, lx2);
    ll rx = min(rx1, rx2);

    if (lx > rx)
        return 0;
    return (rx - lx) / abs(b) + 1;
}

/*Handles all the edge cases*/
ll findAllSolutions(ll a, ll b, ll c, ll x1, ll x2, ll y1, ll y2){
    ll ans=0;
    if(a==0&&b==0){
        if(c==0) ans=(x2-x1+1)*(y2-y1+1);
        else ans=0;
    }
    else if(a==0){
        if(c%b==0&&c/b>=y1&&c/b<=y2) ans=(x2-x1+1);
        else ans=0;
    }
    else if(b==0){
        if(c%a==0&&c/a>=x1&&c/a<=x2) ans=(y2-y1+1);
        else ans=0;
    }
    else{
        ans=find_all_solutions(a,b,c,x1,x2,y1,y2);
    }
    return ans;
}

/*FIND POWER OF PRIME k IN n!*/
/*Time Complexity: O(logn)*/
ll factPowP(ll n, ll k) {
    ll res = 0;
    while (n) {
        n /= k;
        res += n;
    }
    return res;
}

/*FIND POWER OF COMPOSITE k IN n!*/
/*Time Complexity: O(clogn)=O(logn)*/
ll factPow(ll n, ll k){
    ll res=INF;
    map<ll,ll> pf=primeFactorization(n);
    for(auto [prime, power] : pf){
        ll frequency=factPowP(n,prime);
        res=min(res,frequency/power);
    }
    return res;
}

/*MOBIUS FUNCTION*/
/*Time Complexity: O(nlogn)*/
vector<ll> mobius;
void mobiusFunction(ll n){
    mobius.assign(n+1,0);
    mobius[1]=-1;
    for (int i=1;i<=n;i++) {
        if(!mobius[i]) continue;
		mobius[i]=-mobius[i];
		for (int j=2*i;j<=n;j+=i)
            mobius[j]+=mobius[i];
	}
}

/*SOLVES CO-PRIME QUERY IN O(k*2^k) WHERE k=Number of primes (<=9 for n=10^9)*/
/*Co-prime Query: Number of integers in the range [1,r] which are coprime with n*/
/*The same logic can be used to solve the problem: Given an array of n integers determine how many integers in the range [1,r] are multiple of atleast one of the elements in the array*/
ll solveCoprimeQuery(ll n, ll r) {
    vector<ll> p;
    for (int i=2; i*i<=n; ++i)
        if (n % i == 0) {
            p.push_back (i);
            while (n % i == 0)
                n /= i;
        }
    if (n > 1)
        p.push_back (n);

    ll sum = 0;
    for (int msk=1; msk<(1<<p.size()); ++msk) {
        ll mult = 1,bits = 0;
        for (int i=0; i<(int)p.size(); ++i)
            if (msk & (1<<i)) {
                ++bits;
                mult *= p[i];
            }

        ll cur = r / mult;
        if (bits % 2 == 1) sum += cur;
        else sum -= cur;
    }
    return r - sum;
}


/*MODULUS INVERSE WRT ANY MODULUS (NOT NECCESSARILY PRIME)*/
/*Time Complexity: O(logm) or O(logm+sqrt(m)) (Depending on whether phi is precomputed or calculated per query)*/
ll modInv(ll a, ll m){
    auto binPow=[&](ll a, ll b, ll m){
        ll ans = 1;
        while (b){
            if (b & 1) ans = (ans*1ll*a)%m;
            a = (a*1ll*a)%m;
            b >>= 1;
        }
        return ans;
    };
    return binPow(a,phiFunction(m)-1,m);
}

/*CHINESE REMAINDER THEOREM*/
/*Returns solution x mod M where M=product of mi assuming all the mi are coprime and also that M will not over flow out of ll*/
struct Congruence {
    ll a, m;
};

ll chineseRemainderTheorem(vector<Congruence> const& congruences) {
    ll M = 1, solution=0;
    for (auto const& congruence : congruences) M *= congruence.m;
    for (auto const& congruence : congruences) {
        ll a_i = congruence.a,M_i = M / congruence.m,N_i = modInv(M_i, congruence.m);
        solution = (solution + a_i * M_i % M * N_i) % M;
    }
    return solution;
}

int main()
{
    return 0;
}