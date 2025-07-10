#include "macros_debug.hpp"

using namespace std;

const ll M = 1e9 + 7;

// 1<=a,b,M<=1e9
ll binPowRec(ll a, ll b)
{
    if (b == 0)
        return 1;

    ll res = binPowRec(a, b / 2);

    if (b & 1)
    {
        return (a * ((res * 1LL * res) % M)) % M;
    }

    return (res * 1LL * res) % M;
}

// TC=O(logN) SC=O(1)
// 1<=a,b,M<=1e9
ll binPowIter(ll a, ll b)
{
    ll ans = 1;
    while (b)
    {
        if (b & 1)
        {
            ans = (ans * 1LL * a) % M;
        }
        a = (a * 1LL * a) % M;
        b >>= 1;
    }

    return ans;
}

// TC=O(logN) SC=O(1)
// Case1: a is very large, b,M are in [1,1e9]
// Problem: a*1LL*a may overflow
ll binPowIterA(ll a, ll b)
{
    // Make the range of a [1,1e9]
    a = a % M;

    ll ans = 1;
    while (b)
    {
        if (b & 1)
        {
            ans = (ans * 1LL * a) % M;
        }
        a = (a * 1LL * a) % M;
        b >>= 1;
    }

    return ans;
}

// Case2: M is very large, a,b are in [1,1e9]
// Problem: a*1LL*a may over flow out of range of LL since they can now be as large as 10e18
ll binMul(ll a, ll b)
{
    ll ans = 0;
    while (b)
    {
        if (b & 1)
        {
            ans = (ans + a) % M;
        }
        a = (a + a) % M;
        b >>= 1;
    }

    return ans;
}

// TC=O(log^2N)
ll binPowIterM(ll a, ll b)
{
    // Make the range of a [1,1e9]
    a = a % M;

    ll ans = 1;
    while (b)
    {
        if (b & 1)
        {
            ans = binMul(ans, a) % M;
        }
        a = binMul(a, a) % M;
        b >>= 1;
    }

    return ans;
}

// Case3: b is very large, a,M are in the range [1,1e9]
// Euler's theoerem - (a^b)%M = (a^(b%phi(M)))%M, phi(M) is Euler Totient Function
// If M is prime - (a^b)%M = (a^(b%(M-1)))%M
// Look at the example in ll main
ll binPowB(ll a, ll b, ll m)
{
    ll ans = 1;
    while (b)
    {
        if (b & 1)
        {
            ans = (ans * 1LL * a) % m;
        }
        a = (a * 1LL * a) % m;
        b >>= 1;
    }

    return ans;
}

vector<vector<ll>> matAdd(vector<vector<ll>> a, vector<vector<ll>> b)
{
    /*Matrix Addition (MODULO M1=10^9+7)
    Input: matrix a-(n1xm1) and b-(n2xm2)
    Returns: a+b-(n1xm1)
    NOTE: If m1!=n1||n2!=m2, (0x0) matrix is returned
    Time Complexity: O(n1*m1)*/
    ll n1 = sz(a), n2 = sz(b), m1 = sz(a[0]), m2 = sz(b[0]);
    if (n1 != m1 || n2 != m2)
        return {{}};
    vector<vector<ll>> c(n1, vector<ll>(m1, 0));
    for (ll i = 0; i < n1; i++)
    {
        for (ll j = 0; j < m1; j++)
        {
            c[i][j] = add(a[i][j], b[i][j]);
        }
    }
    return c;
}

vector<vector<ll>> matSub(vector<vector<ll>> a, vector<vector<ll>> b)
{
    /*Matrix Addition (MODULO M1=10^9+7)
    Input: matrix a-(n1xm1) and b-(n2xm2)
    Returns: a-b-(n1xm1)
    NOTE: If m1!=n1||n2!=m2, (0x0) matrix is returned
    Time Complexity: O(n1*m1)*/

    ll n1 = a.size(), n2 = b.size(), m1 = a[0].size(), m2 = b[0].size();

    if (n1 != m1 || n2 != m2)
    {
        return {{}};
    }

    vector<vector<ll>> c(n1, vector<ll>(m1, 0));
    for (ll i = 0; i < n1; i++)
    {
        for (ll j = 0; j < m1; j++)
        {
            c[i][j] = sub(a[i][j], b[i][j]);
        }
    }
    return c;
}

vector<vector<ll>> matMul(vector<vector<ll>> a, vector<vector<ll>> b)
{
    /*Matrix Multiplication (MODULO M1=10^9+7)
    Input: matrix a-(n1xm1) and b-(n2xm2)
    Returns: axb-(n1xm2)
    NOTE: If m1!=n2, (0x0) matrix is returned
    Time Complexity: O(n1*m1*m2)*/

    ll n1 = a.size(), n2 = b.size(), m1 = a[0].size(), m2 = b[0].size();

    if (n2 != m1)
    {
        return {{}};
    }

    vector<vector<ll>> c(n1, vector<ll>(m2, 0));
    for (ll i = 0; i < n1; i++)
    {
        for (ll j = 0; j < m2; j++)
        {
            for (ll k = 0; k < n2; k++)
            {
                c[i][j] = add(c[i][j], mul(a[i][k], b[k][j]));
            }
        }
    }
    return c;
}

vector<vector<ll>> matPow(vector<vector<ll>> a, ll k)
{
    /*Matrix Exponentiation
    Input: Matrix a-(nxn), Power k
    Returns: Matrix a^k
    Time Complexity: O(n^3logk)*/
    ll n = a.size();
    vector<vector<ll>> ans(n, vector<ll>(n, 0));

    for (ll i = 0; i < n; i++)
        ans[i][i] = 1;

    while (k)
    {
        if (k & 1)
        {
            ans = matMul(ans, a);
        }
        a = matMul(a, a);
        k >>= 1;
    }
    return ans;
}

int main()
{
    ll a = 562;
    ll b = 223;

    cout << binPowIter(a, b) << endl;
    cout << binPowRec(a, b) << endl;

    cout << binPowB(50, binPowB(32, 64, M - 1), M) << endl;
    return 0;
}