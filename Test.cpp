#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;
using namespace std;
using ll = long long;
using ull = unsigned long long;
using ld = long double;
template <class T>
using oset = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
template <class T>
using omset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

#define all(v) (v).begin(), (v).end()
#define SpeedIsNeeded                 \
    ios_base::sync_with_stdio(false); \
    cin.tie(NULL);                    \
    cout.tie(NULL)

#ifndef ONLINE_JUDGE
#define debug(x)                           \
    cerr << __LINE__ << ": " << #x << " "; \
    _print(x);                             \
    cerr << endl;
#define custom_solve(t, temp)                   \
    cerr << "Test " << temp - t << ":" << endl; \
    solve();                                    \
    cerr << endl;
#include "debug.h"
#else
#define custom_solve(t, temp) solve();
#define debug(x)
#endif

template <typename T>
istream &operator>>(istream &in, vector<T> &v)
{
    for (auto &x : v)
        in >> x;
    return in;
}

template <typename T>
ostream &operator<<(ostream &out, vector<T> &v)
{
    for (auto &x : v)
        out << x << " ";
    return out;
}

namespace Constants
{
    constexpr ll N = 2e6 + 10;
    constexpr ll MOD1 = 1000000007;
    constexpr ll MOD2 = 998244353;
    constexpr ll MOD3 = ll(1e18) + 3;
    constexpr ll MOD4 = ll(1e15) + 37;
    constexpr ll INF = 1e18;
    constexpr double eps = 1e-6;
};

namespace Math
{
    vector<ll> factorials;
    void computeFactorials();

    ll ncr(ll n, ll r);

    void printD(ld a, int prec, bool bankerRound = false); // Prints rounded a to prec decimal places.
};
using namespace Math;

namespace Mod32
{
    constexpr ll mod = Constants::MOD1;
    ll binPow(ll a, ll b, ll mod = mod);
    inline ll add(ll a, ll b, ll mod = mod) { return (a + b) % mod; }
    inline ll sub(ll a, ll b, ll mod = mod) { return (a - b + mod) % mod; }
    inline ll mul(ll a, ll b, ll mod = mod) { return (a * b) % mod; }
    inline ll inv(ll a, ll mod = mod) { return binPow(a, mod - 2, mod) % mod; }
    inline ll divi(ll a, ll b, ll mod = mod) { return (a * 1LL * inv(b, mod)) % mod; }
};
using namespace Mod32;

void solve()
{
    ll n;
    cin >> n;
    vector<ll> a(n);
    cin >> a;
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("Input.txt", "r", stdin);
    freopen("Output.txt", "w", stdout);
    freopen("Error.txt", "w", stderr);
#endif
    SpeedIsNeeded;
    int t = 1;
    cin >> t;
    int temp = t;
    while (t--)
    {
        custom_solve(t, temp);
    }
    return 0;
}

namespace Math
{
    void computeFactorials()
    {
        factorials.assign(Constants::N + 1, 0);
        factorials[0] = factorials[1] = 1;
        for (int i = 2; i <= Constants::N; i++)
        {
            factorials[i] = mul(i, factorials[i - 1]);
        }
    }

    ll ncr(ll n, ll r)
    {
        ll nr = factorials[n];
        ll dr = mul(factorials[n - r], factorials[r]);
        return divi(nr, dr);
    }

    void printD(ld a, int prec, bool bankerRound)
    {
        cout << fixed << setprecision(prec);
        ld num = pow(10, prec);
        if (bankerRound)
        {
            ld temp = floor(a * num);
            if ((ll)temp % 2 == 0 && a * num - temp == 0.5)
                cout << temp / num << endl;
            else
                cout << round(a * num) / num << endl;
        }
        else
        {
            cout << round(a * num) / num << endl;
        }
    }
};

namespace Mod32
{
    ll binPow(ll a, ll b, ll mod)
    {
        ll result = 1;
        while (b)
        {
            if (b & 1)
            {
                result = mul(result, a, mod);
            }
            a = mul(a, a, mod);
            b >>= 1;
        }
        return result;
    }
};