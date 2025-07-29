#include "macros_debug.hpp"

/*SUFFIX MEX ARRAY CREATION*/
/*Time Complexity:O(nlogn)
Returns: Suffix Mex array. suffix_mex[i]=mex of array a from [i...n-1]
0-BASED INDEXING*/
vector<ll> construct_suffix_mex(vector<ll> &a)
{
    ll n = a.size();
    vector<ll> suffix_mex(n);
    set<ll> mex_set = {0};
    set<ll> s;
    int current = 1;
    for (int i = n - 1; i >= 0; i--)
    {
        if (s.find(current) == s.end())
            mex_set.insert(current);
        current++;
        mex_set.erase(a[i]);
        s.insert(a[i]);
        suffix_mex[i] = *mex_set.begin();
    }
    return suffix_mex;
}

/*PREFIX MEX ARRAY CREATION*/
/*Time Complexity:O(nlogn)
Returns: Prefix Mex array. prefix[i]=mex of array a from [0...i] 
0-BASED INDEXING*/
vector<ll> construct_prefix_mex(vector<ll> &a){
    ll n=a.size();
    vector<ll> prefix_mex(n);
    set<ll> mex_set={0};
    set<ll> s;
    ll current=1;
    for(int i=0;i<n;i++){
        if(s.find(current)==s.end())
            mex_set.insert(current);
        current++;
        mex_set.erase(a[i]);
        s.insert(a[i]);
        prefix_mex[i]=*mex_set.begin();
    }
    return prefix_mex;
}

/*RANDOM MAPPING GENERATOR*/
/*Applications: Multiset hashing and Xor hashing problems*/

/*XOR HASHING PROBLEMS:
XOR_HASH(Array) = (hash[a1]^hash[a2]^...^hash[an])
Usually used in determining if an array contains all elements even number of times.
IF XOR_HASH(Array)=0 --> Only even elements are contained*/

/*MULTISET HASHING PROBLEMS:
MULTISET_HASH(Array) = (Sigma(hash(ai)))%INF
Usually used in checking if two arrays are equal multisets.
IF MULTISET_HASH(Array1)=MULTISET_HASH(Array2) --> Both arrays are equal multisets*/

/*In both cases discussed above, the probability of getting wrong answer is 1/2^64 which is negligible*/

mt19937_64 randomNumberGenerator(chrono::steady_clock::now().time_since_epoch().count());
class RandomMapping{
    map<ll,ull> randomMapping;
    set<ull> obtainedHashes;

    public:
    RandomMapping(){}

    /*Creates hashes correponding each array element if not present already*/
    void createHashes(vector<ll> &arr){
        for(auto &element : arr){
            ull hashValue;
            do{
                hashValue=randomNumberGenerator();
            }while(obtainedHashes.count(hashValue));
            randomMapping[element]=hashValue;
            obtainedHashes.insert(hashValue);
        }
    }

    /*Returns the hash corresponding to value*/
    ull hash(ll value){
        return randomMapping[value];
    }
};

/*Returns gcd and initializes x and y with the integers such that ax+by=gcd(a,b)*/
ll gcd(ll a, ll b, ll& x, ll& y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    ll x1, y1, d = gcd(b, a % b, x1, y1);
    x = y1, y = x1 - y1 * (a / b);
    return d;
}

/*Returns true if solution exists and initializes x0 and y0 with any solution of the equation ax+by=c and g with the gcd(a,b).
Returns false if solution does not exist
*/
bool find_any_solution(ll a, ll b, ll c, ll &x0, ll &y0, ll &g) {
    g = gcd(abs(a), abs(b), x0, y0);
    if (c % g) return false;
    x0 *= c / g, y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}

/*MINIMUM QUEUE */
// This is a queue that supports push, pop and getMin operations in O(1) time complexity
class MinQueue {
    stack<pair<ll, ll>> s1, s2;

    public:
    ll getMin() const {
        ll minimum;
        if (s1.empty() || s2.empty())
            minimum = s1.empty() ? s2.top().second : s1.top().second;
        else
            minimum = min(s1.top().second, s2.top().second);
        return minimum;
    }

    void push(ll newElement) {
        ll minimum = s1.empty() ? newElement : min(newElement, s1.top().second);
        s1.push({newElement, minimum});
    }

    ll pop(){
        if (s2.empty())
        {
            while (!s1.empty())
            {
                ll element = s1.top().first;
                s1.pop();
                ll minimum = s2.empty() ? element : min(element, s2.top().second);
                s2.push({element, minimum});
            }
        }
        ll removeElement = s2.top().first;
        s2.pop();
        return removeElement;
    }
};