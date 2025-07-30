#include <bits/stdc++.h>
using namespace std;
#include "macros_debug.hpp"

/*FUNCTION TO COMPUTE THE PREFIX FUNCTION PI/LPS (LONGEST PREFIX SEQUENCE)*/
/* Time Complexity: O(patternLength)
Space Complexity: O(patternLength)
Returns: The prefix pi function */
vector<int> prefixFunction(string &pattern)
{
    int patternLength = int(pattern.size());
    vector<int> pi(patternLength + 1);
    pi[0] = pi[1] = 0;
    for (int i = 1; i < patternLength; i++)
    {
        int j = pi[i];
        while (j > 0 && pattern[i] != pattern[j])
            j = pi[j];
        if (pattern[i] == pattern[j])
            j++;
        pi[i + 1] = j;
    }
    return pi;
}

/*KMP ALGORITHM*/
/* Time Complexity: O(n+m) */
vector<int> occurencesOfPattern(string &pattern, string &text)
{
    vector<int> pi = prefixFunction(pattern);
    vector<int> occurenceIndex;
    int textLength = text.size();
    int patternLength = pattern.size();
    int piPointer = 0;

    for (int textPointer = 0; textPointer < textLength; textPointer++)
    {
        while (piPointer && (piPointer == patternLength || text[textPointer] != pattern[piPointer]))
            piPointer = pi[piPointer];
        if (text[textPointer] == pattern[piPointer])
            piPointer++;
        if (piPointer == patternLength)
            occurenceIndex.push_back(textPointer - patternLength + 1);
    }
    return occurenceIndex;
}

/*Z-FUNCTION COMPUTATION*/
vector<int> zFunction(string s)
{
    int n = s.size();
    vector<int> z(n, 0);
    int l = 0, r = 0;
    for (int i = 1; i < n; i++)
    {
        if (i < r)
            z[i] = min(r - i, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
            z[i]++;
        if (i + z[i] > r)
            l = i, r = i + z[i];
    }
    return z;
}

/*Computing Effective Hash: Double Hash for String*/
class Hash
{
    const ll base = 31;
    vector<ll> rollingHash1;
    vector<ll> rollingHash2;

public:
    Hash() = default;
    Hash(const string &s)
    {
        ll power1 = 1, power2 = 1;
        rollingHash1.clear();
        rollingHash2.clear();
        rollingHash1.push_back(0ll);
        rollingHash2.push_back(0ll);
        for (auto c : s)
        {
            rollingHash1.push_back(add(rollingHash1.back(), mul(c - 'a' + 1, power1)));
            rollingHash2.push_back(add(rollingHash2.back(), mul(c - 'a' + 1, power2, Constants::MOD2), Constants::MOD2));
            power1 = mul(power1, base);
            power2 = mul(power2, base, Constants::MOD2);
        }
    }

    pair<ll, ll> rangeHashValue(ll l, ll r)
    {
        ll firstHash = divi(sub(rollingHash1[r + 1], rollingHash1[l]), binPow(base, l));
        ll secondHash = divi(sub(rollingHash2[r + 1], rollingHash2[l], Constants::MOD2), binPow(base, l, Constants::MOD2), Constants::MOD2);
        return make_pair(firstHash, secondHash);
    }

    pair<ll, ll> getHash() const
    {
        ll firstHash = rollingHash1.back();
        ll secondHash = rollingHash2.back();
        return make_pair(firstHash, secondHash);
    }
};