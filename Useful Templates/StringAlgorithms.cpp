#include <bits/stdc++.h>
using namespace std;
#include "macros_debug.hpp"

/*FUNCTION TO COMPUTE THE PREFIX FUNCTION PI/LPS (LONGEST PREFIX SEQUENCE)*/
/* Time Complexity: O(patternLength)
Space Complexity: O(patternLength)
Returns: The prefix pi function */
vector<int> prefixFunction(string& pattern) {
    int patternLength = int(pattern.size());
    vector<int> pi(patternLength + 1);
    pi[0] = pi[1] = 0;
    for (int i = 1; i < patternLength; i++) {
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
vector<int> occurencesOfPattern(string& pattern, string& text) {
    vector<int> pi = prefixFunction(pattern);
    vector<int> occurenceIndex;
    int textLength = text.size();
    int patternLength = pattern.size();
    int piPointer = 0;

    for (int textPointer = 0; textPointer < textLength; textPointer++) {
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
vector<int> zFunction(string s) {
    int n = s.size();
    vector<int> z(n,0);
    int l = 0, r = 0;
    for (int i = 1; i < n; i++) {
        if (i < r) 
            z[i] = min(r - i, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]])
            z[i]++;
        if (i + z[i] > r)
            l = i, r = i + z[i];
    }
    return z;
}

/*ROLLING HASH FUNCTION*/
/*Hash values are equal => Strings are equal*/
vector<ll> rollingHash;
void createRollingHash(string &s){
    ll len=s.size();
    rollingHash.assign(len+1,0ll);
    for(int i=1;i<=len;i++)
        rollingHash[i]=add(rollingHash[i-1],mul(s[i-1]-'a',binPow(26,i-1)));
}

/*Returns Hash value for the string s[l...r] (0-BASED INDEXING)*/
ll hashVal(ll l, ll r){
    return divi(sub(rollingHash[r+1],rollingHash[l]),binPow(26,l));
}