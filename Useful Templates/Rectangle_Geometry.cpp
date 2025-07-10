#include "macros_debug.hpp"

/*Finds the intersecting length between two intervals [a,b] and [c,d] a<=b and c<=d*/
ll inter_length(ll a, ll b, ll c, ll d){
    ll l1=b-a;
    ll l2=d-c;
    ll inter=max(min(b,d)-max(a,c),0LL);
    return l1+l2-inter;
}

/*Find the area of a rectangle.
bl - bottom left, tr - top right*/
ll area(vll &s) {
    ll bl_x = s[0], bl_y = s[1], tr_x = s[2], tr_y = s[3];
	ll length = tr_y - bl_y;
	ll width = tr_x - bl_x;
	return length * width;
}

/*Check if two rectangle intersect.
Returns true if they intersect.*/
bool intersect(vll &s1, vll &s2) {
	ll bl_a_x = s1[0], bl_a_y = s1[1], tr_a_x = s1[2], tr_a_y = s1[3];
	ll bl_b_x = s2[0], bl_b_y = s2[1], tr_b_x = s2[2], tr_b_y = s2[3];

	if (bl_a_x >= tr_b_x || tr_a_x <= bl_b_x || bl_a_y >= tr_b_y ||
	    tr_a_y <= bl_b_y) {
		return false;
	} else {
		return true;
	}
}

/*Returns intersection area of two rectangles.
-1 is returned if rectangles do not intersect*/
ll inter_area(vll &s1, vll &s2) {
    if(!intersect(s1,s2))
        return -1;
	ll bl_a_x = s1[0], bl_a_y = s1[1], tr_a_x = s1[2], tr_a_y = s1[3];
	ll bl_b_x = s2[0], bl_b_y = s2[1], tr_b_x = s2[2], tr_b_y = s2[3];

	return ((min(tr_a_x, tr_b_x) - max(bl_a_x, bl_b_x)) *
	        (min(tr_a_y, tr_b_y) - max(bl_a_y, bl_b_y)));
}

/*Returns the rectangle formed by intersection of two rectangles.
-1 is returned if there is no overlap between the rectangles*/
vll inter_rectangle(vll &s1, vll &s2){
    if(!intersect(s1,s2))
        return {};

    ll bl_a_x = s1[0], bl_a_y = s1[1], tr_a_x = s1[2], tr_a_y = s1[3];
	ll bl_b_x = s2[0], bl_b_y = s2[1], tr_b_x = s2[2], tr_b_y = s2[3];

    vll rect;
    rect.push_back(max(bl_a_x, bl_b_x));
    rect.push_back(max(bl_a_y, bl_b_y));
    rect.push_back(min(tr_a_x, tr_b_x));
    rect.push_back(min(tr_a_y, tr_b_y));
    return rect;
}

using namespace std;
int main()
{
    return 0;
}