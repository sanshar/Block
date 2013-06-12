/*
 * Simple algorithm for enumerating even permutations of a set.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _EPERMUTE_H
#define _EPERMUTE_H

#include <algorithm>


/* Templates */

// Permutes the range [first,last) into the lexicographically next even
// permutation of the elements.  Starting from a sorted range and calling this
// function repeatedly will enumerate the (N!)/2 even permutations of the
// range.
//
// IMPORTANT: this function assumes that all elements in the range are
// distinct, because only then are even permutations distinct from all
// permutations. The behaviour of this function is undefined if there are
// duplicate elements in the range.
//
// Returns false if there is no lexicographically next permutation, in which
// case the range is transformed into the lexicographically smallest
// permutation; otherwise returns true.
//
// Amortized time complexity is linear; space complexity is O(1).
template <typename BidirectionalIterator>
bool next_even_permutation(BidirectionalIterator first,
			BidirectionalIterator last);

// Same as the previous function, except that the strict weak ordering comp is
// used to determine lexicographic ordering, instead of operator<().
template <typename BidirectionalIterator, typename StrictWeakOrdering>
bool next_even_permutation(BidirectionalIterator first,
			BidirectionalIterator last, StrictWeakOrdering comp);


/* Template implementations */

template <typename BidirectionalIterator>
bool next_even_permutation(BidirectionalIterator first,
			BidirectionalIterator last)
{
	// Base cases: ranges of 0 or 1 element have no distinct permutations.
	// (Neither do 2-element ranges have distinct even permutations; but
	// that's taken care of by the main loop.)
	if (first == last)
		return false;
	BidirectionalIterator i = first;
	++i;
	if (i == last)
		return false;

	// Enumerate all permutations until we reach an even one.
	bool parity = false;
	bool ret = true;
	do {
		// Find last increasing pair of elements in the range.
		i = last;
		BidirectionalIterator j = --i;
		--i;

		long n = 1;
		while (i != first && *j < *i) {
			j = i;
			++n;
			--i;
		}

		if (*j < *i) {
			// Entire range is decreasing: it's the
			// lexicographically greatest. So wrap around.
			j = first;
			n++;
			ret = false;
		} else {
			// Find last element larger than *i and swap them.
			BidirectionalIterator k = last;
			while (!(*i < *--k));
			std::iter_swap(i, k);
			parity = !parity;
		}

		// Reverse last decreasing range to get lexicographically next
		// smallest permutation.
		std::reverse(j, last);
		if ((n / 2) % 2 == 1)
			parity = !parity;
	} while (parity);

	return ret;
}

template <typename BidirectionalIterator, typename StrictWeakOrdering>
bool next_even_permutation(BidirectionalIterator first,
			BidirectionalIterator last, StrictWeakOrdering comp)
{
	// Base cases: ranges of 0 or 1 element have no distinct permutations.
	// (Neither do 2-element ranges have distinct even permutations; but
	// that's taken care of by the main loop.)
	if (first == last)
		return false;
	BidirectionalIterator i = first;
	++i;
	if (i == last)
		return false;

	// Enumerate all permutations until we reach an even one.
	bool parity = false;
	bool ret = true;
	do {
		// Find last increasing pair of elements in the range.
		i = last;
		BidirectionalIterator j = --i;
		--i;

		long n = 1;
		while (i != first && comp(*j, *i)) {
			j = i;
			++n;
			--i;
		}

		if (comp(*j, *i)) {
			// Entire range is decreasing: it's the
			// lexicographically greatest. So wrap around.
			j = first;
			n++;
			ret = false;
		} else {
			// Find last element larger than *i and swap them.
			BidirectionalIterator k = last;
			while (!comp(*i, *--k));
			std::iter_swap(i, k);
			parity = !parity;
		}

		// Reverse last decreasing range to get lexicographically next
		// smallest permutation.
		std::reverse(j, last);
		if ((n / 2) % 2 == 1)
			parity = !parity;
	} while (parity);

	return ret;
}


#endif // _EPERMUTE_H
