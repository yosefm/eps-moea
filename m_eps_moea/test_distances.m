% Make sure the distances function works fine.

pop = [0, 0; 0, 1; 1, 0; 1, 1];
dst = distances(pop);
correct_dist = [
	0, 1, 1, 2; 
	1, 0, 2, 1;
	1, 2, 0, 1;
	2, 1, 1, 0
	];
assert(all(dst == correct_dist));

compare_to = [-1, 0];
dst = distances(pop, compare_to);
assert(dst == [1, 2, 4, 5]')
