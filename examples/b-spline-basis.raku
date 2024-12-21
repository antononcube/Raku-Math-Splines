#!/usr/bin/env raku
use v6.d;

use lib <. lib>;

use Math::Splines;

my $k = 10;
my @xs = 0, 1 / $k ... 1;
my @values = @xs.map({ b-spline-basis(3, 0, $_) });

say (:@values);

#my @knots = 0, 0, 0, 0, 1/10, 1/5, 3/10, 2/5, 1/2, 3/5, 7/10, 4/5, 9/10, 1, 1, 1, 1;
my @knots = 0, 1/10, 1/5, 3/10, 2/5, 1/2, 3/5, 7/10, 4/5, 9/10, 1;
my @knots-values = @xs.map({ b-spline-basis(d => 3, :@knots, :0index, arg => $_) });

say (:@knots-values);