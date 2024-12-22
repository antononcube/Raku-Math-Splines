#!/usr/bin/env raku
use v6.d;

use lib <. lib>;

use Math::Splines;

my $k = 10;
my @xs = 0, 1 / $k ... 1;

my @knots = b-spline-knots(n => 10, d => 3);
say (:@knots);
say @knots.elems;
my @knots-values = @xs.map({ b-spline-basis(:3degree, :@knots, :3index)($_) });

say (:@knots-values);