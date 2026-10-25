#!/usr/bin/env raku
use v6.d;

use lib 'lib';
use Math::Splines;

# Run with: raku -Ilib bench/b-spline-basis-value.raku
# Increase $samples when comparing changes on a local machine.
my $samples = 10_000;
for 3, 8, 16 -> $degree {
    my $control-points = $degree + 10;
    my @knots = b-spline-knots(:$degree, n => $control-points);
    my $index = $degree + 2;
    my $start = now;
    my $sum = 0e0;
    for ^$samples -> $step {
        my $x = $step / ($samples - 1);
        $sum += b-spline-basis-value(:$degree, :$index, :@knots, :$x);
    }
    my $elapsed = now - $start;
    say "degree $degree: $samples evaluations in {$elapsed.fmt('%.3f')}s (checksum $sum)";
}
