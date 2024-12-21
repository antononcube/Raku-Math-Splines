use v6.d;

unit module Math::Splines;

#===========================================================
# B-spline basis
#===========================================================

multi sub b-spline-basis(|) is export {*}

multi sub b-spline-basis(Int $d, Real $x) {
    return b-spline-basis($d, 0, $x);
}

multi sub b-spline-basis(Int $d, Int $n, Real $x) {
    my $k = $d + 1;
    my @knots = 0, 1 / $k ... 1;
    return b-spline-basis(:$d, :@knots, :$n, :$x)
}

multi sub b-spline-basis(UInt:D :degree(:$d)!, :@knots!, UInt:D :index(:$n)!, Numeric:D :arg(:argument(:$x))!) {

    my $m = @knots.elems;
    return 0 if $n < 0 || $n > $m - $d - 2;

    my $t0 = @knots[$n];
    my $t1 = @knots[$n + $d + 1];
    return 1 if $d == 0 && $t0 <= $x < $t1;
    return 0 if $d == 0;

    my $left = ($x - $t0) / (@knots[$n + $d] - $t0) * b-spline-basis(degree => $d - 1, knots => @knots, index => $n, argument => $x);
    my $right = (@knots[$n + $d + 1] - $x) / (@knots[$n + $d + 1] - @knots[$n + 1]) * b-spline-basis(degree => $d - 1, knots => @knots, index => $n + 1, argument => $x);

    return $left + $right;
}

#===========================================================
# Bernstein basis
#===========================================================
