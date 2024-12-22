use v6.d;

unit module Math::Splines;

#===========================================================
# Utilities
#===========================================================

sub div-z(Numeric:D $s, Numeric:D $t) {
    return $t == 0 ?? 0 !! $s / $t;
}

sub power-z($x, $y) {
    return $x == 0 && $y == 0 ?? 1 !! $x ** $y;
}

#-----------------------------------------------------------

proto sub uniform-knots(|) is export {*}

multi sub uniform-knots(UInt:D :degree(:$d), UInt:D :$n) is export {
    my @a = |(0 xx ($d + 1)), |(1 .. $n - $d - 1).map({ ($_ / ($n - $d)) }), |(1 xx ($d + 1));
    return @a;
}

multi sub uniform-knots(UInt:D $d, UInt:D $n) is export {
    return uniform-knots(:$d, :$n);
}

#===========================================================
# B-spline basis
#===========================================================

multi sub b-spline-basis-value(|) is export {*}

multi sub b-spline-basis-value(UInt:D $d, UInt:D $n, Real $x) {
    my $k = $d + 1;
    my @knots = |(0, 1 / $k ... 1), |(1 xx ($d + 1));
    return 0 if $n > @knots.elems - $d - 2;
    return b-spline-basis-value(:$d, :@knots, :$n, :$x)
}

multi sub b-spline-basis-value(UInt:D :degree(:$d)!, :@knots!, UInt:D :index(:$n)!, Numeric:D :arg(:argument(:$x))!) {
    return 0 if $n > @knots.elems - $d - 2;
    return do if $d {
        my $l = $d - 1;
        my $j = $n + 1;
        div-z($x - @knots[$n], @knots[$n + $d] - @knots[$n]) * b-spline-basis-value(d => $l, :@knots, :$n, :$x) +
                div-z(@knots[$n + $d + 1] - $x, @knots[$n + $d + 1] - @knots[$n + 1]) * b-spline-basis-value(d => $l, :@knots, n => $j, :$x);
    } else {
        ((@knots[$n] <= $x < @knots[$n + 1]) || ($x == 1 && @knots[$n + 1] == 1)) ?? 1 !! 0
    }
}

#===========================================================
# B-spline basis functions
#===========================================================

sub b-spline-basis(UInt:D :degree(:$d), :@knots, UInt:D :index(:$n)) is export {
    my $m = @knots.elems - $d - 2;
    die "The basis function index \$n is expected to be an integer between 0 and $m."
    unless $n ≤ $m;

    return do if $d {
        my $l = $d - 1;
        my $j = $n + 1;
        -> $t {
            div-z($t - @knots[$n], @knots[$n + $d] - @knots[$n]) * b-spline-basis(d => $l, :@knots, :$n)($t) +
                    div-z(@knots[$n + $d + 1] - $t, @knots[$n + $d + 1] - @knots[$n + 1]) * b-spline-basis(d => $l, :@knots, n => $j)($t);
        }
    } else {
        -> $t {
            ((@knots[$n] <= $t < @knots[$n + 1]) || ($t == 1 && @knots[$n + 1] == 1)) ?? 1 !! 0
        }
    }
}

sub b-spline-curve(:@knots, :@x) is export {
    my $n = @x.elems;
    my $m = @knots.elems - @x.elems - 1;

    my &res = -> $t {
        [\+] (^$n).map({ b-spline-basis($_, $m, @knots)($t) * @x[$_] });
    }

    return &res;
}

#===========================================================
# Bernstein basis
#===========================================================
