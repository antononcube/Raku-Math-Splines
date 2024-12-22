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

sub b-spline-knots(*@args, *%args) is export {
    my $method = %args<method> // Whatever;

    return do given $method {
        when 'uniform' { uniform-knots(|@args, |%args) }
        when Whatever { uniform-knots(|@args, |%args) }
        default {
            die "Unknown method."
         }
    }
}

proto sub uniform-knots(|) {*}

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

proto sub b-spline-basis(UInt:D :degree(:$d)!, :$knots!, UInt:D :index(:$n)!) is export {*}

multi sub b-spline-basis(UInt:D :degree(:$d)!, UInt:D :$knots!, UInt:D :index(:$n)!) {
    my @knots = b-spline-knots(:$d, n => $knots);
    return b-spline-basis(:$d, :@knots, :$n);
}

multi sub b-spline-basis(UInt:D :degree(:$d)!, :@knots!, UInt:D :index(:$n)!) {
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

#===========================================================
# B-spline curve
#===========================================================
#`[
proto sub b-spline-curve(:@knots, :x(:@points)) is export {*}

multi sub b-spline-curve(:@knots, :x(:@points)) {
    die 'The argument :x(:@points) is expected to be an array of arrays of length 2.'
    unless @points.all ~~ (List:D | Array:D | Seq:D) && @points>>.elems.min == 2 && @points>>.elems.max == 2;

    my $n = @points.elems;
    my $m = @knots.elems - @points.elems - 1;
    my @funcs =  (^$n).map({ b-spline-basis(degree => $m, :@knots, index => $_) });

    my &res = -> $t {
        my @tbl = (^$n).map({ @funcs[$_]($t) });
        (sum(@tbl >>*<< @points>>.head), sum(@tbl >>*<< @points>>.tail))
    }

    return &res;
}
]

sub b-spline-curve-value(@control-points,
                         :$weights is copy = Whatever,
                         :$knots is copy = Whatever,
                         :d(:$degree) is copy = 3,
                         Bool:D :$closed = False,
                         :$argument!
                         ) is export {
    # Process @control-points
    die 'The argument :@contol-points is expected to be an array of arrays of length 2.'
    unless @control-points.all ~~ (List:D | Array:D | Seq:D) && @control-points>>.elems.min == 2 && @control-points>>.elems.max == 2;

    # Process $weights
    if $weights.isa(Whatever) {
        $weights = (1 xx @control-points.elems).Array
    }
    die 'The argument :$weights is expected to be Whatever or a list of numbers with the same length as the control points.'
    unless $weights ~~ (List:D | Array:D | Seq:D) && $weights.all ~~ Numeric:D;

    my $n = @control-points.elems - 1;

    if $knots.isa(Whatever) {
        $knots = uniform-knots($n, $degree)
    }
    die 'The argument :$knots is expected to be Whatever or a list of numbers.'
    unless $knots ~~ (List:D | Array:D | Seq:D) && $knots.all ~~ Numeric:D;

    my $numerator = [0, 0];
    my $denominator = 0;

    for 0..$n -> $i {
        my $basis-value = b-spline-basis-value(index => $i, :$degree, :$knots, :$argument);
        $numerator[0] += $weights[$i] * $basis-value * @control-points[$i][0];
        $numerator[1] += $weights[$i] * $basis-value * @control-points[$i][1];
        $denominator += $weights[$i] * $basis-value;
    }

    return $denominator == 0 ?? [0, 0] !! [$numerator[0] / $denominator, $numerator[1] / $denominator];
}

#===========================================================
# Bernstein basis
#===========================================================
