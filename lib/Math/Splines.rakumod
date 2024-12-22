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
multi sub b-spline-basis(|) is export {*}

multi sub b-spline-basis(UInt:D $d, UInt:D $i, $x, :$knots is copy = Whatever) {
    return do given $x {
        when Numeric:D { b-spline-basis-value(:$d, :$i, :$x) }
        when Whatever { b-spline-basis-function(:$d, :$i, :$knots) }
        default {
            die 'The third argument is expected to be a number or Whatever.'
        }
    }
}

multi sub b-spline-basis(UInt:D :degree(:$d)!,
                         UInt:D :index(:$i)!,
                         :arg(:argument(:$x))!,
                         :$knots is copy = Whatever) {
    return do given $x {
        when Numeric:D { b-spline-basis-value(:$d, :$knots, :$i, :$x) }
        when Whatever { b-spline-basis-function(:$d, :$knots, :$i) }
        default {
            die 'The :arg(:argument(:$x)) is expected to be a number or Whatever.'
        }
    }
}

#===========================================================
# B-spline basis value
#===========================================================

multi sub b-spline-basis-value(|) is export {*}

multi sub b-spline-basis-value(UInt:D $d, UInt:D $i, Numeric:D $x, :$knots is copy = Whatever) {
    if $knots.isa(Whatever) {
        my $k = $d + 1;
        $knots = [|(0, 1 / $k ... 1), |(1 xx ($d + 1))];
    }
    return 0 if $i > $knots.elems - $d - 2;
    return b-spline-basis-value(:$d, :$knots, :$i, :$x)
}

multi sub b-spline-basis-value(UInt:D :degree(:$d)!, :@knots!, UInt:D :index(:$i)!, Numeric:D :arg(:argument(:$x))!) {
    return 0 if $i > @knots.elems - $d - 2;
    return do if $d {
        my $l = $d - 1;
        my $j = $i + 1;
        div-z($x - @knots[$i], @knots[$i + $d] - @knots[$i]) * b-spline-basis-value(d => $l, :@knots, :$i, :$x) +
                div-z(@knots[$i + $d + 1] - $x, @knots[$i + $d + 1] - @knots[$i + 1]) * b-spline-basis-value(d => $l, :@knots, i => $j, :$x);
    } else {
        ((@knots[$i] <= $x < @knots[$i + 1]) || ($x == 1 && @knots[$i + 1] == 1)) ?? 1 !! 0
    }
}

#===========================================================
# B-spline basis function
#===========================================================

proto sub b-spline-basis-function(UInt:D :degree(:$d)!, :$knots!, UInt:D :index(:$i)!) is export {*}

multi sub b-spline-basis-function(UInt:D :degree(:$d)!,
                                  UInt:D :index(:$i)!,
                                  :$knots is copy = Whatever
                                  ) {
    # Process $knots
    $knots = do given $knots {
        when Whatever { b-spline-knots(:$d, n => $d + $i + 2 ) }
        when Numeric:D { b-spline-knots(:$d, n => $knots) }
        when $_ ~~ (Array:D | List:D | Seq:D ) && $_.all ~~ Numeric:D {
            $knots
        }
        default {
            die 'The knots argument is expected to be a positive integer, a list of numbers, or Whatever.'
        }
    }

    # Basis function index processing
    my $m = $knots.elems - $d - 2;
    die "The basis function index \$i is expected to be an integer between 0 and $m."
    unless $i ≤ $m;

    # Main computation
    return do if $d {
        my $l = $d - 1;
        my $j = $i + 1;
        -> $t {
            div-z($t - $knots[$i], $knots[$i + $d] - $knots[$i]) * b-spline-basis-function(d => $l, :$knots, :$i)($t) +
                    div-z($knots[$i + $d + 1] - $t, $knots[$i + $d + 1] - $knots[$i + 1]) * b-spline-basis-function(d => $l, :$knots, i => $j)($t);
        }
    } else {
        -> $t {
            (($knots[$i] <= $t < $knots[$i + 1]) || ($t == 1 && $knots[$i + 1] == 1)) ?? 1 !! 0
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
                         :arg(:$argument)!
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
        my $basis-value = b-spline-basis(index => $i, :$degree, :$knots, :$argument);
        $numerator[0] += $weights[$i] * $basis-value * @control-points[$i][0];
        $numerator[1] += $weights[$i] * $basis-value * @control-points[$i][1];
        $denominator += $weights[$i] * $basis-value;
    }

    return $denominator == 0 ?? [0, 0] !! [$numerator[0] / $denominator, $numerator[1] / $denominator];
}

#===========================================================
# Bernstein basis
#===========================================================

#| Represents the i^th Bernstein basis function of degree d at x.
#| C<$d> -- degree
#| C<$i> -- index
#| C<$x> -- argument
proto sub bernstein-basis(|) is export {*}

multi sub bernstein-basis(UInt:D :degree(:$d), UInt:D :index(:$i), :arg(:argument(:$x))) {
    return bernstein-basis($d, $i, $x);
}

multi sub bernstein-basis(UInt:D $d, UInt:D $i, Real:D $x) {
    return 0 if $d < $i || $x < 0 || $x > 1;
    return ([*] ($d - $i + 1) .. $d) / ([*] (1..$i)) * $x**$i * (1 - $x)**($d - $i);
}

multi sub bernstein-basis(UInt:D $d, UInt:D $i, Whatever) {
    return { 0 } if $d < $i;
    return -> $x {
        if $x < 0 || $x > 1 { 0 }
        else {
           ([*] ($d - $i + 1) .. $d) /  ([*] (1.. $i )) * $x**$i * (1 - $x)**($d - $i)
        }
    }
}