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
proto sub b-spline-basis(|) is export {*}

multi sub b-spline-basis(UInt:D $d, UInt:D $i, $x, :$knots is copy = Whatever) {
    return do given $x {
        when Numeric:D { b-spline-basis-value(:$d, :$i, :$x, :$knots) }
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

## Return the span whose non-zero basis functions contain $x.  A span is
## deliberately found only once per evaluation; the Cox--de Boor recurrence
## below can then calculate all of its ($d + 1) non-zero values together.
sub b-spline-knot-span(@knots, UInt:D $d, Numeric:D $x) {
    my $n = @knots.elems - $d - 2;
    return Nil if $n < 0 || $x < @knots[0] || $x > @knots[$n + 1];
    return $n if $x == @knots[$n + 1];

    my $low = 0;
    my $high = $n + 1;
    while $high - $low > 1 {
        my $mid = ($low + $high) div 2;
        if $x < @knots[$mid] {
            $high = $mid;
        } else {
            $low = $mid;
        }
    }
    return $low;
}

## Calculate the complete local basis instead of expanding the binary
## recursion tree for one basis value.  The lower-degree rows are retained
## only for the indices needed by the current span, giving O(d ** 2) work.
sub b-spline-active-basis-values(@knots, UInt:D $d, Numeric:D $x, Int:D $span) {
    my @previous;
    @previous[$span] = 1;

    for 1 .. $d -> $order {
        my @current;
        my $first = max(0, $span - $order);
        my $last = min($span, @knots.elems - $order - 2);
        for $first .. $last -> $i {
            @current[$i] =
                div-z($x - @knots[$i], @knots[$i + $order] - @knots[$i]) * (@previous[$i] // 0) +
                div-z(@knots[$i + $order + 1] - $x, @knots[$i + $order + 1] - @knots[$i + 1]) * (@previous[$i + 1] // 0);
        }
        @previous = @current;
    }

    my $first = max(0, $span - $d);
    my $last = min($span, @knots.elems - $d - 2);
    return ($first .. $last).map({ @previous[$_] // 0 });
}

sub b-spline-basis-value-core(UInt:D $d, @knots, UInt:D $i, Numeric:D $x) {
    return 0 if $i > @knots.elems - $d - 2;
    my $span = b-spline-knot-span(@knots, $d, $x);
    return 0 unless $span.defined;
    my $first-index = max(0, $span - $d);
    return 0 unless $i >= $first-index && $i <= $span;

    my @values = b-spline-active-basis-values(@knots, $d, $x, $span);
    return @values[$i - $first-index];
}

sub rational-de-boor-value(@control-points, @weights, @knots,
                           UInt:D $degree, Numeric:D $x, Int:D $span) {
    my @work;
    for 0 .. $degree -> $j {
        my $i = $span - $degree + $j;
        my $weight = @weights[$i];
        @work[$j] = [
            $weight * @control-points[$i][0],
            $weight * @control-points[$i][1],
            $weight,
        ];
    }

    for 1 .. $degree -> $r {
        for ($r .. $degree).reverse -> $j {
            my $i = $span - $degree + $j;
            my $alpha = div-z($x - @knots[$i], @knots[$i + $degree - $r + 1] - @knots[$i]);
            my $left = @work[$j - 1];
            my $right = @work[$j];
            @work[$j] = [
                (1 - $alpha) * $left[0] + $alpha * $right[0],
                (1 - $alpha) * $left[1] + $alpha * $right[1],
                (1 - $alpha) * $left[2] + $alpha * $right[2],
            ];
        }
    }

    my $result = @work[$degree];
    return $result[2] == 0 ?? [0, 0] !! [$result[0] / $result[2], $result[1] / $result[2]];
}

proto sub b-spline-basis-value(|) is export {*}

multi sub b-spline-basis-value(UInt:D $d, UInt:D $i, Numeric:D $x, :$knots is copy = Whatever) {
    return b-spline-basis-value(:$d, :$knots, :$i, :$x)
}

multi sub b-spline-basis-value(UInt:D :degree(:$d)!, :$knots is copy = Whatever, UInt:D :index(:$i)!, Numeric:D :arg(:argument(:$x))!) {
    if $knots.isa(Whatever) {
        my $k = $d + 1;
        $knots = [|(0, 1 / $k ... 1), |(1 xx ($d + 1))];
    }
    return b-spline-basis-value-core($d, $knots, $i, $x);
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

    # The knots and index are validated once when making the callable.  Each
    # invocation uses the same iterative evaluator as b-spline-basis-value.
    return -> $t { b-spline-basis-value-core($d, $knots, $i, $t) }
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
        $knots = uniform-knots($degree, @control-points.elems)
    }
    die 'The argument :$knots is expected to be Whatever or a list of numbers.'
    unless $knots ~~ (List:D | Array:D | Seq:D) && $knots.all ~~ Numeric:D;

    die 'The argument :$weights is expected to have the same length as the control points.'
    unless $weights.elems == @control-points.elems;

    my $span = b-spline-knot-span($knots, $degree, $argument);
    return [0, 0] unless $span.defined;
    return rational-de-boor-value(@control-points, $weights, $knots, $degree, $argument, $span)
    if $span >= $degree;

    # Preserve support before the conventional clamped domain for unusual,
    # non-clamped knot vectors.
    my @basis-values = b-spline-active-basis-values($knots, $degree, $argument, $span);
    my $first-index = max(0, $span - $degree);
    my $numerator = [0, 0];
    my $denominator = 0;

    for ^@basis-values.elems -> $offset {
        my $i = $first-index + $offset;
        next if $i < 0 || $i > $n;
        my $weighted-basis = $weights[$i] * @basis-values[$offset];
        $numerator[0] += $weighted-basis * @control-points[$i][0];
        $numerator[1] += $weighted-basis * @control-points[$i][1];
        $denominator += $weighted-basis;
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
