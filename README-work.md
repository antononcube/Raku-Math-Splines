# Math::Splines

Raku package for splines (piecewise polynomials), like, Bezier splines or B-splines.

------

## Installation

From [Zef ecosystem](https://raku.land):

```
zef install Math::Splines
```

From GitHub:

```
zef install https://github.com/antononcube/Raku-Math-Splines.git
```

------

## Usage examples

The subsections below show example usage for the currently implemented functions.

### B-spline basis

```raku
use Math::Splines;
use Text::Plot;
```

Using positional arguments:

```raku
b-spline-basis(3, 0, 0.5)
```

Using named arguments:

```raku
my $degree = 3; 
my @knots = b-spline-knots(:$degree, n => 6);
b-spline-basis(:$degree, :@knots, :1index, argument => 0.5)
```

Plot the basis of degree 3:

```raku
my $m = @knots.elems - $degree - 2;
my @points = (0..$m).map( -> $k { (0, 0.01 ... 1).map({ [$_, b-spline-basis(:$degree, :@knots, index => $k, argument => $_)] }) })».Array;
text-list-plot(@points, width => 80, height => 20, title => 'B-spline basis')
```

### Functions instead of values

If the argument `:arg(:argument(:$x))` of `b-spline-basis` is `Whatever`, then a function is returned:

```raku
my &bf = b-spline-basis(:3degree, :0index, arg => Whatever);

&bf(0.25)
```

Alternatively, the subs `b-spline-basis-value` and `b-spline-basis-function` can be used to get values and functions respectively.

### Bernstein basis

```raku
my @points = (^4).map( -> $k { (0, 0.01 ... 1).map({ [$_, bernstein-basis(3, $k, $_)] }) })».Array;
text-list-plot(@points, width => 80, height => 20, title => 'Bernstein basis')
```

-------

## TODO

- [ ] TODO Implementation
  - [ ] TODO Bezier curve
- [ ] TODO Documentation
  - [X] DONE Notebooks
    - [X] DONE [B-spline basis](./docs/B-spline-basis.ipynb)
    - [X] DONE [B-spline curve](./docs/B-spline-curve.ipynb)
  - [ ] TODO Demo of 2D B-spline curve with control points


-------

## References

[AAp1] Anton Antonov,
[Math::Fitting Raku package](https://github.com/antononcube/Raku-Math-Fitting),
(2024),
[GitHub/antononcube](https://github.com/antononcube).

[AAp2] Anton Antonov,
[Math::Polynomial::Chebyshev Raku package](https://github.com/antononcube/Raku-Math-Polynomial-Chebyshev),
(2024),
[GitHub/antononcube](https://github.com/antononcube).
