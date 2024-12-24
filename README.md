# Math::Splines

Raku package for splines (piecewise polynomials), like, Bezier splines or B-splines.

------

## Installation

From Zef ecosytem:

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
```
# (Any)
```

Using positional arguments:

```raku
b-spline-basis(3, 0, 0.5)
```
```
# 0.666667
```


Using named arguments:

```raku
my $degree = 3; 
my @knots = b-spline-knots(:$degree, n => 6);
b-spline-basis(:$degree, :@knots, :1index, argument => 0.5)
```
```
# 0.03125
```

Plot the basis of degree 3:

```raku
my $m = @knots.elems - $degree - 2;
my @points = (0..$m).map( -> $k { (0, 0.01 ... 1).map({ [$_, b-spline-basis(:$degree, :@knots, index => $k, argument => $_)] }) })».Array;
text-list-plot(@points, width => 120, height => 25, title => 'B-spline basis')
```
```
# B-spline basis                                                     
# +-----+---------------------+--------------------+--------------------+--------------------+--------------------+------+      
# |                                                                                                                      |      
# +     *                                                                                                         ▷      +  1.00
# |                                                                                                                      |      
# |      *                                                                                                       ▷       |      
# |       *                                                                                                     ▷        |      
# +                                                                                                                      +  0.80
# |        *                                                                                                   ▷         |      
# |          *                                                                                                ▷          |      
# |                                                                                                                      |      
# +           *      □□□□□□□               ▽▽▽ ▽▽▽▽▽▽                  ❍❍❍❍❍❍❍❍❍               ◇◇◇ ◇◇◇◇      ▷           +  0.60
# |            *   □□       □ □□       ▽▽▽▽          ▽▽▽▽▽        ❍❍❍❍❍         ❍ ❍❍❍       ◇◇◇        ◇◇   ▷            |      
# |             *□□             □□□ ▽▽▽                   ▽▽▽ ❍❍ ❍                   ❍❍❍ ◇◇◇             ◇◇▷             |      
# |             □*                ▽▽□                      ❍❍❍▽▽                       ◇◇❍                ▷◇             |      
# +            □  *             ▽▽   □□                 ❍❍❍      ▽▽▽                 ◇◇   ❍❍             ▷  ◇            +  0.40
# |           □    *          ▽▽       □□             ❍❍            ▽▽             ◇◇       ❍❍          ▷    ◇           |      
# |          □      *      ▽▽            □□        ❍❍❍                ▽▽▽       ◇ ◇           ❍❍       ▷      ◇          |      
# |        □         **  ▽▽                □□   ❍❍❍                      ▽▽▽  ◇◇                ❍❍   ▷▷        ◇         |      
# +                    ▽▽                   ❍❍ ❍□                          ◇◇◇▽                    ❍▷                    +  0.20
# |       □         ▽▽▽ **               ❍❍❍     □□□                    ◇◇◇    ▽▽ ▽              ▷ ▷ ❍❍❍        ◇        |      
# |      □        ▽▽      ***       ❍❍❍❍❍           □□□□            ◇◇◇◇           ▽▽▽▽▽      ▷▷▷       ❍❍       ◇       |      
# |           ▽▽▽▽         ❍❍ ❍❍❍❍❍❍                    □□□□◇◇◇◇ ◇◇◇                    ▽▽▷▷▷▷▽▽          ❍❍❍❍           |      
# +     ▷▷▷▷ ▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷ ▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷ ▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷ ▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷▷ ▷▷▷▷▷▷▷▷□□□□□□▽▽ ▽▽▽▽▽▽▽▽▽▽▽❍❍❍❍◇      +  0.00
# |                                                                                                                      |      
# +-----+---------------------+--------------------+--------------------+--------------------+--------------------+------+      
#       0.00                  0.20                 0.40                 0.60                 0.80                 1.00
```

### Bernstein basis

```raku
my @points = (^4).map( -> $k { (0, 0.01 ... 1).map({ [$_, bernstein-basis(3, $k, $_)] }) })».Array;
text-list-plot(@points, width => 120, height => 25, title => 'Bernstein basis')
```
```
# Bernstein basis                                                     
# +-----+---------------------+--------------------+--------------------+--------------------+--------------------+------+      
# |                                                                                                                      |      
# +     *                                                                                                         ❍      +  1.00
# |      **                                                                                                     ❍❍       |      
# |        * *                                                                                                ❍❍         |      
# |           **                                                                                            ❍❍           |      
# +             **                                                                                        ❍❍             +  0.80
# |               **                                                                                    ❍❍               |      
# |                 **                                                                                ❍❍                 |      
# |                   **                                                                            ❍❍                   |      
# +                     **                                                                       ❍ ❍                     +  0.60
# |                       ***                                                                 ❍❍❍                        |      
# |                           **                                                            ❍❍                           |      
# |                             *** □□□□□□□□□□ □□□□□□                  ▽▽▽▽▽▽▽▽▽▽ ▽▽▽▽▽▽ ❍❍❍                             |      
# +                           □□□□□□**               □□□□□□□□ ▽▽ ▽▽▽▽▽▽               ❍❍❍▽▽▽▽▽                           +  0.40
# |                      □□□□         ****              ▽▽▽▽▽▽□□ □□□              ❍❍❍❍        ▽▽▽▽                       |      
# |                  □□□□                 ***     ▽▽▽▽▽▽            □□□□□□    ❍❍❍                  ▽▽▽▽                  |      
# |                □□                       ▽▽ ▽▽▽*                      ❍❍❍❍❍□                        ▽▽                |      
# +             □□□                    ▽▽▽▽▽       *****            ❍❍❍❍❍      □□ □□□                    ▽▽▽             +  0.20
# |           □□                  ▽▽▽▽▽                 ***** ❍❍ ❍❍❍                 □□□□□                  ▽▽           |      
# |        □ □             ▽▽ ▽▽▽▽                    ❍❍❍❍❍❍❍❍** *****                    □□□□□□              ▽▽         |      
# |      □□        ▽▽▽▽▽▽▽▽             ❍❍❍❍❍❍ ❍❍❍❍❍❍❍                *********** **            □□ □□□□□□       ▽▽       |      
# +     ❍❍❍❍ ❍❍❍❍❍❍❍❍❍❍❍❍❍❍❍❍ ❍❍❍❍❍❍❍❍❍❍                                            ************** ******□□□□□□□□□▽      +  0.00
# |                                                                                                                      |      
# +-----+---------------------+--------------------+--------------------+--------------------+--------------------+------+      
#       0.00                  0.20                 0.40                 0.60                 0.80                 1.00
```

-------

## TODO

- [ ] Documentation
  - [X] Notebooks
    - [X] [B-spline basis](./docs/B-spline-basis.ipynb)
    - [X] [B-spline curve](./docs/B-spline-curve.ipynb)
  - [ ] 2D B-spline curve with control points


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
