
This is just the Java implementation of XYZ2BLH program which was originally written in C. Regarding <a href="https://github.com/sasamil/XYZ2BLH">the original version</a>, I compared the performancies in various systems, OS's and builds (gcc, msvc, mingw) and the results were always the same. This is just an attempt to make test-comparison in an another language. IMO, Java is the second best choice.

Let me repeat few important things:

New algorithms are based on my (easy to derive) formula: (1+f')·x·tan(E) - b·e'^2·sin(E) - y = 0 . (*)

This formula is transformed afterwards to quartic equaton: t^4 + A·t^3 + B·t - 1 = 0 .

Four new algorithms have been created for solving this quartic equation: 1) direct solution (please, see my <a href="https://github.com/sasamil/Quartic">Quartic</a> repository) 2) optimized direct solution (taylored for this specific case) 3) Newton-Raphson based iterative solution 4) Solution based on so-called semiquadratic interpolation. The main idea behind all these methods is - speed. To simplify procedure and to reduce using of trancedential functions, in order to enhance the performancies.

New algorithms are tested against the most popular and well-known methods: Moritz-Heiskanen, Bowring (both implemented in <a href="https://en.wikipedia.org/wiki/PROJ.4">proj4</a>) and Borkowsky (recommended by <a href="https://www.iers.org">IERS</a>). The tested methods are enumerated and explained in <a href="https://github.com/sasamil/XYZ2BLH_Java/blob/master/methods.txt">'methods.txt'</a>. The comparison results are presented in <a href="https://github.com/sasamil/XYZ2BLH_Java/blob/master/results.txt">'results.txt'</a>. In a word - new methods are more than  promissing. They can cope with the old ones and they win. There is a clear logical explaination why they are so fast and I believe that they will be used everywhere where performance matters. (especially - nr2) 

And, <strong>performances are everything!</strong> They save time and energy; they save the planet's resources. <img src="https://raw.githubusercontent.com/sasamil/WMS-TMS-Maker-Qt-GUI/master/icons/emoticons/healthy_earth.jpeg" alt="healthy earth" height="28" width="42">

<div style="font-size:small;">
(*)<br/>
E  - reduced (parametric) latitude<br/>
f' - second flattening of the Earth ellipsoid<br/>
e' - second excentrity of the Earth ellipsoid<br/>
b  - polar axis<br/>
x  - abscissa of the point at the meridian cross-section: x = sqrt(X*X + Y*Y)<br/>
y  - ordinate of the point at the meridian cross-section: y = Z<br/>
</div>
