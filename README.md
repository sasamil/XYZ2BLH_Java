
This is just the Java implementation of XYZ2BLH program which was originally written in C. Regarding the original version (github.com/sasamil/XYZ2BLH), I compared the performancies in various systems, OS's and builds (gcc, msvc, mingw) and the results were always the same. This is just an attempt to make test-comparison in an another language. IMO, Java is the second best choice.

Let me repeat few important things:

New algorithms are based on my (easy to derive) formula: (1+f')·X·tan(E) - b·e'^2·sin(E) - Y = 0 . 

This formula is transformed afterwards to quartic equaton: t^4 + A·t^3 + B·t - 1 = 0 .

Four algorithms are created for solving this quartic equation: 1) direct solution (please, see my 'Quartic' repository) 2) optimized direct solution (taylored for this specific case) 3) Newton-Raphson based iterative solution 4) Solution based on so-called semiquadratic interpolation. The main idea behind all these methods is - speed. To simplify procedure and to reduce using of trancedential functions, in order to enhance the performancies.

New algorithms are tested against the most popular and well-known methods: Moritz-Heiskanen (implemented in proj), Bowring, Borkowsky (recommended by IERS) <!--img src="http://forum.srpskinacionalisti.com/images/smilies/eusa_naughty.gif" alt="not recommended" height="16" width="20"-->. The tested methods are enumerated and explained in <i>'methods.txt'</i>. The comparison results are presented in <i>'results.txt'</i>. In a word - new methods are more than  promissing. They can cope with the old ones and they win. There is a clear logical explaination why they are so fast(er) and I believe that they will be used everywhere where performance matters. (especially - nr2) 

And, <strong>performances are everything!</strong> They save time and energy; they save the planet's resources. <img src="https://encrypted-tbn2.gstatic.com/images?q=tbn:ANd9GcQO6YPxiFiU933Zvsw6Ihd-jSC1NZMYghcNbxA1G-npmkx0G1nXPQ" alt="it's our home" height="40" width="40">
