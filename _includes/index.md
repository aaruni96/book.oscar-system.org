{% assign versions = site.data.versions %}

{% for mv in versions %}
    {% if page.url contains mv %}
        {% assign pv = mv %}
    {% endif %}
{% endfor %}

# {{ site.title }}

The book _The Computer Algebra System OSCAR_  presents version 1.0 of the new [Computer Algebra
System OSCAR](https://oscar-system.org). The book is an invitation to use OSCAR.  With discussions
of theoretical and algorithmic aspects included, it offers a multitude of explicit code snippets.
These are valuable for interested researchers from graduate students through established experts.

The book is available from Springer at
[https://link.springer.com/book/9783031621260](https://link.springer.com/book/9783031621260).

This website hosts [code samples](examples/) from the book and provides a [list
of errata]({{ site.baseurl }}/errata/). You are currently viewing the information updated for
OSCAR version {{pv}}.0

## Citation

Please cite the book as described below :

```
[OSCAR-book]
    Wolfram Decker, Christian Eder, Claus Fieker, Max Horn, Michael Joswig, eds.
    The Computer Algebra System OSCAR: Algorithms and Examples,
    Algorithms and Computation in Mathematics, Springer, 2025. (https://link.springer.com/book/9783031621260)
```

If you are using BibTeX, you can use the following BibTeX entry to cite the book :

```bibtex
@book{OSCAR-book,
  editor = {Decker, Wolfram and Eder, Christian and Fieker, Claus and Horn, Max and Joswig, Michael},
  title = {The {C}omputer {A}lgebra {S}ystem {OSCAR}: {A}lgorithms and {E}xamples},
  year = {2025},
  publisher = {Springer},
  series = {Algorithms and {C}omputation in {M}athematics},
  volume = {32},
  edition = {1},
  url = {https://link.springer.com/book/9783031621260},
  issn = {1431-1550},
  doi = {10.1007/978-3-031-62127-7},
}
```
