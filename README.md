# Google-PageRank

There are two files, sequential.cpp contains the sequential code and parallel.cpp contains the parallel code. "Google PageRank.pdf" contains the file presented during evaluation.
The codes have been tested on three datasets:-
* 25 pages, 300 links
* 500 pages, 124750 links
* 1000 pages, 499500 links

To compile and run sequential file
* `g++ sequential.cpp -o sequential -fopenmp`
* `./sequential`

To compile and run parallel file
* `g++ parallel.cpp -o parallel -fopenmp`
* `./parallel`
