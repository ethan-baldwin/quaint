# quaint



<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->

<!-- ABOUT THE PROJECT -->
## Overview

This is an R package that uses quartet frequencies from a set of gene trees to characterize introgression across a species tree. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Installation -->
## Installation
   ```r
   library("devtools")
   install_github("ethan-baldwin/quaint")
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Usage -->
## Usage

To use quaint, you need a set of gene trees and a species tree. The tip names in the gene trees should match those in the species tree exactly. The following example assumes you have the gene trees all in one file and a separate file with the species tree, all in newick format.

First, load the trees.
   ```r
   gene_trees <- read.tree("genetrees.tre")
   sp_tree <- read.tree("speciestree.tre")
   ```
Ensure that your species tree is rooted correctly. Optionally assign your outgroup to a variable.
   ```r
   outgroup <- c("your_outgroup")
   sp_tree <- root(sp_tree,outgroup)
   ```
Next, we need to calculate a table with quartet frequencies for every quartet of tips in the tree. For this we the quartetTable function from the package MSCquartets. For trees larger than about 40 tips, this function may take awhile  it is probably worth using the parallel version of this function to make use of multiple cores
   ```r
   # for small datasets, the single core function is fast enough
   quartet_table <- quartetTable(gene_trees,taxonnames = sp_tree$tip.label)
   # for larger datasets, using multiple cores is nice to speed things up
   quartet_table <- quartetTableParallel(gene_trees,taxonnames = sp_tree$tip.label,numCores = 10)
   ```
You can run quaint by either specifying two taxa that you want to test for gene flow between, or by testing for gene flow between all possible pairs of taxa in the tree.
   ```r
   # run quaint on a specific pair of taxa
   test_taxon_names <- c("taxon1","taxon2") # vector containing names of taxa pair you want to test introgression for
   quaint_table <- quaint(test_taxon_names,sp_tree,quartet_table,outgroup)
   # run quaint on all tips
   qt_all_pairs <- quaint_all(sp_tree,quartet_table,outgroup)
   ```
For each pair of taxa, there may be many quartets that can be used to test for gene flow based on the topology of the tree. Quaint conducts a test for every possible quartet and returns a dataframe for each test. It may be informative to look at each test for a single taxon pair, but when running quaint on all taxa you will want to summarize this table, making a dataframe that only has one row for each pair and the summarized results for that pair.
   ```r
   # summarize results by pair, using an adjusted p value cutoff of 0.05 for the chi sq tests
   sum_table <- summarize_quaint_table(qt_all_pairs,alpha = 0.05,use_adjusted_p = TRUE)
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Ethan Baldwin - Ethan.Baldwin@uga.edu

Project Link: [[https://github.com/ethan-baldwin/quaint](https://github.com/ethan-baldwin/quaint/)]

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ethan-baldwin/quaint.svg?style=for-the-badge
[contributors-url]: https://github.com/ethan-baldwin/quaint/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ethan-baldwin/quaint.svg?style=for-the-badge
[forks-url]: https://github.com/ethan-baldwin/quaint/network/members
[stars-shield]: https://img.shields.io/github/stars/ethan-baldwin/quaint.svg?style=for-the-badge
[stars-url]: https://github.com/ethan-baldwin/quaint/stargazers
[issues-shield]: https://img.shields.io/github/issues/ethan-baldwin/quaint.svg?style=for-the-badge
[issues-url]: https://github.com/ethan-baldwin/quaint/issues
[license-shield]: https://img.shields.io/github/license/ethan-baldwin/quaint.svg?style=for-the-badge
[license-url]: https://github.com/ethan-baldwin/quaint/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
