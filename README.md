# Net-Seq analysis tools
These scripts have been developed for parsing the sam files of Net-Seq data, and extracting the information. 
Please see Imashimizu et al. 2015 for more details.

#Requirement
perl<br>

#Usage
"count_error_bp.pl"<br>
*Note: RNA base is complementary to the base detected in the NET-seq read.
<pre><code>$perl count_error_bp.pl -i [sam file] -l [read length] </code></pre>

"parse_samtool_mpileup.pl"<br>
<pre><code>$perl parse_samtool_mpileup.pl -i [pileupfile] -p [0.7] -d [0] </code></pre>

"extractSequence.pl"<br>
<pre><code>$perl extractSequence.pl -i [fasta] -peak [peak position file] -l [10] -comp [false] </code></pre>

#Publication
<a href="http://www.ncbi.nlm.nih.gov/pubmed/25976475" target="_blank">Imashimizu et al. Genome Biol. 2015 May 15;16(1):98.</a>

#Author
Hiroki Takahashi
