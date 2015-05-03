# Net-Seq analysis tools
These scripts have been developed for parsing the sam files of Net-Seq data.
Please see Imashimizu et al. 2015 for more details.

#Requirement
perl<br>

#Usage
"count_error_bp.pl"<br>
<pre><code>$perl count_error_bp.pl -i [sam file] -l [read length] </code></pre>

"parse_samtool_mpileup.pl"<br>
<pre><code>$perl parse_samtool_mpileup.pl -i [pileupfile] -p [0.7] -d [0] </code></pre>

#Publication
Imashimizu et al. 2015

#Author
Hiroki Takahashi
